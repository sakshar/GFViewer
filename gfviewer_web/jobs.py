"""Asynchronous render-job manager.

A job owns a directory under ``DATA_DIR/<job_id>/`` containing:

    uploads/            raw user files
    out/                rendered figures + analytics CSV/JSON
    features.pkl        cached normalised feature table (for fast restyle)
    genome.json         cached {chrom: length}
    meta.json           job state, params, warnings, families, style

Rendering runs on a bounded ``ThreadPoolExecutor``.  State is kept in memory and
mirrored to ``meta.json`` so completed jobs survive a process restart.
"""

import json
import os
import shutil
import threading
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import pandas as pd

from gfviewer import analytics as _analytics
from gfviewer import render as _render
from gfviewer.errors import GFViewerError, InputValidationError
from gfviewer.genome import load_genome
from gfviewer.io import load_features
from gfviewer.palette import build_palette
from gfviewer.style import StyleConfig

BASE_FORMATS = ("svg", "pdf")          # always produced on completion
ALLOWED_EXT = {
    "xlsx", "xls", "csv", "tsv", "txt", "bed", "gff", "gff3", "gtf",
    "fasta", "fa", "fna", "fai", "sizes", "genome", "json", "yaml", "yml",
    "gz",
}


def _utcnow():
    return datetime.now(timezone.utc).replace(microsecond=0)


def _allowed(filename):
    if "." not in filename:
        return False
    ext = filename.rsplit(".", 1)[1].lower()
    if ext == "gz" and filename.count(".") >= 2:
        ext = filename.rsplit(".", 2)[1].lower()
    return ext in ALLOWED_EXT


class JobError(GFViewerError):
    pass


class JobManager:
    def __init__(self, data_dir, workers=2, ttl_hours=24, max_jobs=500,
                 event_hook=None):
        self.data_dir = os.path.abspath(data_dir)
        os.makedirs(self.data_dir, exist_ok=True)
        self.ttl_seconds = ttl_hours * 3600
        self.max_jobs = max_jobs
        # optional callable(event_name) for usage tracking -- "job_done" /
        # "job_error"; must never raise into the worker.
        self._event_hook = event_hook
        self._pool = ThreadPoolExecutor(max_workers=max(1, workers))
        self._lock = threading.Lock()
        self._jobs = {}          # job_id -> meta dict
        self._cleanup_thread = None
        self._load_existing()

    def _event(self, name):
        if self._event_hook is None:
            return
        try:
            self._event_hook(name)
        except Exception:  # noqa: BLE001 - tracking must not break rendering
            pass

    # ------------------------------------------------------------------ #
    # paths
    # ------------------------------------------------------------------ #
    def job_dir(self, job_id):
        return os.path.join(self.data_dir, job_id)

    def _p(self, job_id, *parts):
        return os.path.join(self.job_dir(job_id), *parts)

    # ------------------------------------------------------------------ #
    # submission
    # ------------------------------------------------------------------ #
    def submit(self, data_files, genome_file, params, extra_files=None):
        """Persist uploads and queue a render.  Returns the new job id.

        *data_files* is a list of Werkzeug ``FileStorage`` (>=1); *genome_file*
        one ``FileStorage``; *extra_files* an optional dict with keys
        ``colors`` / ``mapping`` / ``style``.
        """
        from werkzeug.utils import secure_filename

        if not data_files or not any(f and f.filename for f in data_files):
            raise InputValidationError("Please upload at least one annotation file.")
        if not genome_file or not genome_file.filename:
            raise InputValidationError("Please upload a genome / chromosome-lengths file.")

        job_id = uuid.uuid4().hex
        up = self._p(job_id, "uploads")
        os.makedirs(up, exist_ok=True)
        os.makedirs(self._p(job_id, "out"), exist_ok=True)

        data_paths = []
        for f in data_files:
            if not f or not f.filename:
                continue
            name = secure_filename(f.filename)
            if not _allowed(name):
                raise InputValidationError("Unsupported file type: {}".format(f.filename))
            dest = os.path.join(up, name)
            f.save(dest)
            data_paths.append(dest)

        gname = secure_filename(genome_file.filename)
        if not _allowed(gname):
            raise InputValidationError("Unsupported genome file type: {}".format(genome_file.filename))
        genome_path = os.path.join(up, gname)
        genome_file.save(genome_path)

        extra_paths = {}
        for key, f in (extra_files or {}).items():
            if f and f.filename:
                nm = secure_filename(f.filename)
                if not _allowed(nm):
                    raise InputValidationError("Unsupported {} file type: {}".format(key, f.filename))
                p = os.path.join(up, nm)
                f.save(p)
                extra_paths[key] = p

        meta = {
            "id": job_id,
            "state": "queued",
            "created": _utcnow().isoformat(),
            "updated": _utcnow().isoformat(),
            "params": params,
            "data_paths": data_paths,
            "genome_path": genome_path,
            "extra_paths": extra_paths,
            "warnings": [],
            "error": None,
            "families": [],
            "formats": [],
            "has_analytics": bool(params.get("with_analytics")),
        }
        with self._lock:
            self._jobs[job_id] = meta
        self._write_meta(job_id)
        self._pool.submit(self._safe_process, job_id)
        self._enforce_max_jobs()
        return job_id

    # ------------------------------------------------------------------ #
    # processing
    # ------------------------------------------------------------------ #
    def _safe_process(self, job_id):
        try:
            self._process(job_id)
        except InputValidationError as exc:
            self._fail(job_id, str(exc))
        except Exception as exc:  # noqa: BLE001
            traceback.print_exc()
            self._fail(job_id, "Rendering failed: {}".format(exc))

    def _process(self, job_id):
        meta = self._jobs[job_id]
        self._set_state(job_id, "running")
        p = meta["params"]

        style = self._style_from_params(p, meta["extra_paths"].get("style"))
        style.formats = list(dict.fromkeys(list(BASE_FORMATS) + p.get("formats", [])))

        genome = load_genome(meta["genome_path"])
        from gfviewer.genome import chromosome_id_list

        plottable_chroms = chromosome_id_list(genome) or list(genome)

        features, warnings = load_features(
            meta["data_paths"], genome,
            family_attr=p.get("family_attr") or None,
            id_attr=p.get("id_attr") or None,
            gff_types=p.get("gff_types") or None,
            mapping_file=meta["extra_paths"].get("mapping"),
            on_unknown_chrom=p.get("on_unknown_chrom", "error"),
            coord_bounds=p.get("coord_bounds", "clip"),
        )

        genes = features[features["kind"] == "gene"]
        fam_order = list(dict.fromkeys(genes["gene_family"].tolist()))
        fam_counts = genes["gene_family"].value_counts().to_dict()
        color_map, pal_warn, collapsed = build_palette(
            fam_order, color_file=meta["extra_paths"].get("colors"),
            collapse_rare=p.get("collapse_rare", False), family_counts=fam_counts,
        )
        warnings += pal_warn
        if collapsed:
            features = features.copy()
            features.loc[features["gene_family"].isin(collapsed), "gene_family"] = "Other"

        features.to_pickle(self._p(job_id, "features.pkl"))
        with open(self._p(job_id, "genome.json"), "w") as fh:
            json.dump(genome, fh)

        out = self._p(job_id, "out")
        figs = _render.render(features, genome, style, color_map)
        outputs = _render.save_figures(figs, out, "gfviewer", style.formats, dpi=style.dpi)
        _render.close_all(figs)

        if meta["has_analytics"]:
            akw = {"n_permutations": int(p.get("permutations", 500)),
                   "colocalization": bool(p.get("colocalization"))}
            res = _analytics.compute_analytics(features, genome, **akw)
            res.write(out, chart_formats=list(style.formats),
                      color_map=color_map, dpi=style.dpi,
                      chart_titles=style.show_titles)

        with self._lock:
            meta["state"] = "done"
            meta["updated"] = _utcnow().isoformat()
            meta["warnings"] = warnings
            meta["families"] = list(color_map)
            meta["color_map"] = {k: list(v) for k, v in color_map.items()}
            meta["formats"] = sorted(outputs)
            meta["chromosomes"] = plottable_chroms
            meta["style"] = style.to_dict()
        self._write_meta(job_id)
        self._event("job_done")

    def _style_from_params(self, p, style_file):
        style = StyleConfig.load(style_file) if style_file else StyleConfig()
        m = {
            "orientation": ("orientation", str),
            "telomere_length": ("telomere_length", int),
            "chromosomes_per_row": ("columns", int),
            "title": ("title", str),
            "subtitle": ("subtitle", str),
            "dpi": ("dpi", int),
            "tick_style": ("tick_style", str),
            "label_mode": ("label_mode", str),
            "legend_location": ("legend_location", str),
            "legend_columns": ("legend_columns", int),
            "row_height_cm": ("row_height_cm", float),
            "length_cm": ("length_cm", float),
            "page_width_cm": ("page_width_cm", float),
        }
        for attr, (key, cast) in m.items():
            if p.get(key) not in (None, ""):
                try:
                    setattr(style, attr, cast(p[key]))
                except (TypeError, ValueError):
                    pass
        if "per_page" in p and p["per_page"]:
            try:
                cols = style.chromosomes_per_row or 1
                style.rows_per_page = max(1, int(round(int(p["per_page"]) / cols)))
            except (TypeError, ValueError):
                pass
        style.show_centromeres = bool(p.get("show_centromeres"))
        style.show_unplaced = bool(p.get("show_unplaced"))
        style.show_titles = bool(p.get("show_titles", True))
        style.split_by_strand = not p.get("no_split_strand")
        if p.get("no_legend"):
            style.legend_show = False
        return style.validate()

    # ------------------------------------------------------------------ #
    # restyle (no re-upload)
    # ------------------------------------------------------------------ #
    def restyle(self, job_id, overrides):
        meta = self.get(job_id)
        if not meta or meta["state"] != "done":
            raise JobError("Job is not ready to restyle.")
        fpkl = self._p(job_id, "features.pkl")
        gjson = self._p(job_id, "genome.json")
        if not (os.path.exists(fpkl) and os.path.exists(gjson)):
            raise JobError("Cached job data is no longer available; please re-run.")

        features = pd.read_pickle(fpkl)
        with open(gjson) as fh:
            genome = json.load(fh)

        base = dict(meta.get("style") or {})
        base.update({k: v for k, v in (overrides or {}).items()})
        style = StyleConfig.from_dict(base).validate()
        style.formats = list(BASE_FORMATS)

        genes = features[features["kind"] == "gene"]
        fam_order = list(dict.fromkeys(genes["gene_family"].tolist()))
        cmap = {k: tuple(v) for k, v in meta.get("color_map", {}).items()}
        # honour colour overrides sent by the editor
        for fam, rgb in (overrides or {}).get("_colors", {}).items():
            if fam in cmap:
                cmap[fam] = tuple(rgb)
        if set(cmap) != set(fam_order) | ({"Other"} if "Other" in cmap else set()):
            cmap, _, _ = build_palette(fam_order)

        out = self._p(job_id, "out")
        _chart_img_ext = (".png", ".jpg", ".jpeg", ".svg", ".pdf", ".eps",
                          ".tif", ".tiff")
        for f in os.listdir(out):
            if f.startswith("gfviewer.") and not f.endswith((".json",)):
                os.remove(os.path.join(out, f))
            # drop cached analytics-chart images so they re-render with the
            # current style (e.g. the show-titles toggle); keep the data files
            elif f.startswith("analytics_") and f.endswith(_chart_img_ext):
                os.remove(os.path.join(out, f))
        figs = _render.render(features, genome, style, cmap)
        _render.save_figures(figs, out, "gfviewer", style.formats, dpi=style.dpi)
        _render.close_all(figs)

        formats = sorted(
            x.split(".")[-1] for x in os.listdir(out) if x.startswith("gfviewer.")
        )
        with self._lock:
            live = self._jobs.get(job_id)          # self.get() returns a copy --
            if live is None:                       # mutate the real record here
                raise JobError("Job disappeared during restyle.")
            live["style"] = style.to_dict()
            live["color_map"] = {k: list(v) for k, v in cmap.items()}
            live["formats"] = formats
            live["updated"] = _utcnow().isoformat()
        self._write_meta(job_id)
        return self.get(job_id)

    def render_format(self, job_id, fmt):
        """Render *fmt* on demand from cached data; return the file path."""
        meta = self.get(job_id)
        if not meta or meta["state"] != "done":
            raise JobError("Job is not ready.")
        fmt = fmt.lower()
        out = self._p(job_id, "out")
        existing = [
            os.path.join(out, f) for f in os.listdir(out)
            if f.startswith("gfviewer.") and f.endswith("." + fmt)
        ]
        if existing:
            return sorted(existing)[0]

        features = pd.read_pickle(self._p(job_id, "features.pkl"))
        with open(self._p(job_id, "genome.json")) as fh:
            genome = json.load(fh)
        style = StyleConfig.from_dict(meta.get("style") or {})
        style.formats = [fmt]
        cmap = {k: tuple(v) for k, v in meta.get("color_map", {}).items()}
        figs = _render.render(features, genome, style, cmap)
        outputs = _render.save_figures(figs, out, "gfviewer", [fmt], dpi=style.dpi)
        _render.close_all(figs)
        return outputs[fmt][0]

    # analytics charts that can be rendered on demand:  which -> csv source
    CHART_SPECS = {
        "genes_per_family": None,                         # from summary json
        "positional_profile": "analytics_positional_profile.csv",
        "ripley": "analytics_ripley.csv",
        "family_proximity": "analytics_family_proximity.csv",
    }

    def render_chart(self, job_id, fmt, which="genes_per_family"):
        """Render an analytics chart in *fmt* on demand; cache it in ``out/``."""
        meta = self.get(job_id)
        if not meta or meta["state"] != "done":
            raise JobError("Job is not ready.")
        if which not in self.CHART_SPECS:
            raise JobError("Unknown analytics chart: {}".format(which))
        fmt = fmt.lower()
        out = self._p(job_id, "out")
        ext = "jpg" if fmt == "jpeg" else fmt
        cached = os.path.join(out, "analytics_{}.{}".format(which, ext))
        if os.path.exists(cached):
            return cached

        from gfviewer import charts as _charts

        cmap = {k: tuple(v) for k, v in meta.get("color_map", {}).items()}
        style = StyleConfig.from_dict(meta.get("style") or {})

        titles = getattr(style, "show_titles", True)

        if which == "genes_per_family":
            spath = os.path.join(out, "analytics_summary.json")
            if not os.path.exists(spath):
                raise JobError("This job has no analytics; re-run with analytics enabled.")
            with open(spath) as fh:
                rows = json.load(fh).get("genes_per_family") or []
            if not rows:
                raise JobError("No per-family gene counts are available.")
            return _charts.write_genes_per_family(
                rows, out, formats=[fmt], color_map=cmap, dpi=style.dpi, titles=titles
            )[0]

        csv_path = os.path.join(out, self.CHART_SPECS[which])
        if not os.path.exists(csv_path):
            raise JobError("This analytic was not produced for this job.")
        if which == "positional_profile":
            df = pd.read_csv(csv_path)
            return _charts.write_positional_profile(
                df, out, formats=[fmt], color_map=cmap, dpi=style.dpi, titles=titles
            )[0]
        if which == "ripley":
            df = pd.read_csv(csv_path)
            return _charts.write_ripley(
                df, out, formats=[fmt], color_map=cmap, dpi=style.dpi, titles=titles
            )[0]
        # family_proximity -- matrix + clustering order from the summary
        mat = pd.read_csv(csv_path, index_col=0)
        summ = {}
        spath = os.path.join(out, "analytics_summary.json")
        if os.path.exists(spath):
            with open(spath) as fh:
                summ = json.load(fh)
        return _charts.write_family_proximity(
            mat, out, order=summ.get("family_proximity_order"),
            clusters=summ.get("family_proximity_clusters"),
            formats=[fmt], dpi=style.dpi, titles=titles,
        )[0]

    # ------------------------------------------------------------------ #
    _BUNDLE_FIG_FORMATS = ("pdf", "svg", "png", "jpg", "eps")

    def bundle(self, job_id, all_formats=False):
        """Zip the job's ``out/`` directory.  With *all_formats*, first render the
        chromosome figure and every analytics chart in every figure format."""
        meta = self.get(job_id)
        if not meta or meta["state"] != "done":
            raise JobError("Job is not ready.")
        out = self._p(job_id, "out")
        if all_formats:
            for f in self._BUNDLE_FIG_FORMATS:
                try:
                    self.render_format(job_id, f)
                except Exception:  # noqa: BLE001
                    pass
                for which in self.CHART_SPECS:
                    try:
                        self.render_chart(job_id, f, which)
                    except Exception:  # noqa: BLE001
                        pass
        zname = "gfviewer_all_formats.zip" if all_formats else "gfviewer_bundle.zip"
        zpath = os.path.join(out, zname)
        try:
            if os.path.exists(zpath):
                os.remove(zpath)
        except OSError:
            pass
        import zipfile

        with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
            for name in sorted(os.listdir(out)):
                if name.endswith(".zip"):
                    continue
                z.write(os.path.join(out, name), arcname=name)
        return zpath

    # ------------------------------------------------------------------ #
    # accessors
    # ------------------------------------------------------------------ #
    def get(self, job_id):
        with self._lock:
            meta = self._jobs.get(job_id)
            return dict(meta) if meta else None

    def svg_path(self, job_id):
        out = self._p(job_id, "out")
        cand = os.path.join(out, "gfviewer.svg")
        if os.path.exists(cand):
            return cand
        pages = sorted(
            os.path.join(out, f) for f in os.listdir(out)
            if f.startswith("gfviewer.p") and f.endswith(".svg")
        )
        return pages[0] if pages else None

    def analytics_path(self, job_id, name):
        out = self._p(job_id, "out")
        cand = os.path.join(out, name)
        return cand if os.path.exists(cand) and os.path.abspath(cand).startswith(out) else None

    # ------------------------------------------------------------------ #
    # internals
    # ------------------------------------------------------------------ #
    def _set_state(self, job_id, state):
        with self._lock:
            self._jobs[job_id]["state"] = state
            self._jobs[job_id]["updated"] = _utcnow().isoformat()
        self._write_meta(job_id)

    def _fail(self, job_id, message):
        with self._lock:
            if job_id in self._jobs:
                self._jobs[job_id]["state"] = "error"
                self._jobs[job_id]["error"] = message
                self._jobs[job_id]["updated"] = _utcnow().isoformat()
        self._write_meta(job_id)
        self._event("job_error")

    def _write_meta(self, job_id):
        meta = self.get(job_id)
        if not meta:
            return
        try:
            with open(self._p(job_id, "meta.json"), "w") as fh:
                json.dump(meta, fh, indent=2, default=str)
        except OSError:
            pass

    def _load_existing(self):
        if not os.path.isdir(self.data_dir):
            return
        for job_id in os.listdir(self.data_dir):
            mp = os.path.join(self.data_dir, job_id, "meta.json")
            if not os.path.exists(mp):
                continue
            try:
                with open(mp) as fh:
                    meta = json.load(fh)
                if meta.get("state") == "running":
                    meta["state"] = "error"
                    meta["error"] = "Interrupted by a server restart."
                self._jobs[job_id] = meta
            except (OSError, ValueError):
                continue

    def _enforce_max_jobs(self):
        with self._lock:
            if len(self._jobs) <= self.max_jobs:
                return
            ordered = sorted(self._jobs.items(), key=lambda kv: kv[1].get("created", ""))
            for job_id, _ in ordered[: len(self._jobs) - self.max_jobs]:
                self._jobs.pop(job_id, None)
                shutil.rmtree(self.job_dir(job_id), ignore_errors=True)

    # ------------------------------------------------------------------ #
    # cleanup thread
    # ------------------------------------------------------------------ #
    def start_cleanup(self, interval=1800):
        if self._cleanup_thread and self._cleanup_thread.is_alive():
            return

        def loop():
            while True:
                try:
                    self.cleanup()
                except Exception:  # noqa: BLE001
                    traceback.print_exc()
                time.sleep(interval)

        self._cleanup_thread = threading.Thread(target=loop, daemon=True)
        self._cleanup_thread.start()

    def cleanup(self):
        now = time.time()
        removed = 0
        for job_id in list(self._jobs):
            d = self.job_dir(job_id)
            try:
                age = now - os.path.getmtime(d)
            except OSError:
                age = self.ttl_seconds + 1
            if age > self.ttl_seconds:
                shutil.rmtree(d, ignore_errors=True)
                with self._lock:
                    self._jobs.pop(job_id, None)
                removed += 1
        # orphan directories with no in-memory entry
        for name in os.listdir(self.data_dir):
            d = os.path.join(self.data_dir, name)
            if os.path.isdir(d) and name not in self._jobs:
                try:
                    if now - os.path.getmtime(d) > self.ttl_seconds:
                        shutil.rmtree(d, ignore_errors=True)
                        removed += 1
                except OSError:
                    pass
        return removed
