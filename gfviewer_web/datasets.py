"""Bundled example datasets.

Every dataset lives under ``static/tests/`` (built by ``tests/make_fixtures.py``).
The web app can (a) submit a dataset straight to the job queue and (b) hand the
raw files back as a ZIP.
"""

import io
import os
import zipfile

from werkzeug.datastructures import FileStorage

BASE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "static", "tests")

# ordered -- this is also the display order on the home page
DATASETS = [
    {
        "key": "babesia_mo1_18",
        "title": "Babesia sp. MO1 — 18 gene families",
        "organism": "Babesia sp. MO1",
        "families": 18,
        "format": "Excel table (.xlsx) + colour map",
        "data": ["data_test_1.xlsx"],
        "genome": "chrs_test_1.txt",
        "colors": "colors_test_1.txt",
        "params": {},
        "extra": [],
    },
    {
        "key": "babesia_mo1_6",
        "title": "Babesia sp. MO1 — 6 gene families",
        "organism": "Babesia sp. MO1",
        "families": 6,
        "format": "CSV table + colour map",
        "data": ["data_test_2.csv"],
        "genome": "chrs_test_2.txt",
        "colors": "colors_test_2.txt",
        "params": {},
        "extra": [],
    },
    {
        "key": "babesia_duncani_10",
        "title": "Babesia duncani — 10 gene families (with centromeres)",
        "organism": "Babesia duncani",
        "families": 10,
        "format": "TSV table + genome FASTA",
        "data": ["data_test_3.tsv"],
        "genome": "chrs_test_3.fasta",
        "colors": None,
        "params": {"show_centromeres": True},
        "extra": [],
    },
    {
        "key": "formats",
        "title": "Every input format — BED / GFF3 / GTF",
        "organism": "Babesia sp. MO1 (6 families)",
        "families": 6,
        "format": "BED, per-family BED, BED + mapping, GFF3, GTF",
        "data": ["formats/genes.gff3"],
        "genome": "chrs_test_2.txt",
        "colors": None,
        "params": {},
        # the ZIP carries the whole formats/ folder, not just the file we run
        "extra": [
            "formats/genes.bed", "formats/genes.gtf", "formats/genes_named.bed",
            "formats/mapping.tsv", "formats/README.txt",
            "formats/per_family/MGF1.bed", "formats/per_family/MGF2.bed",
            "formats/per_family/MGF3.bed", "formats/per_family/MGF4.bed",
            "formats/per_family/MGF5.bed", "formats/per_family/MGF6.bed",
        ],
    },
    {
        "key": "arabidopsis_10",
        "title": "Arabidopsis thaliana — 10 synthetic gene families",
        "organism": "Arabidopsis thaliana (TAIR10, 5 chromosomes)",
        "families": 10,
        "format": "synthetic TSV table + genome + colour map",
        "data": ["synthetic/arabidopsis_10/genes.tsv"],
        "genome": "synthetic/arabidopsis_10/genome.txt",
        "colors": "synthetic/arabidopsis_10/colors.txt",
        "params": {},
        "extra": ["synthetic/arabidopsis_10/README.txt"],
    },
    {
        "key": "celegans_20",
        "title": "Caenorhabditis elegans — 20 synthetic gene families",
        "organism": "C. elegans (WBcel235, 6 chromosomes)",
        "families": 20,
        "format": "synthetic TSV table + genome + colour map",
        "data": ["synthetic/celegans_20/genes.tsv"],
        "genome": "synthetic/celegans_20/genome.txt",
        "colors": "synthetic/celegans_20/colors.txt",
        "params": {},
        "extra": ["synthetic/celegans_20/README.txt"],
    },
]

_BY_KEY = {d["key"]: d for d in DATASETS}

# defaults applied to every "run"; a dataset's own ``params`` win
_RUN_DEFAULTS = {
    "with_analytics": True,
    "colocalization": False,
    "formats": ["png"],
    "permutations": "200",
    "on_unknown_chrom": "error",
}


def get(key):
    return _BY_KEY.get(key)


def available(entry):
    """True when every file the dataset needs is on disk."""
    need = list(entry["data"]) + [entry["genome"]]
    if entry.get("colors"):
        need.append(entry["colors"])
    return all(os.path.exists(os.path.join(BASE, p)) for p in need)


def listing():
    """Datasets for the template, each annotated with ``available``."""
    return [dict(d, available=available(d)) for d in DATASETS]


# --------------------------------------------------------------------------- #
def _fs(rel):
    path = os.path.join(BASE, rel)
    return FileStorage(
        stream=open(path, "rb"),
        filename=os.path.basename(path),
        content_type="application/octet-stream",
    )


def open_files(entry):
    """Return ``(data_files, genome_file, extra_files)`` of FileStorage objects.
    The caller must close every stream (``close_files``)."""
    data_files = [_fs(p) for p in entry["data"]]
    genome_file = _fs(entry["genome"])
    extra = {}
    if entry.get("colors"):
        extra["colors"] = _fs(entry["colors"])
    if entry.get("mapping"):
        extra["mapping"] = _fs(entry["mapping"])
    return data_files, genome_file, extra


def close_files(data_files, genome_file, extra):
    for f in list(data_files) + [genome_file] + list(extra.values()):
        try:
            f.close()
        except Exception:  # noqa: BLE001
            pass


def run_params(entry):
    p = dict(_RUN_DEFAULTS)
    p.update(entry.get("params") or {})
    return p


# --------------------------------------------------------------------------- #
def _zip(paths, arc_prefix=""):
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", zipfile.ZIP_DEFLATED) as z:
        for rel in paths:
            src = os.path.join(BASE, rel)
            if os.path.exists(src):
                z.write(src, arcname=os.path.join(arc_prefix, rel) if arc_prefix else rel)
    bio.seek(0)
    return bio


def zip_dataset(entry):
    paths = list(entry["data"]) + [entry["genome"]] + list(entry.get("extra") or [])
    if entry.get("colors"):
        paths.append(entry["colors"])
    if entry.get("mapping"):
        paths.append(entry["mapping"])
    # de-dup, keep order
    seen, uniq = set(), []
    for p in paths:
        if p not in seen:
            seen.add(p)
            uniq.append(p)
    return _zip(uniq, arc_prefix=entry["key"])


def zip_all():
    """Every dataset file in one archive (skips the large legacy render PDFs
    and the old per-test zips)."""
    keep = []
    for root, _dirs, files in os.walk(BASE):
        for name in files:
            if name.endswith((".mgf.pdf", ".zip", ".DS_Store")):
                continue
            keep.append(os.path.relpath(os.path.join(root, name), BASE))
    return _zip(sorted(keep), arc_prefix="gfviewer_datasets")
