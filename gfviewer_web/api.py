"""JSON + file API for render jobs (everything under ``/api``)."""

import os

from flask import Blueprint, current_app, jsonify, request, send_file

from gfviewer.errors import InputValidationError

bp = Blueprint("api", __name__, url_prefix="/api")

_FORM_LIST_FIELDS = {"formats", "gff_types"}
_BOOL_FIELDS = {
    "with_analytics", "collapse_rare", "colocalization",
    "show_centromeres", "show_unplaced", "no_split_strand", "no_legend",
    "show_titles",
}


def _jobs():
    return current_app.extensions["jobs"]


def _usage():
    return current_app.extensions.get("usage")


def _track(event):
    u = _usage()
    if u is not None:
        u.incr(event)


def _stats_authorized():
    token = current_app.config.get("STATS_TOKEN") or ""
    if not token:
        return True
    given = request.args.get("token") or request.headers.get("X-Stats-Token", "")
    return given == token


def _parse_params(form):
    params = {}
    for key in form:
        if key in _FORM_LIST_FIELDS:
            params[key] = [v for v in form.getlist(key) if v]
        else:
            params[key] = form.get(key)
    for b in _BOOL_FIELDS:
        params[b] = form.get(b) in ("1", "true", "on", "yes")
    params.setdefault("formats", [])
    return params


@bp.post("/jobs")
def create_job():
    data_files = request.files.getlist("data")
    genome_file = request.files.get("genome")
    extra = {
        "colors": request.files.get("colors"),
        "mapping": request.files.get("mapping"),
        "style": request.files.get("style"),
    }
    params = _parse_params(request.form)
    job_id = _jobs().submit(data_files, genome_file, params, extra_files=extra)
    _track("jobs_submitted")
    return jsonify({"job_id": job_id, "status_url": "/api/jobs/{}/status".format(job_id)}), 202


@bp.get("/jobs/<job_id>/status")
def job_status(job_id):
    meta = _jobs().get(job_id)
    if not meta:
        return jsonify({"error": "Unknown job id."}), 404
    return jsonify({
        "id": meta["id"],
        "state": meta["state"],
        "error": meta.get("error"),
        "warnings": meta.get("warnings", []),
        "families": meta.get("families", []),
        "chromosomes": meta.get("chromosomes", []),
        "color_map": meta.get("color_map", {}),
        "formats": meta.get("formats", []),
        "has_analytics": meta.get("has_analytics", False),
        "style": meta.get("style", {}),
    })


@bp.get("/jobs/<job_id>/svg")
def job_svg(job_id):
    path = _jobs().svg_path(job_id)
    if not path:
        return jsonify({"error": "No figure available yet."}), 404
    resp = send_file(path, mimetype="image/svg+xml")
    resp.headers["Cache-Control"] = "no-store"
    return resp


@bp.post("/jobs/<job_id>/restyle")
def job_restyle(job_id):
    overrides = request.get_json(silent=True) or {}
    if not isinstance(overrides, dict):
        raise InputValidationError("Restyle payload must be a JSON object.")
    meta = _jobs().restyle(job_id, overrides)
    return jsonify({"ok": True, "style": meta.get("style", {}),
                    "formats": meta.get("formats", [])})


_FIG_FORMATS = {"pdf", "svg", "png", "jpg", "jpeg", "tif", "tiff", "eps"}


@bp.get("/jobs/<job_id>/download/<fmt>")
def job_download(job_id, fmt):
    if fmt.lower() not in _FIG_FORMATS:
        return jsonify({"error": "Unsupported format."}), 400
    path = _jobs().render_format(job_id, fmt)
    _track("downloads")
    return send_file(path, as_attachment=True,
                     download_name="gfviewer.{}".format("jpg" if fmt == "jpeg" else fmt))


_CHART_NAMES = {"genes_per_family", "positional_profile", "ripley", "family_proximity"}


@bp.get("/jobs/<job_id>/chart/<fmt>")
def job_chart(job_id, fmt):
    if fmt.lower() not in _FIG_FORMATS:
        return jsonify({"error": "Unsupported format."}), 400
    which = request.args.get("name", "genes_per_family")
    if which not in _CHART_NAMES:
        return jsonify({"error": "Unknown chart."}), 400
    path = _jobs().render_chart(job_id, fmt, which)
    _track("downloads")
    return send_file(
        path, as_attachment=True,
        download_name="{}.{}".format(which, "jpg" if fmt == "jpeg" else fmt),
    )


@bp.get("/jobs/<job_id>/bundle")
def job_bundle(job_id):
    all_fmt = request.args.get("all", "").lower() in ("1", "true", "yes", "on")
    path = _jobs().bundle(job_id, all_formats=all_fmt)
    _track("downloads")
    return send_file(path, as_attachment=True, mimetype="application/zip",
                     download_name=os.path.basename(path))


@bp.get("/jobs/<job_id>/analytics/<name>")
def job_analytics(job_id, name):
    if "/" in name or "\\" in name or not name.startswith("analytics_"):
        return jsonify({"error": "Bad file name."}), 400
    path = _jobs().analytics_path(job_id, name)
    if not path:
        return jsonify({"error": "Not found."}), 404
    return send_file(path, as_attachment=True)


@bp.get("/health")
def health():
    return jsonify({"status": "ok"})


@bp.get("/stats")
def stats():
    """Aggregate usage counters (no personal data). Gated by STATS_TOKEN if set."""
    if not _stats_authorized():
        return jsonify({"error": "A valid ?token= is required for usage stats."}), 403
    u = _usage()
    if u is None:
        return jsonify({"error": "Usage tracking is not enabled."}), 404
    return jsonify(u.snapshot())
