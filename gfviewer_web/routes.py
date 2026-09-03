"""HTML pages: home / result / help / about / install, plus dataset endpoints."""

import os

from flask import (
    Blueprint, abort, current_app, redirect, render_template, request, send_file,
    url_for,
)

from gfviewer_web import datasets as _datasets

bp = Blueprint("main", __name__)


def _bump(event):
    u = current_app.extensions.get("usage")
    if u is not None:
        u.incr(event)


_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_INSTALL_MD = os.path.join(_REPO_ROOT, "docs", "INSTALL.md")

# Where users can get the tool.  Kept here so the template and the About page
# stay in sync.
DOWNLOAD_LINKS = {
    "pypi": "https://pypi.org/project/gfviewer/",
    "github": "https://github.com/sakshar/GFViewer",
    "releases": "https://github.com/sakshar/GFViewer/releases",
    "guide": "https://github.com/sakshar/GFViewer/blob/main/docs/INSTALL.md",
    "conda": "https://anaconda.org/bioconda/gfviewer",
}


@bp.get("/")
def index():
    return render_template("index.html", datasets=_datasets.listing(),
                           links=DOWNLOAD_LINKS)


@bp.get("/jobs/<job_id>")
def result(job_id):
    meta = current_app.extensions["jobs"].get(job_id)
    if not meta:
        abort(404, description="Unknown job id.")
    return render_template("result.html", job=meta)


@bp.get("/help")
def help_page():
    return render_template("help.html", links=DOWNLOAD_LINKS)


@bp.get("/about_us")
def about_us():
    return render_template("about_us.html", links=DOWNLOAD_LINKS)


@bp.get("/install")
def install():
    """pip-installation guideline and download links for the CLI tool."""
    return render_template("install.html", links=DOWNLOAD_LINKS,
                           datasets_zip=url_for("main.datasets_download_all"))


@bp.get("/install/guide.md")
def install_guide_md():
    """The raw Markdown installation guide (docs/INSTALL.md)."""
    if not os.path.exists(_INSTALL_MD):
        abort(404, description="Installation guide not found in this deployment.")
    return send_file(_INSTALL_MD, mimetype="text/markdown", as_attachment=True,
                     download_name="GFViewer-INSTALL.md")


@bp.get("/stats")
def stats_page():
    """Human-readable usage dashboard (aggregate counts only, no personal data)."""
    token = current_app.config.get("STATS_TOKEN") or ""
    if token and request.args.get("token") != token:
        abort(403, description="This deployment protects its usage stats. "
              "Append ?token=… to view them.")
    usage = current_app.extensions.get("usage")
    data = usage.snapshot() if usage else None
    return render_template("stats.html", stats=data,
                           api_url=url_for("api.stats"))


# --------------------------------------------------------------------------- #
# example datasets
# --------------------------------------------------------------------------- #
@bp.get("/datasets/<key>/run")
def dataset_run(key):
    entry = _datasets.get(key)
    if not entry:
        abort(404, description="Unknown dataset.")
    if not _datasets.available(entry):
        abort(404, description="This dataset is not installed on the server. "
              "Run  python tests/make_fixtures.py .")
    data_files, genome_file, extra = _datasets.open_files(entry)
    try:
        job_id = current_app.extensions["jobs"].submit(
            data_files, genome_file, _datasets.run_params(entry), extra_files=extra
        )
    finally:
        _datasets.close_files(data_files, genome_file, extra)
    _bump("dataset_runs")
    _bump("jobs_submitted")
    return redirect(url_for("main.result", job_id=job_id))


@bp.get("/datasets/<key>/download")
def dataset_download(key):
    entry = _datasets.get(key)
    if not entry:
        abort(404, description="Unknown dataset.")
    bio = _datasets.zip_dataset(entry)
    _bump("downloads")
    return send_file(bio, mimetype="application/zip", as_attachment=True,
                     download_name="{}.zip".format(key))


@bp.get("/datasets/download-all")
def datasets_download_all():
    bio = _datasets.zip_all()
    _bump("downloads")
    return send_file(bio, mimetype="application/zip", as_attachment=True,
                     download_name="gfviewer_datasets.zip")
