"""Flask application factory for the GFViewer web portal.

The heavy lifting (validation, rendering, analytics) is done **in-process** by
the :mod:`gfviewer` package -- there is no ``subprocess`` and no per-chromosome
GenBank round-trip.  Long renders run on a small thread pool behind an
asynchronous job API so a slow upload never ties up a request worker.
"""

import os

from flask import Flask, request

from gfviewer_web.config import Config
from gfviewer_web.jobs import JobManager
from gfviewer_web.usage import UsageStats


def _client_ip(req):
    """Best-effort client address (only ever hashed, never stored)."""
    fwd = req.headers.get("X-Forwarded-For", "")
    if fwd:
        return fwd.split(",")[0].strip()
    return req.remote_addr or "-"


def create_app(config=None):
    app = Flask(__name__, static_folder="../static", template_folder="../templates")
    app.config.from_object(Config())
    if config:
        app.config.update(config)

    app.config["MAX_CONTENT_LENGTH"] = app.config["MAX_UPLOAD_MB"] * 1024 * 1024

    usage_file = app.config.get("USAGE_FILE") or os.path.join(
        os.path.dirname(os.path.abspath(app.config["DATA_DIR"])), "usage.json"
    )
    usage = UsageStats(usage_file, secret=app.config.get("SECRET_KEY", ""))
    app.extensions["usage"] = usage

    manager = JobManager(
        data_dir=app.config["DATA_DIR"],
        workers=app.config["WORKERS"],
        ttl_hours=app.config["JOB_TTL_HOURS"],
        max_jobs=app.config["MAX_JOBS"],
        event_hook=usage.job_event,
    )
    app.extensions["jobs"] = manager
    manager.start_cleanup()

    @app.after_request
    def _count_visit(response):
        try:
            if (request.method == "GET"
                    and response.status_code == 200
                    and (response.mimetype or "").startswith("text/html")
                    and not request.path.startswith(("/static", "/api"))):
                usage.record_page(request.path)
                usage.record_visit(_client_ip(request),
                                   request.headers.get("User-Agent", ""))
        except Exception:  # noqa: BLE001 - never break a response over stats
            pass
        return response

    from gfviewer_web.routes import bp as main_bp
    from gfviewer_web.api import bp as api_bp

    app.register_blueprint(main_bp)
    app.register_blueprint(api_bp)

    from gfviewer_web.errors import register_error_handlers

    register_error_handlers(app)
    return app
