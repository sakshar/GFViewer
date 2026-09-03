"""Environment-driven configuration for the GFViewer web app.

Every value can be overridden with an environment variable of the same name
(prefixed ``GFVIEWER_`` for the app-specific ones).
"""

import os


def _int(name, default):
    try:
        return int(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


class Config:
    SECRET_KEY = os.environ.get("SECRET_KEY", "dev-only-change-me")
    DEBUG = os.environ.get("FLASK_DEBUG", "0") == "1"

    # storage
    DATA_DIR = os.environ.get(
        "GFVIEWER_DATA_DIR",
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "instance", "jobs"),
    )
    EXAMPLES_FOLDER = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), "static", "tests"
    )

    # limits
    MAX_UPLOAD_MB = _int("GFVIEWER_MAX_UPLOAD_MB", 25)
    WORKERS = _int("GFVIEWER_WORKERS", 2)
    JOB_TTL_HOURS = _int("GFVIEWER_JOB_TTL_HOURS", 24)
    MAX_JOBS = _int("GFVIEWER_MAX_JOBS", 500)
    RENDER_TIMEOUT_S = _int("GFVIEWER_RENDER_TIMEOUT_S", 120)

    # usage monitor
    # where the aggregate counters live (default: next to DATA_DIR, resolved in
    # create_app so a custom DATA_DIR keeps its own file)
    USAGE_FILE = os.environ.get("GFVIEWER_USAGE_FILE") or None
    # if set, /stats and /api/stats require  ?token=<this>  (or an X-Stats-Token
    # header); leave empty to expose the aggregate counts publicly
    STATS_TOKEN = os.environ.get("GFVIEWER_STATS_TOKEN", "")

    # analytics defaults
    PERMUTATIONS = _int("GFVIEWER_PERMUTATIONS", 500)
