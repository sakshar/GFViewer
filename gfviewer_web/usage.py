"""Lightweight, privacy-respecting usage counter for the GFViewer web portal.

Tracks *how many people use the tool* without storing anything that identifies
them:

* **page views** (total and per page),
* **unique visitors per day** -- counted from a salted daily hash of
  ``client-ip + user-agent``; the raw IP / UA are never written anywhere and the
  hash set lives only in memory for the current day,
* **jobs** submitted / completed / failed, **example-dataset runs**, and
  **downloads**,
* a rolling **daily series** (kept ~120 days) for a simple dashboard.

State is a small JSON file (``instance/usage.json`` by default), flushed at most
once every couple of seconds.  Losing it, or a mid-day restart, only costs a
little dedup accuracy -- never correctness of the running totals.
"""

import hashlib
import json
import os
import re
import threading
import time
from datetime import datetime, timezone

# collapse per-resource paths so the "pages" map stays meaningful and bounded
_PATH_NORMALISE = (
    (re.compile(r"^/jobs/[^/]+$"), "/jobs/<id>"),
    (re.compile(r"^/datasets/[^/]+/(run|download)$"), r"/datasets/<key>/\1"),
)

_TOTAL_KEYS = (
    "page_views", "visits",
    "jobs_submitted", "jobs_completed", "jobs_failed",
    "dataset_runs", "downloads",
)
_DAILY_KEYS = ("page_views", "visits", "jobs_submitted", "jobs_completed")
_MAX_DAILY = 120
_MAX_PAGES = 60
_FLUSH_EVERY_S = 2.0
# job lifecycle events emitted by JobManager -> the total they bump
_JOB_EVENTS = {"job_done": "jobs_completed", "job_error": "jobs_failed"}


def _utcnow():
    return datetime.now(timezone.utc).replace(microsecond=0)


class UsageStats:
    def __init__(self, path, secret=""):
        self.path = os.path.abspath(path)
        self._salt = (secret or "gfviewer") + "|usage|v1"
        self._lock = threading.Lock()
        self._seen_date = None          # date string the seen-set belongs to
        self._seen = set()              # visitor hashes seen today (in memory only)
        self._last_flush = 0.0
        self._data = self._load()

    # ------------------------------------------------------------------ #
    def _load(self):
        base = {
            "since": _utcnow().isoformat(),
            "updated": _utcnow().isoformat(),
            "totals": {k: 0 for k in _TOTAL_KEYS},
            "pages": {},
            "daily": {},
        }
        try:
            with open(self.path) as fh:
                disk = json.load(fh)
            if isinstance(disk, dict):
                base["since"] = disk.get("since", base["since"])
                base["totals"].update(
                    {k: int(v) for k, v in (disk.get("totals") or {}).items()
                     if k in _TOTAL_KEYS and isinstance(v, (int, float))}
                )
                base["pages"] = {str(k): int(v) for k, v in (disk.get("pages") or {}).items()
                                 if isinstance(v, (int, float))}
                base["daily"] = {str(k): v for k, v in (disk.get("daily") or {}).items()
                                 if isinstance(v, dict)}
        except (OSError, ValueError):
            pass
        return base

    def _flush(self, force=False):
        now = time.time()
        if not force and now - self._last_flush < _FLUSH_EVERY_S:
            return
        self._last_flush = now
        self._data["updated"] = _utcnow().isoformat()
        tmp = self.path + ".tmp"
        try:
            os.makedirs(os.path.dirname(self.path), exist_ok=True)
            with open(tmp, "w") as fh:
                json.dump(self._data, fh, indent=2)
            os.replace(tmp, self.path)
        except OSError:
            pass

    # ------------------------------------------------------------------ #
    def _today(self):
        return _utcnow().date().isoformat()

    def _day_row(self, day):
        row = self._data["daily"].get(day)
        if row is None:
            row = {k: 0 for k in _DAILY_KEYS}
            self._data["daily"][day] = row
            if len(self._data["daily"]) > _MAX_DAILY:
                for old in sorted(self._data["daily"])[:-_MAX_DAILY]:
                    self._data["daily"].pop(old, None)
        for k in _DAILY_KEYS:
            row.setdefault(k, 0)
        return row

    # ------------------------------------------------------------------ #
    def record_page(self, path):
        """Count one HTML page view of *path* (per-resource ids collapsed)."""
        path = (path or "/")[:120]
        for rx, repl in _PATH_NORMALISE:
            if rx.match(path):
                path = rx.sub(repl, path)
                break
        with self._lock:
            self._data["totals"]["page_views"] += 1
            self._day_row(self._today())["page_views"] += 1
            pages = self._data["pages"]
            pages[path] = pages.get(path, 0) + 1
            if len(pages) > _MAX_PAGES:
                # drop the least-visited path so a crawler can't blow up the map
                least = min(pages, key=pages.get)
                pages.pop(least, None)
            self._flush()

    def record_visit(self, ip, user_agent):
        """Count a unique visitor for today (salted daily hash; no PII stored)."""
        day = self._today()
        raw = "{}|{}|{}|{}".format(day, self._salt, ip or "-", (user_agent or "-")[:200])
        h = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]
        with self._lock:
            if self._seen_date != day:
                self._seen_date = day
                self._seen = set()
            if h in self._seen:
                return
            self._seen.add(h)
            self._data["totals"]["visits"] += 1
            self._day_row(day)["visits"] += 1
            self._flush(force=True)   # one write per visitor/day -- cheap, worth persisting

    def incr(self, key, n=1):
        """Bump a total counter (also mirrors job counts into the daily series)."""
        key = _JOB_EVENTS.get(key, key)
        if key not in _TOTAL_KEYS:
            return
        with self._lock:
            self._data["totals"][key] += n
            if key in _DAILY_KEYS:
                self._day_row(self._today())[key] += n
            self._flush(force=True)   # jobs / downloads are low-frequency, high-value

    def job_event(self, event):
        """Hook target for :class:`JobManager` (``"job_done"`` / ``"job_error"``)."""
        self.incr(event)

    # ------------------------------------------------------------------ #
    def snapshot(self, days=45):
        """A JSON-serialisable copy: totals, top pages, and the recent daily series."""
        with self._lock:
            daily = self._data["daily"]
            recent = sorted(daily)[-days:]
            pages = sorted(self._data["pages"].items(), key=lambda kv: -kv[1])[:20]
            return {
                "since": self._data["since"],
                "updated": self._data["updated"],
                "totals": dict(self._data["totals"]),
                "pages": [{"path": p, "views": v} for p, v in pages],
                "daily": [dict({"date": d}, **{k: daily[d].get(k, 0) for k in _DAILY_KEYS})
                          for d in recent],
            }
