"""GFViewer web portal entry point.

The application is built by :func:`gfviewer_web.create_app`.  For production use
a WSGI server, e.g.::

    gunicorn -w 1 --threads 4 -b 0.0.0.0:5001 "gfviewer_web:create_app()"

(one worker + threads keeps the in-memory job registry coherent; scale with
threads / RENDER workers, not processes).
"""

import os

from gfviewer_web import create_app

app = create_app()

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5001))  # 5000 is taken by AirPlay on macOS
    app.run(host="0.0.0.0", port=port, debug=app.config["DEBUG"])
