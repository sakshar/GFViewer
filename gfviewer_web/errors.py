"""Uniform error handling: never leak a traceback to the client."""

from flask import jsonify, render_template, request
from werkzeug.exceptions import HTTPException

from gfviewer.errors import GFViewerError, InputValidationError


def _wants_json():
    if request.path.startswith("/api/"):
        return True
    accept = request.accept_mimetypes
    return accept["application/json"] >= accept["text/html"]


def register_error_handlers(app):
    @app.errorhandler(InputValidationError)
    def _bad_input(exc):
        payload = {"error": exc.message, "hints": exc.hints}
        if _wants_json():
            return jsonify(payload), 400
        return render_template("error.html", error_message=str(exc), hints=exc.hints), 400

    @app.errorhandler(GFViewerError)
    def _gfv_error(exc):
        if _wants_json():
            return jsonify({"error": str(exc)}), 400
        return render_template("error.html", error_message=str(exc), hints=[]), 400

    @app.errorhandler(413)
    def _too_large(exc):
        msg = "Upload exceeds the {} MB limit.".format(app.config["MAX_UPLOAD_MB"])
        if _wants_json():
            return jsonify({"error": msg}), 413
        return render_template("error.html", error_message=msg, hints=[]), 413

    @app.errorhandler(HTTPException)
    def _http(exc):
        if _wants_json():
            return jsonify({"error": exc.description}), exc.code
        return render_template("error.html", error_message=exc.description, hints=[]), exc.code

    @app.errorhandler(Exception)
    def _unhandled(exc):
        app.logger.exception("Unhandled error")
        msg = "An unexpected server error occurred. Please try again."
        if _wants_json():
            return jsonify({"error": msg}), 500
        return render_template("error.html", error_message=msg, hints=[]), 500
