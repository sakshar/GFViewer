"""Exception types used across GFViewer.

``InputValidationError`` carries a user-facing message that is safe to show
directly in the web UI or print to the terminal -- it never contains a
stack trace or internal path noise.
"""


class GFViewerError(Exception):
    """Base class for all expected GFViewer failures."""


class InputValidationError(GFViewerError):
    """Raised when user-supplied files or parameters are invalid.

    Parameters
    ----------
    message:
        Short, plain-language description of what is wrong.
    hints:
        Optional list of concrete suggestions for fixing the problem.
    """

    def __init__(self, message, hints=None):
        super().__init__(message)
        self.message = message
        self.hints = list(hints) if hints else []

    def __str__(self):
        if not self.hints:
            return self.message
        bullets = "\n".join("  - " + h for h in self.hints)
        return "{}\n{}".format(self.message, bullets)
