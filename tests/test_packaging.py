"""Guards on the pip-packaging metadata so the CLI stays installable."""

import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

try:  # py311+
    import tomllib as _toml
except ModuleNotFoundError:  # py38-310
    try:
        import tomli as _toml
    except ModuleNotFoundError:  # pragma: no cover - only if the parser is absent
        _toml = None


def _pyproject():
    if _toml is None:
        pytest.skip("no TOML parser available")
    with open(os.path.join(ROOT, "pyproject.toml"), "rb") as fh:
        return _toml.load(fh)


def test_console_script_declared_and_importable():
    cfg = _pyproject()
    target = cfg["project"]["scripts"]["gfviewer"]
    assert target == "gfviewer.cli:main"

    module, func = target.split(":")
    mod = __import__(module, fromlist=[func])
    assert callable(getattr(mod, func))


def test_version_matches_package():
    cfg = _pyproject()
    import gfviewer

    assert gfviewer.__version__.count(".") >= 2
    assert cfg["project"]["version"] == gfviewer.__version__, (
        "pyproject.toml version and gfviewer.__version__ have drifted apart"
    )


def test_core_dependencies_have_no_flask():
    cfg = _pyproject()
    core = " ".join(cfg["project"]["dependencies"]).lower()
    assert "flask" not in core and "gunicorn" not in core
    assert "flask" in " ".join(cfg["project"]["optional-dependencies"]["web"]).lower()
    for pkg in ("numpy", "pandas", "matplotlib", "biopython", "pyyaml"):
        assert pkg in core


def test_packaging_side_files_present():
    for name in ("LICENSE", "MANIFEST.in", os.path.join("docs", "INSTALL.md"),
                 os.path.join("scripts", "build_wheel.sh")):
        assert os.path.exists(os.path.join(ROOT, name)), name


def test_cli_entrypoint_runs_as_module():
    out = subprocess.run(
        [sys.executable, "-m", "gfviewer.cli", "--version"],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert out.returncode == 0
    assert "GFViewer" in (out.stdout + out.stderr)
