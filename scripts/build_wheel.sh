#!/usr/bin/env bash
#
# Build the GFViewer command-line tool as a wheel + source distribution,
# validate the artifacts, and print the commands to install or publish them.
#
# Usage:
#   scripts/build_wheel.sh                # build + check into ./dist
#   PYTHON=python3.11 scripts/build_wheel.sh
#   scripts/build_wheel.sh --install      # also pip-install the fresh wheel
#   scripts/build_wheel.sh --testpypi     # also upload to TestPyPI (needs a token)
#   scripts/build_wheel.sh --publish      # also upload to PyPI          (needs a token)
#
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

PYTHON="${PYTHON:-python3}"
do_install=0 ; do_testpypi=0 ; do_publish=0
for arg in "$@"; do
    case "$arg" in
        --install)  do_install=1 ;;
        --testpypi) do_testpypi=1 ;;
        --publish)  do_publish=1 ;;
        -h|--help)  sed -n '3,16p' "$0" ; exit 0 ;;
        *) echo "unknown option: $arg" >&2 ; exit 2 ;;
    esac
done

echo ">> Python:      $("$PYTHON" --version 2>&1)  ($("$PYTHON" -c 'import sys;print(sys.executable)'))"
echo ">> Project:     $("$PYTHON" -c 'import tomllib,pathlib;d=tomllib.loads(pathlib.Path("pyproject.toml").read_text());print(d["project"]["name"])' 2>/dev/null || echo gfviewer)"

echo ">> Cleaning build/ dist/ *.egg-info"
rm -rf build dist ./*.egg-info src/*.egg-info

echo ">> Ensuring build backend + twine are available"
"$PYTHON" -m pip install --quiet --upgrade build twine

echo ">> Building sdist + wheel"
"$PYTHON" -m build

echo ">> Validating artifacts with twine"
"$PYTHON" -m twine check dist/*

echo
echo ">> Built:"
ls -la dist

wheel="$(ls dist/gfviewer-*.whl | head -n1)"

if [ "$do_install" -eq 1 ]; then
    echo
    echo ">> Installing $wheel"
    "$PYTHON" -m pip install --force-reinstall "$wheel"
    "$PYTHON" -m gfviewer.cli --version || gfviewer --version
fi

if [ "$do_testpypi" -eq 1 ]; then
    echo ">> Uploading to TestPyPI"
    "$PYTHON" -m twine upload --repository testpypi dist/*
fi

if [ "$do_publish" -eq 1 ]; then
    echo ">> Uploading to PyPI"
    "$PYTHON" -m twine upload dist/*
fi

cat <<EOF

Next steps
----------
  Install the wheel locally     :  pip install "$wheel"
  Smoke-test the CLI            :  scripts/smoke_cli.sh
  Upload to TestPyPI           :  $PYTHON -m twine upload --repository testpypi dist/*
  Install from TestPyPI       :  pip install --index-url https://test.pypi.org/simple/ \\
                                   --extra-index-url https://pypi.org/simple gfviewer
  Upload to PyPI              :  $PYTHON -m twine upload dist/*
EOF
