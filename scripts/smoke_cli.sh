#!/usr/bin/env bash
#
# Post-install sanity check for the `gfviewer` command-line tool.
# Renders one bundled dataset (if present) and confirms the outputs exist.
#
# Usage:
#   scripts/smoke_cli.sh                 # uses the `gfviewer` entry point on PATH
#   GFVIEWER="python -m gfviewer.cli" scripts/smoke_cli.sh
#
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GFVIEWER="${GFVIEWER:-gfviewer}"
export MPLBACKEND="${MPLBACKEND:-Agg}"

echo ">> which: $(command -v gfviewer || true)"
$GFVIEWER --version
$GFVIEWER --color-guide >/dev/null && echo ">> --color-guide ok"

data="$repo_root/static/tests/data_test_2.csv"
genome="$repo_root/static/tests/chrs_test_2.txt"
if [ ! -f "$data" ]; then
    data="$repo_root/static/tests/synthetic/arabidopsis_10/genes.tsv"
    genome="$repo_root/static/tests/synthetic/arabidopsis_10/genome.txt"
fi
if [ ! -f "$data" ]; then
    echo ">> no bundled dataset found (run: python tests/make_fixtures.py) -- skipping render check"
    exit 0
fi

out="$(mktemp -d)"
trap 'rm -rf "$out"' EXIT
echo ">> rendering $data -> $out"
$GFVIEWER -d "$data" -g "$genome" -o "$out" -f pdf svg png \
          --analytics --permutations 50 -cen --title "GFViewer smoke test"

for f in gfviewer.pdf gfviewer.svg analytics_summary.json; do
    if [ -s "$out/$f" ]; then echo "   ok   $f"; else echo "   MISSING $f"; exit 1; fi
done
echo ">> smoke test passed"
