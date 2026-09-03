# GFViewer — pip installation & command-line guide

GFViewer ships a single console command, **`gfviewer`**, that turns gene
coordinates + a genome into publication-quality chromosome ideograms and a
battery of localization statistics. This document covers installing that command
with `pip` and driving it from the terminal.

- [1. Requirements](#1-requirements)
- [2. Install](#2-install)
  - [2a. From PyPI (recommended)](#2a-from-pypi-recommended)
  - [2b. From GitHub](#2b-from-github)
  - [2c. From a source checkout (editable)](#2c-from-a-source-checkout-editable)
  - [2d. With conda / mamba](#2d-with-conda--mamba)
  - [Isolated environments](#isolated-environments)
- [3. Verify the installation](#3-verify-the-installation)
- [4. Upgrade / uninstall](#4-upgrade--uninstall)
- [5. Running `gfviewer` from the command line](#5-running-gfviewer-from-the-command-line)
  - [Synopsis](#synopsis)
  - [Input files](#input-files)
  - [Worked examples](#worked-examples)
  - [Full option reference](#full-option-reference)
  - [What gets written](#what-gets-written)
  - [Exit codes](#exit-codes)
- [6. Bundled example datasets](#6-bundled-example-datasets)
- [7. Running the web portal](#7-running-the-web-portal)
- [8. Troubleshooting](#8-troubleshooting)

---

## 1. Requirements

| | |
|---|---|
| **Python** | 3.8 – 3.12 (CPython) |
| **pip** | 21.3 or newer (`python -m pip install --upgrade pip`) |
| **OS** | Linux, macOS, Windows |
| **Compiler** | none — every dependency installs from a pre-built wheel |
| **Disk** | ~150 MB including NumPy / pandas / matplotlib |

Installing `gfviewer` pulls in `biopython`, `matplotlib`, `pandas`, `numpy`,
`openpyxl`, `Pillow`, `reportlab`, `PyPDF2` and `PyYAML`. The web portal
(`Flask`, `gunicorn`) is an **optional** extra and is not needed for the
command line.

---

## 2. Install

### 2a. From PyPI (recommended)

```bash
python -m pip install --upgrade pip
python -m pip install gfviewer
```

This installs the library and puts the **`gfviewer`** executable on your `PATH`.

### 2b. From GitHub

Install the latest `main`, or pin a tag / commit:

```bash
# latest development version
python -m pip install "git+https://github.com/sakshar/GFViewer.git"

# a specific release tag
python -m pip install "git+https://github.com/sakshar/GFViewer.git@v2.0.0"
```

Or grab a release tarball from
<https://github.com/sakshar/GFViewer/releases> and:

```bash
python -m pip install gfviewer-2.0.0.tar.gz
```

### 2c. From a source checkout (editable)

For development, or to run the bundled tests and web portal:

```bash
git clone https://github.com/sakshar/GFViewer.git
cd GFViewer

python -m pip install -e .            # CLI only, editable
python -m pip install -e ".[web]"     # + Flask web portal
python -m pip install -e ".[dev]"     # + pytest
python -m pip install -e ".[web,dev,build]"   # everything
```

`-e` (editable) means changes to the source are picked up without reinstalling.

### 2d. With conda / mamba

The dependency stack is easiest to get from `conda-forge`:

```bash
conda create -n gfviewer -c conda-forge python=3.11 \
    biopython matplotlib pandas numpy openpyxl pillow reportlab pypdf2 pyyaml
conda activate gfviewer
python -m pip install gfviewer          # or: pip install -e .  from a checkout
```

A ready-made spec is in [`environment.yml`](../environment.yml):

```bash
conda env create -f environment.yml
conda activate gfviewer
python -m pip install -e .
```

### Isolated environments

Installing into a throwaway environment keeps GFViewer's dependencies away from
your system Python.

**venv (standard library):**

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install gfviewer
```

**pipx (one command, isolated automatically):**

```bash
pipx install gfviewer
```

---

## 3. Verify the installation

```bash
gfviewer --version           # -> GFViewer 2.0.0
gfviewer --help              # full option list
gfviewer --color-guide       # the built-in 20-colour palette
```

If `gfviewer` is not found but the install succeeded, your Python scripts
directory is not on `PATH` — call it through the interpreter instead:

```bash
python -m gfviewer.cli --version
```

or add the directory `python -m site --user-base` reports (its `bin` /
`Scripts` sub-folder) to `PATH`.

**End-to-end check** (from a source checkout, after building the example data):

```bash
python tests/make_fixtures.py        # once, builds static/tests/**
scripts/smoke_cli.sh                 # renders one dataset and checks the outputs
```

---

## 4. Upgrade / uninstall

```bash
python -m pip install --upgrade gfviewer
python -m pip show gfviewer          # version, location, dependencies
python -m pip uninstall gfviewer
```

---

## 5. Running `gfviewer` from the command line

### Synopsis

```
gfviewer -d DATA [DATA ...] -g GENOME -o OUTDIR [options]
```

`-d/--data`, `-g/--genome` and `-o/--output` are required for a render.
`--color-guide` and `--save-style` are the only actions that do not need them.

### Input files

**Annotation (`-d`)** — one or more of, detected by extension:

| Type | Extensions | Coordinates | Family comes from |
|---|---|---|---|
| Table | `.xlsx` `.xls` `.csv` `.tsv` `.txt` | 1-based, inclusive | `gene_family` column |
| BED | `.bed` (`.gz` ok) | 0-based, half-open | column 4 — *or* the file's base name (several files) — *or* `-m` mapping |
| GFF3 | `.gff` `.gff3` | 1-based, inclusive | `--family-attr` (default search: `gene_family`, `Family`, `gene_biotype`, `Name`, …) |
| GTF | `.gtf` | 1-based, inclusive | as GFF3 |

Table columns (header synonyms accepted, case-insensitive):
`gene_id, gene_family, chromosome, start, end, strand`. A row whose
`gene_family` is `centromere` and `strand` is `0` marks a centromere.

**Genome (`-g`)** — a FASTA (`.fasta/.fa/.fna`), a samtools `.fai`, a UCSC
`chrom.sizes`, or a plain `chromosome,length` table (one per line). Chromosome
names must match the annotation exactly.

**Colour map (`-c`, optional)** — one family per line:

```
MGF1,1              # palette index 1-20
MGF2,#1b9e77        # hex
MGF3,0.2,0.6,0.9    # RGB in 0-1 (or 0-255)
```

### Worked examples

```bash
# 1. Minimal render (PDF + SVG are always produced)
gfviewer -d genes.csv -g genome.fasta -o out/

# 2. Extra raster/vector formats and a title
gfviewer -d genes.csv -g chrom.sizes -o out/ -f pdf svg png jpg eps \
         --title "Arabidopsis thaliana"

# 3. Draw centromeres, vertical layout, lollipop marks, label each cluster
gfviewer -d genes.tsv -g genome.fai -o out/ \
         -cen --orientation vertical --tick-style lollipop --label-mode family

# 4. Full analytics with a fixed permutation seed
gfviewer -d genes.csv -g genome.fasta -o out/ \
         --analytics --permutations 1000 --seed 42 -cen

# 5. Analytics only, no figure
gfviewer -d genes.csv -g genome.fasta -o out/ --analytics-only

# 6. Several BED files -> one family per file
gfviewer -d families/*.bed -g genome.fasta -o out/

# 7. One BED + an explicit gene->family map
gfviewer -d all_genes.bed -m gene_to_family.tsv -g genome.fasta -o out/

# 8. GFF3 with the family in a custom attribute
gfviewer -d annotation.gff3 -g genome.fasta -o out/ --family-attr locus_type

# 9. Many families: fold rare ones into "Other" above the 40-family cap
gfviewer -d big.csv -g genome.fasta -o out/ --collapse-rare --keep MGF1 MGF7

# 10. Save a style, reuse it later
gfviewer -d genes.csv -g genome.fasta -o out/ --tick-style triangle \
         --font "DejaVu Sans" --legend-location "outside right" \
         --save-style mystyle.yaml
gfviewer -d other.csv -g other.fasta -o out2/ --style mystyle.yaml
```

### Full option reference

Run `gfviewer --help` for the authoritative list. Grouped summary:

**Actions**

| Flag | Effect |
|---|---|
| `--version` | print version and exit |
| `--color-guide` | print the palette and exit |
| `--save-style FILE` | write the fully-resolved style to YAML/JSON and exit |

**Input / output**

| Flag | Meaning |
|---|---|
| `-d, --data FILE [FILE ...]` | annotation file(s) |
| `-g, --genome FILE` | genome FASTA / `.fai` / `chrom.sizes` / `chrom,length` table |
| `-o, --output DIR` | output directory (created if needed) |
| `-c, --colors FILE` | colour-map file |
| `-m, --mapping FILE` | `gene_id,gene_family` map (single annotation file only) |
| `--family-attr KEY` | GFF3/GTF attribute holding the family |
| `--id-attr KEY` | GFF3/GTF attribute holding the gene id |
| `--gff-types T [T ...]` | GFF3/GTF feature types to keep (default: gene-like) |
| `--basename NAME` | output base name (default `gfviewer`) |
| `-f, --format FMT [FMT ...]` | export formats: `pdf svg png jpg tif eps` |
| `--on-unknown-chrom {error,drop,unplaced}` | genes on sequences missing from the genome file |
| `--coord-bounds {clip,error,keep}` | genes whose coordinates fall outside the chromosome |

**Gene families**

| Flag | Meaning |
|---|---|
| `--collapse-rare` | merge infrequent families into `Other` above the cap |
| `--keep FAMILY [FAMILY ...]` | families to protect from the `Other` bucket |

**Analytics**

| Flag | Default | Meaning |
|---|---|---|
| `--analytics` | off | also compute statistics, write `analytics_*.csv/.json/.bed` + 4 figures |
| `--analytics-only` | off | analytics without rendering the ideogram |
| `--subtelomere-fraction F` | `0.10` | fraction of arm length counted as sub-telomeric |
| `--subtelomere-bp N` | — | absolute sub-telomere size (overrides the fraction) |
| `--cluster-gap N` | `50000` | max bp between paralogues in a tandem array |
| `--proximal-window N` | `5 × cluster-gap` | max bp for the `proximal` duplication class |
| `--ripley-scales BP [BP ...]` | `1kb … 1Mb` | scales for the 1-D Ripley K/L test |
| `--hotspot-window N` | `100000` | window (bp) for the multigene-family hotspot scan |
| `--hotspot-step N` | `window / 2` | step (bp) for the hotspot scan |
| `--proximity-clusters K` | auto | clusters to cut the family-proximity tree into |
| `--permutations N` | `1000` | permutations per null model |
| `--seed N` | `0` | RNG seed (reproducible nulls) |
| `--colocalization` | off | pairwise family co-localization test (slower) |

**Style** (each overrides the `--style` file and the defaults)

| Flag | Meaning |
|---|---|
| `--style FILE` | load a YAML/JSON style file |
| `--orientation {horizontal,vertical}` | chromosome direction |
| `-t, --telomere N` | telomere cap length (bp) |
| `-p, --per-page N` | chromosomes per page (rows × columns) |
| `--columns N` | chromosomes side by side (`0` = auto) |
| `--row-height CM` | space per chromosome across its thickness |
| `--length-cm CM` | length-axis budget (vertical) |
| `--single-page` | put every chromosome on one page |
| `--show-unplaced` | also draw unplaced / stray contigs |
| `--title` / `--subtitle` | figure captions |
| `--no-titles` | write the chromosome figure and the analytics-chart images with no title text (for a manuscript that captions them itself) |
| `--dpi N` | raster resolution |
| `--page-size W_CM H_CM` | fixed page size |
| `--background COLOR` | page background |
| `--tick-style {line,lollipop,triangle,box,arrow}` | gene mark shape |
| `--no-split-strand` | draw all genes on one side of the axis |
| `-cen, --centromeres` | draw centromeres |
| `--label-mode {none,family,gene_id}` | gene labels |
| `--label-size PT` | label font size |
| `--font NAME` | font family for labels / legend / chromosome names |
| `-l, --legend-location LOC` | `upper/lower/center` × `left/center/right`, `outside right`, `outside bottom`, `none` |
| `--legend-columns N` | legend columns |
| `--legend-size PT` | legend font size |
| `--legend-title TEXT` | legend heading |
| `--legend-frame` | draw a box around the legend |
| `--legend-separate-page` | legend on its own page |
| `--no-legend` | omit the legend |

### What gets written

Into `OUTDIR/` (base name from `--basename`, default `gfviewer`):

- `gfviewer.pdf`, `gfviewer.svg`, and any extra `-f` formats — the chromosome
  ideogram. Multi-page output adds `.p1`, `.p2`, … suffixes for raster/SVG;
  PDF stays a single multi-page file.
- With `--analytics`: `analytics_*.csv`, `analytics_hotspots.bed`,
  `analytics_summary.json`, and the figures `analytics_genes_per_family.*`,
  `analytics_positional_profile.*`, `analytics_ripley.*`,
  `analytics_family_proximity.*` (one per requested figure format).

The column-by-column layout of every analytics file is documented on the
**Help** page of the web portal (the *Downloadable outputs* section) and in
[`README.md`](../README.md).

`gfviewer` prints each file it writes and any non-fatal warnings to stderr.

### Exit codes

| Code | Meaning |
|---|---|
| `0` | success |
| `1` | a GFViewer error (bad input, unreadable file, validation failure) — message on stderr |
| `2` | usage error (missing `-d/-g/-o`, unknown flag) |

---

## 6. Bundled example datasets

The repository carries ready-to-run datasets under `static/tests/`. Build them
with:

```bash
python tests/make_fixtures.py
```

This produces:

- the three original *Babesia* sets (`data_test_1.xlsx` … `data_test_3.tsv`
  with matching `chrs_test_*`),
- `formats/` — the same 6-family set expressed as BED, per-family BED,
  BED + `mapping.tsv`, GFF3 and GTF,
- `synthetic/arabidopsis_10/` — 10 random gene families on the *A. thaliana*
  TAIR10 chromosomes (`genes.tsv` + `genome.txt`),
- `synthetic/celegans_20/` — 20 random gene families on *C. elegans* WBcel235.

Run one:

```bash
gfviewer -d static/tests/synthetic/arabidopsis_10/genes.tsv \
         -g static/tests/synthetic/arabidopsis_10/genome.txt \
         -o out/arabidopsis --analytics -cen
```

You can also **download any dataset (or all of them) from the web portal's home
page**, or from the *Install* page, without cloning the repository.

---

## 7. Running the web portal

The Flask portal is optional and is run from a source checkout (its templates
and static assets are not part of the wheel).

```bash
python -m pip install -e ".[web]"

# development
python flaskapp.py                       # http://localhost:5001

# production
gunicorn -w 1 --threads 4 -b 0.0.0.0:5001 "gfviewer_web:create_app()"
```

Configuration is via environment variables: `GFVIEWER_DATA_DIR`,
`GFVIEWER_MAX_UPLOAD_MB` (25), `GFVIEWER_WORKERS` (2),
`GFVIEWER_JOB_TTL_HOURS` (24), `GFVIEWER_USAGE_FILE`, `GFVIEWER_STATS_TOKEN`,
`SECRET_KEY`, `PORT`.

A privacy-respecting **usage monitor** is served at `/stats` (HTML) and
`/api/stats` (JSON): page views, unique visitors per day (salted daily hash — no
IPs or cookies stored), job counts, dataset runs and downloads. Set
`GFVIEWER_STATS_TOKEN` to require `?token=…` on those endpoints; counters live in
`GFVIEWER_USAGE_FILE` (default `instance/usage.json`).

---

## 8. Troubleshooting

| Symptom | Fix |
|---|---|
| `gfviewer: command not found` | The scripts dir is not on `PATH`. Use `python -m gfviewer.cli …`, or add the folder from `python -m site --user-base` (its `bin`/`Scripts`). `pipx install gfviewer` avoids this. |
| `ModuleNotFoundError` right after install | You installed into a different interpreter. Run `python -m pip show gfviewer` and compare its `Location` with `python -c "import sys;print(sys.executable)"`. Use a venv. |
| Matplotlib complains about a display / `TclError` on a headless server | Force the non-interactive backend: `export MPLBACKEND=Agg` (GFViewer already selects Agg internally, but a stray `matplotlibrc` can override it). |
| `pip` tries to build NumPy/pandas from source (very old pip) | `python -m pip install --upgrade pip setuptools wheel` first. |
| Fonts differ from the paper | Pass `--font "DejaVu Sans"` (bundled with matplotlib) for a portable result. |
| `error: missing required argument(s): --data, --genome, --output` | Exit code 2 — supply `-d`, `-g` and `-o` (or use `--color-guide` / `--save-style`). |
| Excel file won't read | `openpyxl` handles `.xlsx`; legacy `.xls` needs `xlrd` (`pip install xlrd`). Prefer CSV/TSV. |
| Build a wheel yourself | `scripts/build_wheel.sh` (needs `pip install build twine`). |

For anything else, open an issue at
<https://github.com/sakshar/GFViewer/issues> with the command you ran and the
full error.
