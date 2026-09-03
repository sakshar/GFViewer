"""Unified feature loader for GFViewer.

``load_features`` accepts one or more annotation files in any mix of these
formats (auto-detected by extension, then by content):

===============  =====================================================
Format           Notes
===============  =====================================================
Table            ``.xlsx`` / ``.csv`` / ``.tsv`` / ``.txt`` with the
                 columns ``gene_id, gene_family, chromosome, start,
                 end, strand`` (header names are matched
                 case-insensitively and a few synonyms are accepted).
                 Coordinates are 1-based inclusive.
BED              ``.bed`` (BED3-BED12).  0-based half-open coordinates
                 are converted to 1-based inclusive.
GFF3             ``.gff`` / ``.gff3``.  Attribute column parsed as
                 ``key=value;`` pairs.
GTF              ``.gtf``.  Attribute column parsed as
                 ``key "value";`` pairs.
===============  =====================================================

Family assignment (the "support all three" BED workflow, generalised to
GFF/GTF):

* **Multiple files** -> each file's rows get ``gene_family`` = the file's
  base name, unless the file itself already carries family information.
* **One file + a mapping file** (``gene_id,gene_family`` per line) -> family
  looked up from the map.
* **One file, nothing else** -> for BED, column 4 (``name``) is the family and
  ``gene_id`` is synthesised; for GFF/GTF the ``--family-attr`` attribute (or a
  sensible default search) supplies the family.

Every loader funnels into one normalised :class:`pandas.DataFrame` with columns::

    gene_id  gene_family  chromosome  start  end  strand  kind

where ``start <= end`` are 1-based inclusive, ``strand`` is one of ``+ - .`` and
``kind`` is ``"gene"`` or ``"centromere"``.
"""

import gzip
import os
import re
from urllib.parse import unquote

import pandas as pd

from gfviewer.errors import InputValidationError
from gfviewer.genome import classify_sequences

NORMALISED_COLUMNS = [
    "gene_id",
    "gene_family",
    "chromosome",
    "start",
    "end",
    "strand",
    "kind",
    "placed",   # True: on a chromosome; False: on an unplaced / stray contig
]

_CENTROMERE_NAMES = {"centromere", "cen", "centromere_region", "centromeric_repeat"}

# header synonyms -> canonical name (matched lower-case, non-alnum stripped)
_HEADER_SYNONYMS = {
    "geneid": "gene_id",
    "id": "gene_id",
    "name": "gene_id",
    "gene": "gene_id",
    "genefamily": "gene_family",
    "family": "gene_family",
    "mgf": "gene_family",
    "chromosome": "chromosome",
    "chrom": "chromosome",
    "chr": "chromosome",
    "seqid": "chromosome",
    "seqname": "chromosome",
    "scaffold": "chromosome",
    "start": "start",
    "chromstart": "start",
    "begin": "start",
    "end": "end",
    "chromend": "end",
    "stop": "end",
    "strand": "strand",
    "orientation": "strand",
    "kind": "kind",
    "type": "kind",
    "feature": "kind",
}

_BED_EXT = {".bed"}
_GFF3_EXT = {".gff", ".gff3"}
_GTF_EXT = {".gtf", ".gff2"}
_TABLE_EXT = {".xlsx", ".xls", ".csv", ".tsv", ".txt"}

_DEFAULT_FAMILY_ATTRS = (
    "gene_family",
    "genefamily",
    "family",
    "Family",
    "gene_biotype",
    "biotype",
    "Name",
    "gene_name",
)
_DEFAULT_ID_ATTRS = ("gene_id", "ID", "Name", "transcript_id", "locus_tag")
_DEFAULT_GFF_TYPES = ("gene", "pseudogene", "mRNA", "transcript")


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def load_features(
    paths,
    genome,
    family_attr=None,
    id_attr=None,
    gff_types=None,
    mapping_file=None,
    on_unknown_chrom="error",
    coord_bounds="clip",
):
    """Load and normalise one or more annotation files.

    Parameters
    ----------
    paths:
        A single path or a list of paths.
    genome:
        ``{seq_id: length}`` mapping from :func:`gfviewer.genome.load_genome`.
    family_attr, id_attr:
        Attribute keys to use for GFF3/GTF family and id.  ``None`` triggers a
        default search (see module docstring).
    gff_types:
        Iterable of GFF/GTF ``type`` values to keep.  ``None`` keeps a default
        set and, if that matches nothing, falls back to every type with a warning.
    mapping_file:
        Optional ``gene_id,gene_family`` lookup table (used only when exactly
        one annotation file is supplied).
    on_unknown_chrom:
        What to do with a feature whose sequence is absent from *genome*:
        ``"error"`` (default), ``"drop"`` it, or ``"unplaced"`` -- keep it and
        mark it as sitting on a stray / unplaced contig (``placed == False``).
        Sequences that *are* in *genome* are classified automatically:
        scaffold/contig/unplaced-looking names (and very short sequences when
        other sequences look like real chromosomes) get ``placed == False``.
    coord_bounds:
        ``"clip"`` (default, warn), ``"error"`` or ``"keep"`` for features that
        extend past the end of their chromosome.

    Returns
    -------
    (DataFrame, list[str])
        The normalised table and a list of non-fatal warning messages.
    """
    if isinstance(paths, (str, bytes, os.PathLike)):
        paths = [paths]
    paths = [os.fspath(p) for p in paths]
    if not paths:
        raise InputValidationError("No annotation file was provided.")
    for p in paths:
        if not os.path.isfile(p):
            raise InputValidationError("Annotation file not found: {!r}".format(p))

    multi = len(paths) > 1
    warnings = []
    frames = []
    for path in paths:
        fmt = _detect_format(path)
        if fmt == "table":
            df = _read_table(path)
            name_col = None
        elif fmt == "bed":
            df, name_col = _read_bed(path)
        elif fmt == "gff3":
            df = _read_gff(path, "gff3", family_attr, id_attr, gff_types, warnings)
            name_col = None
        elif fmt == "gtf":
            df = _read_gff(path, "gtf", family_attr, id_attr, gff_types, warnings)
            name_col = None
        else:  # pragma: no cover - _detect_format never returns anything else
            raise InputValidationError("Unrecognised annotation format: {!r}".format(path))

        if fmt == "bed":
            _assign_bed_family(df, name_col, path, multi, bool(mapping_file))
        else:
            _fill_family_from_basename(df, path, fmt, multi)
        frames.append(df)

    features = pd.concat(frames, ignore_index=True)

    if mapping_file:
        if multi:
            warnings.append(
                "Mapping file ignored because more than one annotation file was given."
            )
        else:
            features = _apply_mapping(features, mapping_file)

    features = _finalise(features, genome, on_unknown_chrom, coord_bounds, warnings)
    return features, warnings


# --------------------------------------------------------------------------- #
# Format detection
# --------------------------------------------------------------------------- #
def _detect_format(path):
    base = path[:-3] if path.endswith(".gz") else path
    ext = os.path.splitext(base)[1].lower()
    if ext in _BED_EXT:
        return "bed"
    if ext in _GFF3_EXT:
        return "gff3"
    if ext in _GTF_EXT:
        return "gtf"
    if ext in (".xlsx", ".xls"):
        return "table"
    if ext in _TABLE_EXT:
        # .csv/.tsv/.txt: decide from the header line
        return _sniff_text_format(path)
    return _sniff_text_format(path)


def _sniff_text_format(path):
    with _open_text(path) as fh:
        first = ""
        for line in fh:
            if line.strip() and not line.startswith("#"):
                first = line.rstrip("\n")
                break
    if first.startswith("##gff-version 3"):
        return "gff3"
    low = first.lower()
    canon = {_canon_header(tok) for tok in re.split(r"[,\t]", low)}
    if "gene_id" in canon or "gene_family" in canon or {"chromosome", "start", "end"} <= canon:
        return "table"
    cols = first.split("\t")
    if len(cols) >= 8 and _looks_int(cols[3]) and _looks_int(cols[4]):
        # 1-based coord columns 4/5 -> GFF/GTF; attribute column decides sub-type
        return "gtf" if '"' in first else "gff3"
    if len(cols) >= 3 and _looks_int(cols[1]) and _looks_int(cols[2]):
        return "bed"
    # last resort: treat as a table and let column validation complain clearly
    return "table"


# --------------------------------------------------------------------------- #
# Table reader (.xlsx / .csv / .tsv / .txt)
# --------------------------------------------------------------------------- #
def _read_table(path):
    base = path[:-3] if path.endswith(".gz") else path
    ext = os.path.splitext(base)[1].lower()
    try:
        if ext in (".xlsx", ".xls"):
            raw = pd.read_excel(path, sheet_name=0, dtype=str)
        elif ext == ".csv":
            raw = pd.read_csv(path, dtype=str, sep=",")
        else:
            raw = pd.read_csv(path, dtype=str, sep=None, engine="python")
    except Exception as exc:  # noqa: BLE001 - surface a clean message
        raise InputValidationError(
            "Could not read {!r} as a table: {}".format(os.path.basename(path), exc)
        )

    rename = {}
    for col in raw.columns:
        canon = _canon_header(str(col))
        target = _HEADER_SYNONYMS.get(canon)
        if target:
            rename.setdefault(target, col)
    missing = [c for c in ("gene_id", "gene_family", "chromosome", "start", "end") if c not in rename]
    if missing:
        raise InputValidationError(
            "{!r} is missing required column(s): {}".format(
                os.path.basename(path), ", ".join(missing)
            ),
            hints=[
                "Expected header: gene_id, gene_family, chromosome, start, end, strand",
                "Found: " + ", ".join(str(c) for c in raw.columns),
            ],
        )

    df = pd.DataFrame(
        {
            "gene_id": raw[rename["gene_id"]],
            "gene_family": raw[rename["gene_family"]],
            "chromosome": raw[rename["chromosome"]],
            "start": raw[rename["start"]],
            "end": raw[rename["end"]],
        }
    )
    df["strand"] = raw[rename["strand"]] if "strand" in rename else "."
    df["kind"] = None
    return df


# --------------------------------------------------------------------------- #
# BED reader
# --------------------------------------------------------------------------- #
def _read_bed(path):
    """Return ``(df, name_col)`` where *name_col* is the BED column-4 values
    (or ``None`` if the file has no 4th column).  ``gene_id`` is always the
    synthesised ``chrom:start-end``; ``gene_family`` is left blank for the
    caller to fill according to the single/multi/mapping workflow."""
    rows = []
    names = []
    with _open_text(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            s = raw.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue
            parts = s.split("\t") if "\t" in s else s.split()
            if len(parts) < 3:
                raise InputValidationError(
                    "BED line {} in {!r} has fewer than 3 columns: {!r}".format(
                        lineno, os.path.basename(path), s
                    )
                )
            chrom = parts[0]
            start0 = _req_int(parts[1], path, lineno, "start")
            end = _req_int(parts[2], path, lineno, "end")
            name = parts[3] if len(parts) >= 4 and parts[3] not in (".", "") else None
            strand = parts[5] if len(parts) >= 6 else "."
            names.append(name)
            rows.append(
                {
                    "gene_id": "{}:{}-{}".format(chrom, start0 + 1, end),
                    "gene_family": None,
                    "chromosome": chrom,
                    "start": start0 + 1,  # 0-based half-open -> 1-based inclusive
                    "end": end,
                    "strand": strand,
                    "kind": None,
                }
            )
    if not rows:
        raise InputValidationError("No features found in BED file {!r}".format(os.path.basename(path)))
    df = pd.DataFrame(rows, columns=NORMALISED_COLUMNS)
    name_col = pd.Series(names, dtype="object") if any(n is not None for n in names) else None
    return df, name_col


def _assign_bed_family(df, name_col, path, multi, has_mapping):
    """Fill ``gene_family`` / ``gene_id`` for a BED file per the workflow:

    * multiple files          -> family = file base name, id = column 4 (or synth)
    * one file + mapping file  -> id = column 4 (or synth), family from mapping later
    * one file, nothing else   -> family = column 4, id stays ``chrom:start-end``
    """
    if multi:
        df["gene_family"] = _basename_family(path)
        if name_col is not None:
            df["gene_id"] = name_col.where(name_col.notna(), df["gene_id"]).values
    elif has_mapping:
        if name_col is not None:
            df["gene_id"] = name_col.where(name_col.notna(), df["gene_id"]).values
    else:
        if name_col is None or name_col.isna().all():
            raise InputValidationError(
                "{!r} has fewer than 4 columns, so no family name is available.".format(
                    os.path.basename(path)
                ),
                hints=[
                    "Add a 4th 'name' column used as the family, or",
                    "supply one BED per family, or a gene_id->family mapping file.",
                ],
            )
        df["gene_family"] = name_col.values


def _fill_family_from_basename(df, path, fmt, multi):
    needs = df["gene_family"].isna() | (df["gene_family"].astype(str).str.strip() == "")
    if not needs.any():
        return
    if multi:
        df.loc[needs, "gene_family"] = _basename_family(path)
    elif fmt in ("gff3", "gtf"):
        raise InputValidationError(
            "Could not determine a gene family for records in {!r}.".format(
                os.path.basename(path)
            ),
            hints=[
                "Pass --family-attr <key> naming the attribute that holds the family, or",
                "supply one file per family, or a gene_id->family mapping file.",
            ],
        )


# --------------------------------------------------------------------------- #
# GFF3 / GTF reader
# --------------------------------------------------------------------------- #
def _read_gff(path, flavour, family_attr, id_attr, gff_types, warnings):
    wanted_types = set(t.lower() for t in gff_types) if gff_types else None
    rows = []
    skipped_types = set()
    with _open_text(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            col = raw.rstrip("\n").split("\t")
            if len(col) < 8:
                raise InputValidationError(
                    "{} line {} in {!r} has {} columns (need at least 8).".format(
                        flavour.upper(), lineno, os.path.basename(path), len(col)
                    )
                )
            seqid, _src, ftype, start, end, _score, strand = col[:7]
            attr_str = col[8] if len(col) > 8 else ""
            attrs = _parse_gtf_attrs(attr_str) if flavour == "gtf" else _parse_gff3_attrs(attr_str)

            ftype_l = ftype.lower()
            is_cen = ftype_l in _CENTROMERE_NAMES
            if not is_cen:
                if wanted_types is not None:
                    if ftype_l not in wanted_types:
                        skipped_types.add(ftype)
                        continue
                elif ftype_l not in _DEFAULT_GFF_TYPES:
                    skipped_types.add(ftype)
                    continue

            fam = _first_attr(attrs, [family_attr] if family_attr else _DEFAULT_FAMILY_ATTRS)
            gid = _first_attr(attrs, [id_attr] if id_attr else _DEFAULT_ID_ATTRS)
            rows.append(
                {
                    "gene_id": gid or "{}:{}-{}".format(seqid, start, end),
                    "gene_family": "centromere" if is_cen else fam,
                    "chromosome": seqid,
                    "start": _req_int(start, path, lineno, "start"),
                    "end": _req_int(end, path, lineno, "end"),
                    "strand": strand,
                    "kind": "centromere" if is_cen else None,
                }
            )

    if not rows and wanted_types is None and skipped_types:
        # Nothing matched the default type set -- retry keeping everything.
        warnings.append(
            "No {} records of type {} found in {!r}; keeping all feature types "
            "instead ({}).".format(
                flavour.upper(),
                "/".join(_DEFAULT_GFF_TYPES),
                os.path.basename(path),
                ", ".join(sorted(skipped_types)),
            )
        )
        return _read_gff(
            path, flavour, family_attr, id_attr, sorted(skipped_types), warnings
        )
    if not rows:
        raise InputValidationError(
            "No usable records found in {!r}.".format(os.path.basename(path))
        )
    if skipped_types:
        warnings.append(
            "{}: ignored feature types {}".format(
                os.path.basename(path), ", ".join(sorted(skipped_types))
            )
        )
    return pd.DataFrame(rows, columns=NORMALISED_COLUMNS)


def _parse_gff3_attrs(text):
    out = {}
    for chunk in text.strip().strip(";").split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        out[key.strip()] = unquote(value.strip())
    return out


_GTF_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _parse_gtf_attrs(text):
    out = {}
    for key, value in _GTF_ATTR_RE.findall(text):
        out[key] = value
    if not out:  # some GTFs use key=value or key value; be forgiving
        for chunk in text.strip().strip(";").split(";"):
            chunk = chunk.strip()
            m = re.match(r"(\w+)\s*[= ]\s*(.+)", chunk)
            if m:
                out[m.group(1)] = m.group(2).strip().strip('"')
    return out


def _first_attr(attrs, keys):
    lower = {k.lower(): v for k, v in attrs.items()}
    for key in keys:
        if key is None:
            continue
        if key in attrs:
            return attrs[key]
        if key.lower() in lower:
            return lower[key.lower()]
    return None


# --------------------------------------------------------------------------- #
# Mapping file
# --------------------------------------------------------------------------- #
def _apply_mapping(df, mapping_file):
    if not os.path.isfile(mapping_file):
        raise InputValidationError("Mapping file not found: {!r}".format(mapping_file))
    mapping = {}
    with _open_text(mapping_file) as fh:
        for lineno, raw in enumerate(fh, start=1):
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            parts = re.split(r"[,\t]", s) if ("," in s or "\t" in s) else s.split()
            if len(parts) < 2:
                raise InputValidationError(
                    "Mapping line {} is not '<gene_id><sep><gene_family>': {!r}".format(lineno, s)
                )
            if lineno == 1 and _canon_header(parts[1]) in ("genefamily", "family"):
                continue  # header row
            mapping[parts[0].strip()] = parts[1].strip()
    if not mapping:
        raise InputValidationError("Mapping file {!r} contained no entries.".format(mapping_file))

    mapped = df["gene_id"].map(mapping)
    unmatched = df.loc[mapped.isna() & (df["kind"] != "centromere"), "gene_id"].unique()
    if len(unmatched):
        raise InputValidationError(
            "{} gene id(s) are missing from the mapping file, e.g. {}".format(
                len(unmatched), ", ".join(map(str, unmatched[:8]))
            )
        )
    df = df.copy()
    df.loc[mapped.notna(), "gene_family"] = mapped[mapped.notna()]
    return df


# --------------------------------------------------------------------------- #
# Normalisation + validation
# --------------------------------------------------------------------------- #
def _finalise(df, genome, on_unknown_chrom, coord_bounds, warnings):
    df = df.copy()
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df["gene_family"] = df["gene_family"].astype(str).str.strip()
    df["chromosome"] = df["chromosome"].astype(str).str.strip()

    # ---- numeric coordinates ------------------------------------------------
    for col in ("start", "end"):
        coerced = pd.to_numeric(df[col], errors="coerce")
        if coerced.isna().any():
            bad = df.loc[coerced.isna(), col].astype(str).unique()[:5]
            raise InputValidationError(
                "Non-numeric {} coordinate(s): {}".format(col, ", ".join(bad))
            )
        df[col] = coerced.astype("int64")

    # ---- strand -----------------------------------------------------------
    df["strand"] = df["strand"].map(_normalise_strand)

    # ---- kind / centromere ---------------------------------------------------
    fam_lower = df["gene_family"].str.lower()
    is_cen = df["kind"].eq("centromere") | fam_lower.isin(_CENTROMERE_NAMES) | (
        (df["strand"] == ".") & fam_lower.isin(_CENTROMERE_NAMES)
    )
    # historical convention: strand token "0" + family "centromere"
    df.loc[is_cen, "kind"] = "centromere"
    df.loc[~is_cen, "kind"] = "gene"
    df.loc[df["kind"] == "centromere", "gene_family"] = "centromere"
    df.loc[df["kind"] == "centromere", "strand"] = "."

    # ---- start <= end -----------------------------------------------------
    swapped = df["start"] > df["end"]
    if swapped.any():
        warnings.append("{} feature(s) had start > end and were swapped.".format(int(swapped.sum())))
        lo = df[["start", "end"]].min(axis=1)
        hi = df[["start", "end"]].max(axis=1)
        df["start"], df["end"] = lo, hi

    if (df["start"] < 1).any():
        n = int((df["start"] < 1).sum())
        warnings.append("{} feature(s) had start < 1 and were clamped to 1.".format(n))
        df.loc[df["start"] < 1, "start"] = 1

    # ---- chromosome membership -------------------------------------------
    known = set(genome)
    unknown = sorted(set(df.loc[~df["chromosome"].isin(known), "chromosome"]))
    if unknown:
        preview = ", ".join(unknown[:10])
        if on_unknown_chrom == "drop":
            warnings.append(
                "Dropped features on {} sequence(s) not in the genome file: {}".format(
                    len(unknown), preview
                )
            )
            df = df[df["chromosome"].isin(known)].copy()
        elif on_unknown_chrom == "unplaced":
            warnings.append(
                "{} sequence(s) not in the genome file are treated as unplaced "
                "contigs: {}".format(len(unknown), preview)
            )
        else:
            raise InputValidationError(
                "Feature file references chromosome(s) absent from the genome "
                "file: {}".format(preview),
                hints=[
                    "Genome file defines: " + ", ".join(list(genome)[:10]),
                    "Make the chromosome identifiers match exactly (case-sensitive), or",
                    "pass on_unknown_chrom='unplaced' to count them as stray contigs.",
                ],
            )

    # ---- chromosome vs. unplaced-contig classification ------------------
    chrom_flags = classify_sequences(genome)
    chrom_set = {s for s, is_chrom in chrom_flags.items() if is_chrom}
    df["placed"] = df["chromosome"].isin(chrom_set)
    n_unplaced = int((~df["placed"]).sum())
    if n_unplaced:
        seqs = sorted(set(df.loc[~df["placed"], "chromosome"]))
        warnings.append(
            "{} gene feature(s) sit on {} sequence(s) treated as unplaced / "
            "stray contigs ({}{}); they are counted separately in the analytics "
            "and are hidden from the figure unless 'show unplaced' is enabled.".format(
                n_unplaced, len(seqs), ", ".join(seqs[:8]),
                ", ..." if len(seqs) > 8 else "",
            )
        )

    # ---- coordinate bounds (placed features only) ----------------------
    lengths = df["chromosome"].map(genome)
    over = df["placed"] & (df["end"] > lengths)
    if over.any():
        n = int(over.sum())
        if coord_bounds == "error":
            ex = df.loc[over].iloc[0]
            raise InputValidationError(
                "Feature {!r} ends at {} but chromosome {!r} is only {} bp.".format(
                    ex["gene_id"], ex["end"], ex["chromosome"], genome[ex["chromosome"]]
                )
            )
        if coord_bounds == "clip":
            warnings.append("{} feature(s) extended past the chromosome end and were clipped.".format(n))
            df.loc[over, "end"] = lengths[over]
            df.loc[df["start"] > df["end"], "start"] = df["end"]

    # ---- empty result -----------------------------------------------------
    if df.empty:
        raise InputValidationError("No features remain after validation.")

    genes = df[df["kind"] == "gene"]
    if genes.empty:
        raise InputValidationError(
            "The annotation contains centromeres but no gene-family features."
        )

    df = df.reset_index(drop=True)
    return df[NORMALISED_COLUMNS]


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _canon_header(text):
    return re.sub(r"[^a-z0-9]", "", str(text).strip().lower())


def _basename_family(path):
    base = os.path.basename(path)
    for ext in (".gz",):
        if base.endswith(ext):
            base = base[: -len(ext)]
    base = os.path.splitext(base)[0]
    return re.sub(r"\s+", "_", base.strip()) or "family"


def _normalise_strand(value):
    s = str(value).strip()
    if s in ("+", "-"):
        return s
    if s in ("1", "+1"):
        return "+"
    if s in ("-1",):
        return "-"
    return "."


def _looks_int(text):
    try:
        int(str(text).replace(",", "").strip())
        return True
    except (ValueError, AttributeError):
        return False


def _req_int(text, path, lineno, what):
    try:
        return int(str(text).replace(",", "").strip())
    except ValueError:
        raise InputValidationError(
            "Line {} of {!r}: {} value {!r} is not an integer.".format(
                lineno, os.path.basename(path), what, text
            )
        )
