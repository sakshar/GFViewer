"""Read chromosome / scaffold lengths from a variety of common formats.

Supported inputs (detected by extension, then by content):

* FASTA (``.fasta`` / ``.fa`` / ``.fna``) -- lengths measured without loading
  whole sequences into memory.
* FASTA index (``.fai``) -- column 1 = name, column 2 = length.
* UCSC ``chrom.sizes`` / ``.genome`` / ``.sizes`` -- ``<name>\\t<length>``.
* Delimited table (``.txt`` / ``.tsv`` / ``.csv``) -- ``<name><sep><length>``
  per line, separator auto-detected from ``,`` ``\\t`` or whitespace.  A header
  row (non-numeric second field) is skipped.

The public entry point is :func:`load_genome`, which returns an ordered
``dict`` mapping sequence id -> length (``int``), preserving file order.
"""

import gzip
import os
import re

from gfviewer.errors import InputValidationError

# Names that clearly denote an unplaced / unlocalized sequence rather than a
# whole chromosome.  ``unloc`` already covers "unlocalized" / "unlocalised";
# both spellings are kept explicit for readability.
_UNPLACED_RE = re.compile(
    r"(?i)(scaffold|contig|\bscf\b|\bctg\b|unplaced|unloc|unlocalized|unlocalised"
    r"|_random\b|\brandom_|chrun|chr_un|_alt\b|_fix\b|\bpatch\b|debris|mito"
    r"|chrM\b|\bMT\b|plastid|chloroplast|apicoplast)"
)
_CHROMOSOME_RE = re.compile(
    r"(?i)^(chr|chrom|chromosome|linkage_?group|lg)?[\s_-]?"
    r"(2L|2R|3L|3R|[0-9]{1,3}|[IVXLCDM]{1,6}|[XYZW])$"
)

_FASTA_EXT = {".fasta", ".fa", ".fna", ".ffn", ".frn"}
_TABLE_EXT = {".txt", ".tsv", ".csv", ".sizes", ".genome", ".chromsizes"}


def _open_text(path):
    """Open ``path`` for text reading, transparently handling ``.gz``."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _strip_ext(path):
    base = path[:-3] if path.endswith(".gz") else path
    return os.path.splitext(base)


def load_genome(path):
    """Return an ordered ``{seq_id: length}`` mapping for *path*.

    Raises :class:`InputValidationError` with an actionable message if the file
    is missing, empty, or cannot be parsed.
    """
    if not path or not os.path.isfile(path):
        raise InputValidationError(
            "Genome / chromosome-lengths file not found: {!r}".format(path)
        )

    root, ext = _strip_ext(path)
    ext = ext.lower()

    if ext == ".fai":
        lengths = _parse_fai(path)
    elif ext in _FASTA_EXT:
        lengths = _parse_fasta_lengths(path)
    elif ext in _TABLE_EXT:
        lengths = _parse_table(path)
    else:
        # Unknown extension: sniff the first non-blank line.
        lengths = _sniff_and_parse(path)

    if not lengths:
        raise InputValidationError(
            "No chromosome lengths could be read from {!r}.".format(os.path.basename(path)),
            hints=[
                "For a FASTA genome use .fasta/.fa/.fna",
                "For a lengths table use lines of '<chromosome>,<length>' or "
                "'<chromosome><TAB><length>'",
            ],
        )

    bad = [name for name, n in lengths.items() if n <= 0]
    if bad:
        raise InputValidationError(
            "These sequences have a non-positive length: {}".format(", ".join(bad[:10]))
        )
    return lengths


# --------------------------------------------------------------------------- #
# Format-specific parsers
# --------------------------------------------------------------------------- #
def _parse_fasta_lengths(path):
    lengths = {}
    name = None
    count = 0
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = count
                # header up to first whitespace, '>' stripped
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                count = 0
            else:
                count += len(line.strip())
    if name is not None:
        lengths[name] = count
    if "" in lengths:
        raise InputValidationError("FASTA file contains a record with an empty header.")
    return lengths


def _parse_fai(path):
    lengths = {}
    with _open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            lengths[parts[0]] = _to_int(parts[1], path, line)
    return lengths


def _split_row(line):
    line = line.rstrip("\n")
    if "\t" in line:
        return line.split("\t")
    if "," in line:
        return line.split(",")
    return line.split()


def _parse_table(path):
    lengths = {}
    with _open_text(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            parts = _split_row(raw)
            if len(parts) < 2:
                raise InputValidationError(
                    "Line {} of {!r} does not look like '<name><sep><length>': {!r}".format(
                        lineno, os.path.basename(path), raw.strip()
                    )
                )
            name, value = parts[0].strip(), parts[1].strip()
            if lineno == 1 and not _looks_int(value):
                # header row -- skip once
                continue
            lengths[name] = _to_int(value, path, raw)
    return lengths


def _sniff_and_parse(path):
    with _open_text(path) as fh:
        for raw in fh:
            if raw.strip():
                first = raw.strip()
                break
        else:
            return {}
    if first.startswith(">"):
        return _parse_fasta_lengths(path)
    return _parse_table(path)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _looks_int(text):
    try:
        int(text.replace(",", "").replace("_", ""))
        return True
    except ValueError:
        return False


def _to_int(text, path, raw):
    try:
        return int(float(text.replace(",", "").replace("_", "")))
    except ValueError:
        raise InputValidationError(
            "Expected an integer length in {!r} but found {!r}".format(
                os.path.basename(path), raw.strip()
            )
        )


# --------------------------------------------------------------------------- #
# chromosome vs. unplaced-contig classification
# --------------------------------------------------------------------------- #
def classify_sequences(lengths, chromosome_ids=None):
    """Return ``{seq_id: is_chromosome}`` for the sequences in *lengths*.

    A sequence is treated as an **unplaced contig** (``is_chromosome == False``)
    when its name matches a scaffold/contig/unplaced pattern, or -- if some other
    sequences look like real chromosomes -- when its own name does not look like
    a chromosome *and* it is much shorter than the median chromosome.

    Pass *chromosome_ids* (an explicit iterable) to override the heuristic; any
    sequence not listed is then considered unplaced.
    """
    ids = list(lengths)
    if chromosome_ids is not None:
        wanted = set(map(str, chromosome_ids))
        return {s: (s in wanted) for s in ids}

    named_unplaced = {s for s in ids if _UNPLACED_RE.search(str(s))}
    named_chrom = {
        s for s in ids
        if s not in named_unplaced and _CHROMOSOME_RE.match(str(s).strip())
    }

    result = {}
    if named_chrom:
        # some names are unambiguous chromosomes -> everything else that is not a
        # recognisable chromosome name and is small is unplaced
        med = _median([lengths[s] for s in named_chrom]) or 1
        for s in ids:
            if s in named_chrom:
                result[s] = True
            elif s in named_unplaced:
                result[s] = False
            else:
                result[s] = lengths[s] >= 0.2 * med
    else:
        # no obvious chromosome names: only demote the explicitly-flagged ones
        for s in ids:
            result[s] = s not in named_unplaced
    if not any(result.values()):          # never classify everything away
        return {s: True for s in ids}
    return result


def chromosome_id_list(lengths, chromosome_ids=None):
    """The ordered list of sequence ids classified as chromosomes."""
    flags = classify_sequences(lengths, chromosome_ids)
    return [s for s in lengths if flags.get(s, True)]


def _median(values):
    vs = sorted(values)
    n = len(vs)
    if n == 0:
        return 0
    return vs[n // 2] if n % 2 else (vs[n // 2 - 1] + vs[n // 2]) / 2.0
