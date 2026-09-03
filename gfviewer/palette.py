"""Colour-palette generation for gene families.

The first 20 colours are the hand-picked, maximally-distinct set that GFViewer
has always used.  Beyond that the palette is extended programmatically with
perceptually spaced colours (golden-angle hue rotation with alternating
lightness/saturation), so there is no hard 19/20-family ceiling.

A *soft* cap keeps output readable:

* more than ``WARN_ABOVE`` (25) families           -> warning
* more than ``HARD_ABOVE`` (40) families           -> error, unless the caller
  opts into collapsing rare families into a single "Other" bucket.

The centromere colour is fixed (mid grey) and never handed out to a family.
"""

import colorsys

from gfviewer.errors import InputValidationError

WARN_ABOVE = 25
HARD_ABOVE = 40

CENTROMERE_COLOR = (0.50, 0.50, 0.50)
OTHER_LABEL = "Other"
OTHER_COLOR = (0.72, 0.72, 0.72)

# 20 distinct base colours (the classic GFViewer set, grey/black removed since
# grey is reserved for centromeres).
_BASE_HEX = [
    "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
    "#008080", "#dcbeff", "#aa6e28", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000080", "#e6beff",
]

# Named colours accepted in a 2-column colour file (<gf>,<code|name>).
COLOR_NAMES = [
    "red", "green", "yellow", "blue", "orange",
    "purple", "cyan", "magenta", "lime", "pink",
    "teal", "lavender", "brown", "beige", "maroon",
    "mint", "olive", "apricot", "navy", "grape",
]


def _hex_to_rgb(value):
    value = value.lstrip("#")
    return tuple(round(int(value[i : i + 2], 16) / 255.0, 4) for i in (0, 2, 4))


BASE_RGB = [_hex_to_rgb(h) for h in _BASE_HEX]
NAME_TO_RGB = dict(zip(COLOR_NAMES, BASE_RGB))


def _generated_rgb(k):
    """Return the *k*-th generated colour (k >= 0) beyond the base set."""
    hue = (0.5 + k * 0.61803398875) % 1.0  # golden-angle rotation
    lightness = 0.45 + 0.18 * ((k // 2) % 3) / 2.0
    saturation = 0.85 if k % 2 == 0 else 0.6
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return (round(r, 4), round(g, 4), round(b, 4))


def distinct_colors(n):
    """Return *n* RGB triples, each in ``[0, 1]``, as distinct as practical."""
    out = list(BASE_RGB[: max(0, n)])
    k = 0
    while len(out) < n:
        out.append(_generated_rgb(k))
        k += 1
    return out


def color_guide():
    """Return ``[(index, name, rgb), ...]`` for the documented palette."""
    return [(i + 1, name, rgb) for i, (name, rgb) in enumerate(zip(COLOR_NAMES, BASE_RGB))]


# --------------------------------------------------------------------------- #
# colour-file parsing
# --------------------------------------------------------------------------- #
def parse_color_file(path):
    """Parse a colour-map file into ``{family: rgb}``.

    Accepted per-line formats::

        <family>,<1-20>                # palette index (1-based)
        <family>,<colour name>         # e.g. "red"
        <family>,#rrggbb               # hex
        <family>,<r>,<g>,<b>           # floats 0-1 or ints 0-255
    """
    import os
    import re

    color_map = {}
    with open(path, "r") as fh:
        for lineno, raw in enumerate(fh, start=1):
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            tokens = [t.strip() for t in re.split(r"[,\t]", s)]
            fam = tokens[0]
            try:
                if len(tokens) == 2:
                    color_map[fam] = _resolve_single(tokens[1])
                elif len(tokens) == 4:
                    color_map[fam] = _resolve_triple(tokens[1:])
                else:
                    raise ValueError("expected 2 or 4 fields")
            except ValueError as exc:
                raise InputValidationError(
                    "Colour file {!r}, line {}: {} ({!r})".format(
                        os.path.basename(path), lineno, exc, s
                    )
                )
    if not color_map:
        raise InputValidationError("Colour file {!r} contained no entries.".format(path))
    return color_map


def _resolve_single(token):
    token = token.strip()
    if token.startswith("#") and len(token) == 7:
        return _hex_to_rgb(token)
    if token.lower() in NAME_TO_RGB:
        return NAME_TO_RGB[token.lower()]
    idx = int(token)  # raises ValueError -> caught by caller
    if not 1 <= idx <= len(BASE_RGB):
        raise ValueError("palette index must be 1-{}".format(len(BASE_RGB)))
    return BASE_RGB[idx - 1]


def _resolve_triple(parts):
    vals = [float(p) for p in parts]
    if any(v > 1.0 for v in vals):
        vals = [v / 255.0 for v in vals]
    if any(not 0.0 <= v <= 1.0 for v in vals):
        raise ValueError("RGB values must be 0-1 (or 0-255)")
    return tuple(round(v, 4) for v in vals)


# --------------------------------------------------------------------------- #
# top-level builder
# --------------------------------------------------------------------------- #
def build_palette(
    families,
    color_file=None,
    collapse_rare=False,
    family_counts=None,
    keep=None,
):
    """Build the ``{family: rgb}`` mapping used by the renderer.

    Parameters
    ----------
    families:
        Ordered list of family names actually present (centromere excluded).
    color_file:
        Optional path to a colour-map file.  When given it must cover exactly
        the families present.
    collapse_rare:
        If ``True`` and there are more than ``HARD_ABOVE`` families, the least
        frequent ones are merged into an "Other" bucket until the count is
        within the cap.
    family_counts:
        ``{family: n_genes}`` -- required when ``collapse_rare`` is set.
    keep:
        Optional explicit set/list of families to keep out of the "Other"
        bucket regardless of frequency.

    Returns
    -------
    (color_map, warnings, collapsed)
        ``color_map`` maps every *displayed* family (possibly including
        ``"Other"``) to an RGB triple; ``collapsed`` is the set of original
        family names folded into "Other" (empty if none).
    """
    families = list(dict.fromkeys(families))  # de-dupe, keep order
    warnings = []
    collapsed = set()

    if color_file:
        cmap = parse_color_file(color_file)
        missing = [f for f in families if f not in cmap]
        extra = [f for f in cmap if f not in families]
        if missing or extra:
            raise InputValidationError(
                "Colour file does not match the families in the data.",
                hints=(
                    (["Missing colours for: " + ", ".join(missing)] if missing else [])
                    + (["Unused entries: " + ", ".join(extra)] if extra else [])
                ),
            )
        return {f: cmap[f] for f in families}, warnings, collapsed

    n = len(families)
    if n > HARD_ABOVE:
        if not collapse_rare:
            raise InputValidationError(
                "{} gene families exceeds the readable limit of {}.".format(n, HARD_ABOVE),
                hints=[
                    "Enable 'collapse rare families' to merge infrequent families "
                    "into a single 'Other' category, or",
                    "provide a colour file and accept reduced legibility, or",
                    "pre-filter the input to the families of interest.",
                ],
            )
        if not family_counts:
            raise InputValidationError(
                "collapse_rare requires per-family gene counts."
            )
        keep = set(keep or ())
        ranked = sorted(
            families, key=lambda f: (f in keep, family_counts.get(f, 0)), reverse=True
        )
        survivors = ranked[: HARD_ABOVE - 1]  # leave room for the "Other" slot
        survivors = [f for f in families if f in set(survivors)]  # restore order
        collapsed = set(families) - set(survivors)
        warnings.append(
            "Collapsed {} rare families into '{}' (kept the {} most frequent).".format(
                len(collapsed), OTHER_LABEL, len(survivors)
            )
        )
        families = survivors + [OTHER_LABEL]
        n = len(families)
    elif n > WARN_ABOVE:
        warnings.append(
            "{} gene families: the legend and colours may be hard to tell apart "
            "above ~{}.".format(n, WARN_ABOVE)
        )

    colors = distinct_colors(n)
    color_map = {}
    for fam, rgb in zip(families, colors):
        color_map[fam] = OTHER_COLOR if fam == OTHER_LABEL else rgb
    return color_map, warnings, collapsed
