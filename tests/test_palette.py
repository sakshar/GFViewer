import pytest

from gfviewer.errors import InputValidationError
from gfviewer.palette import (
    HARD_ABOVE,
    WARN_ABOVE,
    build_palette,
    distinct_colors,
    parse_color_file,
)


def test_distinct_colors_length_and_range():
    cols = distinct_colors(120)
    assert len(cols) == 120
    for r, g, b in cols:
        assert 0.0 <= r <= 1.0 and 0.0 <= g <= 1.0 and 0.0 <= b <= 1.0
    assert len(set(cols)) > 100  # not degenerate


def test_build_palette_basic():
    fams = ["A", "B", "C"]
    cmap, warn, collapsed = build_palette(fams)
    assert set(cmap) == set(fams)
    assert not warn and not collapsed


def test_warn_above_threshold():
    fams = ["F%d" % i for i in range(WARN_ABOVE + 2)]
    _, warn, collapsed = build_palette(fams)
    assert warn and not collapsed


def test_hard_cap_without_collapse_raises():
    fams = ["F%d" % i for i in range(HARD_ABOVE + 5)]
    with pytest.raises(InputValidationError):
        build_palette(fams)


def test_hard_cap_with_collapse():
    fams = ["F%d" % i for i in range(HARD_ABOVE + 5)]
    counts = {f: (100 - i) for i, f in enumerate(fams)}  # F0 most frequent
    cmap, warn, collapsed = build_palette(
        fams, collapse_rare=True, family_counts=counts
    )
    assert "Other" in cmap
    assert len(cmap) <= HARD_ABOVE
    assert "F0" in cmap and "F0" not in collapsed
    assert len(collapsed) >= 5


def test_keep_families_survive_collapse():
    fams = ["F%d" % i for i in range(HARD_ABOVE + 5)]
    counts = {f: 1 for f in fams}
    rare = "F%d" % (HARD_ABOVE + 4)
    cmap, _, collapsed = build_palette(
        fams, collapse_rare=True, family_counts=counts, keep=[rare]
    )
    assert rare in cmap and rare not in collapsed


def test_color_file_index(tmp_path):
    p = tmp_path / "c.txt"
    p.write_text("A,1\nB,2\n")
    cmap = parse_color_file(str(p))
    assert cmap["A"] == pytest.approx((0.902, 0.098, 0.294), abs=1e-2)


def test_color_file_rgb_and_hex(tmp_path):
    p = tmp_path / "c.txt"
    p.write_text("A,#ff0000\nB,0,1,0\nC,255,255,0\n")
    cmap = parse_color_file(str(p))
    assert cmap["A"] == (1.0, 0.0, 0.0)
    assert cmap["B"] == (0.0, 1.0, 0.0)
    assert cmap["C"] == (1.0, 1.0, 0.0)


def test_color_file_must_match_families(tmp_path):
    p = tmp_path / "c.txt"
    p.write_text("A,1\nB,2\n")
    with pytest.raises(InputValidationError):
        build_palette(["A", "B", "C"], color_file=str(p))
