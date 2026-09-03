import json

import pytest

from gfviewer.cli import build_parser, main
from gfviewer.errors import InputValidationError
from gfviewer.style import DEFAULT_STYLE, StyleConfig


def test_style_roundtrip_json(tmp_path):
    s = StyleConfig(title="X", dpi=150, tick_style="lollipop", formats=["svg", "png"])
    p = tmp_path / "s.json"
    s.save(str(p))
    loaded = StyleConfig.load(str(p))
    assert loaded.title == "X"
    assert loaded.tick_style == "lollipop"
    assert loaded.formats == ["svg", "png"]


def test_style_roundtrip_yaml(tmp_path):
    yaml = pytest.importorskip("yaml")
    s = StyleConfig(legend_location="outside right", label_mode="family")
    p = tmp_path / "s.yaml"
    s.save(str(p))
    loaded = StyleConfig.load(str(p))
    assert loaded.legend_location == "outside right"
    assert loaded.label_mode == "family"


def test_style_validation_rejects_bad_values():
    with pytest.raises(InputValidationError):
        StyleConfig(legend_location="nowhere").validate()
    with pytest.raises(InputValidationError):
        StyleConfig(formats=["docx"]).validate()
    with pytest.raises(InputValidationError):
        StyleConfig(dpi=5).validate()


def test_from_dict_ignores_unknown_keys():
    cfg = StyleConfig.from_dict({"title": "T", "bogus": 1})
    assert cfg.title == "T"


def test_cli_color_guide(capsys):
    assert main(["--color-guide"]) == 0
    out = capsys.readouterr().out
    assert "red" in out and "palette" in out.lower()


def test_cli_missing_args_returns_2():
    assert main([]) == 2


def test_cli_end_to_end(tmp_path, data_csv, genome_txt):
    out = tmp_path / "out"
    rc = main([
        "-d", data_csv, "-g", genome_txt, "-o", str(out),
        "-f", "svg", "--title", "T", "--tick-style", "lollipop",
        "--analytics", "--permutations", "50",
    ])
    assert rc == 0
    assert (out / "gfviewer.svg").exists()
    assert (out / "analytics_summary.json").exists()


def test_cli_save_style(tmp_path):
    p = tmp_path / "s.json"
    rc = main(["--save-style", str(p), "-d", "x", "-g", "y", "-o", "z",
               "--tick-style", "triangle"])
    assert rc == 0
    assert json.loads(p.read_text())["tick_style"] == "triangle"


def test_show_titles_default_and_no_titles_flag(tmp_path):
    assert DEFAULT_STYLE.show_titles is True
    assert StyleConfig().validate().show_titles is True

    p = tmp_path / "s.json"
    assert main(["--save-style", str(p), "-d", "x", "-g", "y", "-o", "z",
                 "--no-titles"]) == 0
    saved = json.loads(p.read_text())
    assert saved["show_titles"] is False
    assert StyleConfig.load(str(p)).show_titles is False


def test_no_titles_strips_figure_and_chart_titles(tmp_path, data_csv, genome_txt):
    plain = tmp_path / "plain"
    rc = main([
        "-d", data_csv, "-g", genome_txt, "-o", str(plain), "-f", "svg",
        "--title", "My organism", "--analytics", "--permutations", "30",
        "--no-titles",
    ])
    assert rc == 0
    fig_svg = (plain / "gfviewer.svg").read_text()
    assert "gfv-title" not in fig_svg and "My organism" not in fig_svg

    chart = (plain / "analytics_genes_per_family.svg").read_text()
    assert "Genes per family" not in chart

    # and with titles on (default) they are present
    titled = tmp_path / "titled"
    main([
        "-d", data_csv, "-g", genome_txt, "-o", str(titled), "-f", "svg",
        "--title", "My organism", "--analytics", "--permutations", "30",
    ])
    assert "My organism" in (titled / "gfviewer.svg").read_text()
    assert "Genes per family" in (titled / "analytics_genes_per_family.svg").read_text()
