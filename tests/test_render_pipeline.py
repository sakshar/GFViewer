import os

import pytest

from gfviewer.pipeline import run
from gfviewer.render import figure_to_svg, render
from gfviewer.style import StyleConfig


def _run(data, genome, outdir, **kw):
    return run(data_paths=data, genome_path=genome, outdir=str(outdir), **kw)


def test_pipeline_csv_multi_format(data_csv, genome_txt, tmp_path):
    res = _run(data_csv, genome_txt, tmp_path, formats=["pdf", "svg", "png"])
    for fmt in ("pdf", "svg", "png"):
        assert res.outputs[fmt]
        for p in res.outputs[fmt]:
            assert os.path.getsize(p) > 0


def test_pipeline_xlsx_18_families(data_xlsx, tmp_path):
    g = os.path.join(os.path.dirname(data_xlsx), "chrs_test_1.txt")
    res = _run(data_xlsx, g, tmp_path, formats=["svg"])
    assert len(res.families) == 18


def test_pipeline_tsv_with_centromeres_and_analytics(data_tsv, genome_fasta, tmp_path):
    style = StyleConfig(show_centromeres=True, formats=["svg"])
    res = _run(
        data_tsv, genome_fasta, tmp_path, style=style, with_analytics=True,
        analytics_kwargs={"n_permutations": 50},
    )
    assert res.analytics is not None
    assert any(f.endswith(".json") for f in res.analytics_files)


def test_pipeline_bed(bed_file, genome_txt, tmp_path):
    res = _run(bed_file, genome_txt, tmp_path, formats=["svg"])
    assert set(res.families) == {"MGF1", "MGF2", "MGF3", "MGF4", "MGF5", "MGF6"}


def test_pipeline_gff3_and_gtf(gff3_file, gtf_file, genome_txt, tmp_path):
    for f in (gff3_file, gtf_file):
        res = _run(f, genome_txt, tmp_path / os.path.basename(f), formats=["svg"])
        assert len(res.families) == 6


def test_svg_has_interactive_gids(data_csv, genome_txt, tmp_path):
    from gfviewer.genome import load_genome
    from gfviewer.io import load_features
    from gfviewer.palette import build_palette

    g = load_genome(genome_txt)
    feats, _ = load_features(data_csv, g)
    fams = list(dict.fromkeys(feats[feats.kind == "gene"]["gene_family"]))
    cmap, _, _ = build_palette(fams)
    figs = render(feats, g, StyleConfig(), cmap)
    svg = figure_to_svg(figs[0])
    assert 'id="gfv-legend"' in svg
    assert 'id="gfv-chrom-chr1"' in svg
    assert 'id="gfv-fam-MGF1"' in svg


def test_only_chromosomes_subset(data_csv, genome_txt):
    import matplotlib.pyplot as plt

    from gfviewer.genome import load_genome
    from gfviewer.io import load_features
    from gfviewer.palette import build_palette

    g = load_genome(genome_txt)
    feats, _ = load_features(data_csv, g)
    fams = list(dict.fromkeys(feats[feats.kind == "gene"]["gene_family"]))
    cmap, _, _ = build_palette(fams)

    figs = render(feats, g, StyleConfig(only_chromosomes=["chr1", "chr3"]), cmap)
    svg = figure_to_svg(figs[0])
    assert 'id="gfv-chrom-chr1"' in svg
    assert 'id="gfv-chrom-chr3"' in svg
    assert 'id="gfv-chrom-chr2"' not in svg
    for f in figs:
        plt.close(f)

    # an all-unknown / stale selection must not yield an empty figure
    figs = render(feats, g, StyleConfig(only_chromosomes=["nope"]), cmap)
    svg = figure_to_svg(figs[0])
    assert 'id="gfv-chrom-chr1"' in svg and 'id="gfv-chrom-chr2"' in svg
    for f in figs:
        plt.close(f)


def test_only_families_subset(data_csv, genome_txt):
    import matplotlib.pyplot as plt

    from gfviewer.genome import load_genome
    from gfviewer.io import load_features
    from gfviewer.palette import build_palette

    g = load_genome(genome_txt)
    feats, _ = load_features(data_csv, g)
    fams = list(dict.fromkeys(feats[feats.kind == "gene"]["gene_family"]))
    cmap, _, _ = build_palette(fams)

    style = StyleConfig(only_families=["MGF1", "MGF3"])
    figs = render(feats, g, style, cmap)
    svg = figure_to_svg(figs[0])
    assert 'id="gfv-legend"' in svg
    assert 'id="gfv-fam-MGF1"' in svg and 'id="gfv-fam-MGF3"' in svg
    for gone in ("MGF2", "MGF4", "MGF5", "MGF6"):
        assert 'id="gfv-fam-{}"'.format(gone) not in svg

    # the legend is built straight from color_map, which render() has filtered
    from gfviewer.render import _legend_handles
    handles = _legend_handles(style, {k: v for k, v in cmap.items()
                                      if k in style.only_families})
    assert [h.get_label() for h in handles] == ["MGF1", "MGF3"]
    for f in figs:
        plt.close(f)

    # an all-unknown / stale selection must not blank the figure
    figs = render(feats, g, StyleConfig(only_families=["nope"]), cmap)
    svg = figure_to_svg(figs[0])
    assert 'id="gfv-fam-MGF1"' in svg and 'id="gfv-fam-MGF2"' in svg
    for f in figs:
        plt.close(f)


def test_both_orientations_render(data_tsv, genome_fasta, tmp_path):
    for orient in ("horizontal", "vertical"):
        style = StyleConfig(orientation=orient, show_centromeres=True, formats=["svg", "png"])
        res = _run(data_tsv, genome_fasta, tmp_path / orient, style=style)
        for fmt in ("svg", "png"):
            for p in res.outputs[fmt]:
                assert os.path.getsize(p) > 0


def test_vertical_svg_keeps_gids(data_csv, genome_txt, tmp_path):
    from gfviewer.genome import load_genome
    from gfviewer.io import load_features
    from gfviewer.palette import build_palette

    g = load_genome(genome_txt)
    feats, _ = load_features(data_csv, g)
    fams = list(dict.fromkeys(feats[feats.kind == "gene"]["gene_family"]))
    cmap, _, _ = build_palette(fams)
    figs = render(feats, g, StyleConfig(orientation="vertical", show_centromeres=True), cmap)
    svg = figure_to_svg(figs[0])
    assert 'id="gfv-fam-MGF1"' in svg and 'id="gfv-legend"' in svg


def test_labels_widen_chromosome_spacing(data_csv, genome_txt):
    import matplotlib.pyplot as plt

    from gfviewer.genome import load_genome
    from gfviewer.io import load_features
    from gfviewer.palette import build_palette

    g = load_genome(genome_txt)
    feats, _ = load_features(data_csv, g)
    fams = list(dict.fromkeys(feats[feats.kind == "gene"]["gene_family"]))
    cmap, _, _ = build_palette(fams)

    # horizontal: labelled page is taller (more room between stacked chromosomes)
    plain = render(feats, g, StyleConfig(label_mode="none"), cmap)[0]
    lab = render(feats, g, StyleConfig(label_mode="gene_id"), cmap)[0]
    assert lab.get_size_inches()[1] > plain.get_size_inches()[1] + 0.3

    # vertical: labelled page is wider (more room between side-by-side chromosomes).
    # legend off so the page width tracks the content and isn't pinned to the
    # minimum canvas a bottom legend needs.
    vp = render(feats, g, StyleConfig(orientation="vertical", label_mode="none",
                                      legend_show=False), cmap)[0]
    vl = render(feats, g, StyleConfig(orientation="vertical", label_mode="family",
                                      legend_show=False), cmap)[0]
    assert vl.get_size_inches()[0] > vp.get_size_inches()[0] + 0.2

    for f in (plain, lab, vp, vl):
        plt.close(f)


def test_vertical_gene_id_labels_stay_in_view(data_xlsx):
    """18 families, gene_id labels, vertical: every label must sit inside the
    canvas and no chromosome's dense cluster may push labels off the page."""
    import matplotlib.pyplot as plt

    from gfviewer.genome import load_genome
    from gfviewer.io import load_features
    from gfviewer.palette import build_palette

    g = load_genome(os.path.join(os.path.dirname(data_xlsx), "chrs_test_1.txt"))
    feats, _ = load_features(data_xlsx, g)
    fams = list(dict.fromkeys(feats[feats.kind == "gene"]["gene_family"]))
    cmap, _, _ = build_palette(fams)

    figs = render(feats, g, StyleConfig(orientation="vertical", label_mode="gene_id"), cmap)
    for fig in figs:
        fig.canvas.draw()
        ax = fig.axes[0]
        inv = ax.transData.inverted()
        w, h = ax.get_xlim()[1], ax.get_ylim()[1]
        placed = [t for t in ax.texts if (t.get_gid() or "").startswith("gfv-label-")]
        assert placed
        for t in placed:
            bb = t.get_window_extent()
            dx0, dy0 = inv.transform((bb.x0, bb.y0))
            dx1, dy1 = inv.transform((bb.x1, bb.y1))
            assert dx0 >= -0.2 and dx1 <= w + 0.2, (t.get_text(), dx0, dx1, w)
            assert dy0 >= -0.2 and dy1 <= h + 0.2, (t.get_text(), dy0, dy1, h)
        plt.close(fig)


def test_synthetic_datasets_render_and_analyse(synthetic_dataset, tmp_path):
    data, genome = synthetic_dataset
    res = _run(data, genome, tmp_path, formats=["svg"], with_analytics=True,
              analytics_kwargs={"n_permutations": 40})
    assert res.outputs["svg"] and os.path.getsize(res.outputs["svg"][0]) > 0
    assert res.analytics is not None
    assert not res.analytics.family_proximity.empty
    assert not res.analytics.duplication_modes.empty
    assert len(res.families) >= 10


def test_collapse_rare_end_to_end(genome_txt, tmp_path):
    import pandas as pd
    from gfviewer.genome import load_genome

    g = load_genome(genome_txt)
    rows = []
    for fam in range(45):
        for k in range(2 + (fam % 3)):
            s = 1000 + fam * 5000 + k * 300
            rows.append(("g%d_%d" % (fam, k), "F%02d" % fam, "chr1", s, s + 200, "+"))
    csv = tmp_path / "many.csv"
    pd.DataFrame(
        rows, columns=["gene_id", "gene_family", "chromosome", "start", "end", "strand"]
    ).to_csv(csv, index=False)

    res = _run(str(csv), genome_txt, tmp_path / "out", formats=["svg"], collapse_rare=True)
    assert "Other" in res.color_map
    assert len(res.color_map) <= 40
