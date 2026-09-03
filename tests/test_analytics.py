import numpy as np
import pandas as pd
import pytest

from gfviewer.analytics import _bh_fdr, _binom_test, _poisson_sf, compute_analytics


@pytest.fixture
def rich():
    """3 families with distinct signatures + a centromere track."""
    genome = {"chr1": 2_000_000, "chr2": 2_000_000}
    rows = []
    for i in range(9):                       # TAND: tight tandem array, chr1 start
        s = 1000 + i * 1500
        rows.append(("t%d" % i, "TAND", "chr1", s, s + 400, "+"))
    for i in range(8):                       # SPREAD: evenly spread over chr2
        s = 100_000 + i * 230_000
        rows.append(("s%d" % i, "SPREAD", "chr2", s, s + 400, "-" if i % 2 else "+"))
    for i in range(6):                       # PERI: hugging the chr1 centromere
        s = 980_000 + i * 7_000
        rows.append(("p%d" % i, "PERI", "chr1", s, s + 400, "+"))
    rows.append(("cenA", "centromere", "chr1", 990_000, 1_010_000, "0"))
    rows.append(("cenB", "centromere", "chr2", 995_000, 1_005_000, "0"))
    df = pd.DataFrame(
        rows, columns=["gene_id", "gene_family", "chromosome", "start", "end", "strand"]
    )
    df["kind"] = np.where(df["gene_family"] == "centromere", "centromere", "gene")
    return df, genome


@pytest.fixture
def synthetic():
    genome = {"chr1": 1_000_000, "chr2": 1_000_000}
    rows = []
    # family TEL: tight tandem array right at the start of chr1 (telomere-proximal)
    for i in range(6):
        s = 1000 + i * 2000
        rows.append(("tel%d" % i, "TEL", "chr1", s, s + 500, "+"))
    # family MID: spread around the middle of chr2 (interior)
    for i in range(6):
        s = 400_000 + i * 20_000
        rows.append(("mid%d" % i, "MID", "chr2", s, s + 500, "-"))
    df = pd.DataFrame(
        rows, columns=["gene_id", "gene_family", "chromosome", "start", "end", "strand"]
    )
    df["kind"] = "gene"
    return df, genome


def test_family_summary(synthetic):
    df, genome = synthetic
    res = compute_analytics(df, genome, n_permutations=200, seed=1)
    fs = res.family_summary.set_index("gene_family")
    assert fs.loc["TEL", "n_genes"] == 6
    assert fs.loc["MID", "n_genes"] == 6
    assert fs.loc["TEL", "frac_subtelomeric"] == 1.0
    assert fs.loc["MID", "frac_subtelomeric"] == 0.0


def test_tandem_array_detection(synthetic):
    df, genome = synthetic
    res = compute_analytics(df, genome, cluster_gap=5000, n_permutations=50)
    arr = res.tandem_arrays
    tel = arr[arr["gene_family"] == "TEL"]
    assert len(tel) == 1
    assert tel.iloc[0]["n_genes"] == 6
    summ = res.array_summary.set_index("gene_family")
    assert summ.loc["TEL", "frac_clustered"] == 1.0


def test_telomere_bias_direction(synthetic):
    df, genome = synthetic
    res = compute_analytics(df, genome, n_permutations=500, seed=0)
    tb = res.telomere_bias.set_index("gene_family")
    assert tb.loc["TEL", "observed_mean_norm_dist"] < tb.loc["TEL", "null_mean_norm_dist"]
    assert tb.loc["TEL", "direction"] == "telomere-proximal"
    assert tb.loc["TEL", "p_toward_telomere"] < 0.05


def test_distance_metrics_columns(synthetic):
    df, genome = synthetic
    res = compute_analytics(df, genome, n_permutations=50)
    gm = res.gene_metrics
    assert {"dist_to_telomere", "norm_dist_to_telomere", "subtelomeric"} <= set(gm.columns)
    assert (gm["dist_to_telomere"] >= 0).all()


def test_colocalization_runs(synthetic):
    df, genome = synthetic
    res = compute_analytics(
        df, genome, n_permutations=100, colocalization=True, coloc_window=10_000
    )
    assert not res.colocalization.empty
    assert {"family_a", "family_b", "p_value"} <= set(res.colocalization.columns)


def test_summary_is_json_safe(synthetic, tmp_path):
    import json

    df, genome = synthetic
    res = compute_analytics(df, genome, n_permutations=50)
    files = res.write(str(tmp_path))
    assert any(f.endswith("analytics_summary.json") for f in files)
    summary = json.loads((tmp_path / "analytics_summary.json").read_text())
    assert summary["n_families"] == 2


def test_most_clustered_family_ignores_tiny_families():
    genome = {"chr1": 1_000_000}
    rows = [("a1", "A", "chr1", 1000, 1200, "+"), ("a2", "A", "chr1", 1300, 1500, "+")]
    for i in range(5):
        s = 5000 + i * 40000
        rows.append(("b%d" % i, "B", "chr1", s, s + 200, "+"))
    df = pd.DataFrame(
        rows, columns=["gene_id", "gene_family", "chromosome", "start", "end", "strand"]
    )
    df["kind"] = "gene"
    res = compute_analytics(df, genome, n_permutations=50, cluster_gap=1000)
    # family A is a perfect 2-gene "cluster" but must be excluded (< 3 genes)
    assert res.summary["most_clustered_family"]["gene_family"] == "B"


def test_genes_per_family_placement_split():
    genome = {"chr1": 1_000_000, "chr2": 1_000_000}
    rows = []
    for i in range(4):
        rows.append(("c%d" % i, "F1", "chr1", 1000 + i * 5000, 1200 + i * 5000, "+"))
    rows.append(("u1", "F1", "scaffold_9", 10, 200, "+"))
    rows.append(("u2", "F2", "scaffold_9", 300, 400, "+"))
    df = pd.DataFrame(
        rows, columns=["gene_id", "gene_family", "chromosome", "start", "end", "strand"]
    )
    df["kind"] = "gene"
    df["placed"] = [True] * 4 + [False, False]
    res = compute_analytics(df, genome, n_permutations=50)
    gpf = {r["gene_family"]: r for r in res.summary["genes_per_family"]}
    assert gpf["F1"]["on_chromosomes"] == 4 and gpf["F1"]["on_unplaced"] == 1
    assert gpf["F2"]["on_chromosomes"] == 0 and gpf["F2"]["on_unplaced"] == 1
    assert res.summary["n_genes_on_unplaced"] == 2
    assert res.summary["n_unplaced_contigs"] == 1
    fs = res.family_summary.set_index("gene_family")
    assert fs.loc["F2", "n_on_unplaced"] == 1


def test_most_clustered_family_none_when_all_families_tiny():
    genome = {"chr1": 1_000_000}
    rows = [
        ("a1", "A", "chr1", 1000, 1200, "+"), ("a2", "A", "chr1", 1300, 1500, "+"),
        ("b1", "B", "chr1", 9000, 9200, "+"), ("b2", "B", "chr1", 9300, 9500, "+"),
    ]
    df = pd.DataFrame(
        rows, columns=["gene_id", "gene_family", "chromosome", "start", "end", "strand"]
    )
    df["kind"] = "gene"
    res = compute_analytics(df, genome, n_permutations=50)
    assert res.summary["most_clustered_family"] is None


def test_write_emits_genes_per_family_csv_and_chart(synthetic, tmp_path):
    df, genome = synthetic
    res = compute_analytics(df, genome, n_permutations=30)
    files = res.write(
        str(tmp_path), chart_formats=["pdf", "svg", "png"],
        color_map={"TEL": (0.9, 0.1, 0.1), "MID": (0.1, 0.1, 0.9)},
    )
    assert (tmp_path / "analytics_genes_per_family.csv").exists()
    for ext in ("pdf", "svg", "png"):
        p = tmp_path / ("analytics_genes_per_family." + ext)
        assert p.exists() and p.stat().st_size > 0
    assert any(f.endswith("analytics_genes_per_family.pdf") for f in files)
    rows = (tmp_path / "analytics_genes_per_family.csv").read_text().splitlines()
    assert rows[0] == "gene_family,on_chromosomes,on_unplaced,total"


def test_chart_figure_directly(tmp_path):
    from gfviewer.charts import genes_per_family_figure, save_figure

    rows = [
        {"gene_family": "A", "on_chromosomes": 5, "on_unplaced": 2, "total": 7},
        {"gene_family": "B", "on_chromosomes": 3, "on_unplaced": 0, "total": 3},
    ]
    fig = genes_per_family_figure(rows, color_map={"A": (0.2, 0.6, 0.9)})
    out = save_figure(fig, str(tmp_path / "chart.svg"))
    assert out.endswith(".svg") and (tmp_path / "chart.svg").stat().st_size > 0


def test_summary_has_telomere_bias_table(synthetic):
    df, genome = synthetic
    res = compute_analytics(df, genome, n_permutations=200, seed=0)
    table = res.summary["telomere_bias_table"]
    assert {r["gene_family"] for r in table} == {"TEL", "MID"}
    tel = next(r for r in table if r["gene_family"] == "TEL")
    assert set(tel) >= {
        "n_genes", "observed_mean_norm_dist", "null_mean_norm_dist",
        "direction", "p_toward_telomere", "p_toward_interior", "significant",
        "q_value", "significant_fdr",
    }
    assert tel["direction"] == "telomere-proximal"


# --------------------------------------------------------------------------- #
# new analytics (round 8)
# --------------------------------------------------------------------------- #
def test_bh_fdr_monotone_and_bounded():
    q = _bh_fdr([0.001, 0.5, 0.01, 0.9, 0.02])
    assert (q >= 0).all() and (q <= 1).all()
    # sorting by p and by q must agree (BH preserves the order)
    p = np.array([0.001, 0.5, 0.01, 0.9, 0.02])
    assert list(np.argsort(p)) == list(np.argsort(q))
    assert _bh_fdr([]).size == 0


def test_binom_and_poisson_helpers():
    assert _binom_test(5, 10, 0.5) == pytest.approx(1.0, abs=1e-9)
    assert _binom_test(10, 10, 0.5) == pytest.approx(2 * 0.5 ** 10, rel=1e-6)
    assert 0.0 < _poisson_sf(8, 2.0) < 0.02
    assert _poisson_sf(0, 5.0) == 1.0


def test_new_tables_present_and_shaped(rich):
    df, genome = rich
    res = compute_analytics(df, genome, n_permutations=150, seed=0)

    for name in ("ripley", "duplication_modes", "chromosome_enrichment",
                 "strand_bias", "chromosome_richness", "positional_profile",
                 "hotspots", "centromere_bias", "arm_bias"):
        t = getattr(res, name)
        assert t is not None and not t.empty, name
        assert "q_value" in t.columns or "q_clustered" in t.columns or name in (
            "duplication_modes", "chromosome_richness", "positional_profile", "hotspots"
        )

    # FDR columns on the per-family / per-pair tests
    assert "q_value" in res.telomere_bias.columns
    assert "q_value" in res.centromere_bias.columns
    assert {"significant_fdr"} <= set(res.strand_bias.columns)

    # gene-level enrichments
    gm = res.gene_metrics
    assert "dup_mode" in gm.columns
    assert {"arm", "dist_to_centromere"} <= set(gm.columns)

    # chromosome_richness: one row per chromosome
    assert set(res.chromosome_richness["chromosome"]) == set(genome)

    # positional_profile: each family's fractions sum to ~1
    prof = res.positional_profile
    for fam, g in prof.groupby("gene_family"):
        assert g["frac_of_family"].sum() == pytest.approx(1.0, abs=1e-6)

    # family_proximity: square, symmetric, ordered
    mat = res.family_proximity
    assert list(mat.index) == list(mat.columns)
    m = mat.to_numpy(float)
    assert np.allclose(np.nan_to_num(m), np.nan_to_num(m.T))
    assert sorted(res.summary["family_proximity_order"]) == sorted(mat.index)
    assert res.summary["family_proximity_clusters"]


def test_ripley_flags_the_tandem_family(rich):
    df, genome = rich
    res = compute_analytics(df, genome, n_permutations=200, seed=0)
    rip = res.ripley
    # from a few kb up to ~50 kb the tight array is strongly clustered
    tand = rip[(rip["gene_family"] == "TAND") & rip["scale_bp"].between(5000, 50000)]
    assert (tand["L_minus_t_obs"] > 0).all()
    assert tand["significant_fdr"].any()
    peak = next(r for r in res.summary["ripley_clustering"] if r["gene_family"] == "TAND")
    assert peak["significant_clustered_scales_bp"]


def test_duplication_modes_and_hotspots(rich):
    df, genome = rich
    res = compute_analytics(df, genome, n_permutations=100, seed=0)
    dm = res.duplication_modes.set_index("gene_family")
    assert dm.loc["TAND", "predominant_mode"] == "tandem"
    assert dm.loc["SPREAD", "n_tandem"] == 0

    hs = res.hotspots
    assert (hs["chromosome"] == "chr1").any()          # the TAND array is a hotspot
    assert (hs["min_q_value"] <= 0.05).all()


def test_centromere_and_arm_bias(rich):
    df, genome = rich
    res = compute_analytics(df, genome, n_permutations=300, seed=0)
    cb = res.centromere_bias.set_index("gene_family")
    assert cb.loc["PERI", "direction"] == "pericentromeric"
    assert cb.loc["PERI", "p_pericentromeric"] < 0.05
    assert "arm" in res.gene_metrics.columns
    assert set(res.arm_bias["gene_family"]) <= {"TAND", "SPREAD", "PERI"}


def test_write_emits_new_analytics_and_charts(rich, tmp_path):
    df, genome = rich
    res = compute_analytics(df, genome, n_permutations=60, seed=0)
    written = res.write(str(tmp_path), chart_formats=["png"])
    names = {__import__("os").path.basename(p) for p in written}
    for f in ("analytics_ripley.csv", "analytics_duplication_modes.csv",
              "analytics_chromosome_enrichment.csv", "analytics_strand_bias.csv",
              "analytics_chromosome_richness.csv", "analytics_positional_profile.csv",
              "analytics_family_proximity.csv", "analytics_hotspots.csv",
              "analytics_hotspots.bed", "analytics_centromere_bias.csv",
              "analytics_arm_bias.csv"):
        assert f in names, f
    for f in ("analytics_positional_profile.png", "analytics_ripley.png",
              "analytics_family_proximity.png"):
        assert (tmp_path / f).exists() and (tmp_path / f).stat().st_size > 0
    bed = (tmp_path / "analytics_hotspots.bed").read_text().splitlines()
    assert bed[0].startswith("#chrom\t")
    assert len(bed[1].split("\t")) == 6
