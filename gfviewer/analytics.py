"""Quantitative analytics on multigene-family localization.

Given the normalised feature table (from :func:`gfviewer.io.load_features`) and
the genome lengths, :func:`compute_analytics` returns a bundle of tidy
``pandas`` tables plus a JSON-friendly ``summary`` dict.

Tables produced (see the *Help* page for column-by-column documentation):

``genes_per_family``      genes on chromosomes vs. unplaced, per family
``family_summary``        counts, density, gene length, strand fraction
``family_by_chromosome``  one row per (family, chromosome)
``gene_metrics``          one row per gene: distances, arm, duplication mode
``telomere_bias``         permutation test of telomere proximity per family
``centromere_bias``       permutation test of centromere proximity (needs a
                          centromere track)
``arm_bias``              p-arm vs. q-arm occupancy test (needs centromeres)
``tandem_arrays``         detected tandem arrays (>= 2 neighbouring paralogues)
``array_summary``         per-family clustering summary
``duplication_modes``     tandem / proximal / dispersed split per family
``ripley``                1-D Ripley's K/L multi-scale clustering per family
``chromosome_enrichment`` per (family, chromosome) over/under-representation, FDR
``strand_bias``           per-family +/- strand binomial test, FDR
``chromosome_richness``   family diversity of each chromosome
``positional_profile``    binned gene density along the normalised chromosome
``family_proximity``      family x family mean cross-nearest-neighbour distance
``hotspots``              merged windows unusually dense in multigene-family genes
``colocalization``        optional pairwise family co-localization test

Every per-family / per-pair / per-window test carries a Benjamini-Hochberg
``q_value`` (and a ``significant_fdr`` flag) alongside the raw p-value.

Only :mod:`numpy` / :mod:`pandas` are used -- no SciPy dependency.
"""

import json
import math
import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class AnalyticsResult:
    family_summary: pd.DataFrame
    family_by_chromosome: pd.DataFrame
    gene_metrics: pd.DataFrame
    telomere_bias: pd.DataFrame
    tandem_arrays: pd.DataFrame
    array_summary: pd.DataFrame
    genes_per_family: pd.DataFrame = field(default_factory=pd.DataFrame)
    colocalization: pd.DataFrame = field(default_factory=pd.DataFrame)
    ripley: pd.DataFrame = field(default_factory=pd.DataFrame)
    duplication_modes: pd.DataFrame = field(default_factory=pd.DataFrame)
    chromosome_enrichment: pd.DataFrame = field(default_factory=pd.DataFrame)
    strand_bias: pd.DataFrame = field(default_factory=pd.DataFrame)
    chromosome_richness: pd.DataFrame = field(default_factory=pd.DataFrame)
    centromere_bias: pd.DataFrame = field(default_factory=pd.DataFrame)
    arm_bias: pd.DataFrame = field(default_factory=pd.DataFrame)
    positional_profile: pd.DataFrame = field(default_factory=pd.DataFrame)
    family_proximity: pd.DataFrame = field(default_factory=pd.DataFrame)
    hotspots: pd.DataFrame = field(default_factory=pd.DataFrame)
    summary: dict = field(default_factory=dict)

    # tables written as flat CSV (family_proximity is a matrix -> handled apart)
    def tables(self):
        return {
            "genes_per_family": self.genes_per_family,
            "family_summary": self.family_summary,
            "family_by_chromosome": self.family_by_chromosome,
            "gene_metrics": self.gene_metrics,
            "telomere_bias": self.telomere_bias,
            "centromere_bias": self.centromere_bias,
            "arm_bias": self.arm_bias,
            "tandem_arrays": self.tandem_arrays,
            "array_summary": self.array_summary,
            "duplication_modes": self.duplication_modes,
            "ripley": self.ripley,
            "chromosome_enrichment": self.chromosome_enrichment,
            "strand_bias": self.strand_bias,
            "chromosome_richness": self.chromosome_richness,
            "positional_profile": self.positional_profile,
            "hotspots": self.hotspots,
            "colocalization": self.colocalization,
        }

    def write(self, outdir, chart_formats=None, color_map=None, dpi=200,
              chart_titles=True):
        """Write every non-empty table as CSV plus ``analytics_summary.json``.

        Also writes ``analytics_family_proximity.csv`` (a matrix, with its family
        index kept) and ``analytics_hotspots.bed``.  If *chart_formats* is given,
        renders the analytics figures (genes-per-family bar chart, positional
        density profile, Ripley's L curve, family-proximity heat map) as
        ``analytics_<name>.<fmt>`` for each format.
        """
        os.makedirs(outdir, exist_ok=True)
        written = []
        for name, df in self.tables().items():
            if df is None or getattr(df, "empty", True):
                continue
            path = os.path.join(outdir, "analytics_{}.csv".format(name))
            df.to_csv(path, index=False)
            written.append(path)

        if self.family_proximity is not None and not self.family_proximity.empty:
            path = os.path.join(outdir, "analytics_family_proximity.csv")
            self.family_proximity.to_csv(path)          # keep the family index
            written.append(path)

        if self.hotspots is not None and not self.hotspots.empty:
            bed = os.path.join(outdir, "analytics_hotspots.bed")
            _write_hotspot_bed(self.hotspots, bed)
            written.append(bed)

        spath = os.path.join(outdir, "analytics_summary.json")
        with open(spath, "w") as fh:
            json.dump(self.summary, fh, indent=2, default=_json_default)
        written.append(spath)

        if chart_formats:
            from gfviewer import charts as _charts

            if self.summary.get("genes_per_family"):
                written += _charts.write_genes_per_family(
                    self.summary["genes_per_family"], outdir,
                    formats=chart_formats, color_map=color_map, dpi=dpi,
                    titles=chart_titles,
                )
            if not self.positional_profile.empty:
                written += _charts.write_positional_profile(
                    self.positional_profile, outdir,
                    formats=chart_formats, color_map=color_map, dpi=dpi,
                    titles=chart_titles,
                )
            if not self.ripley.empty:
                written += _charts.write_ripley(
                    self.ripley, outdir,
                    formats=chart_formats, color_map=color_map, dpi=dpi,
                    titles=chart_titles,
                )
            if not self.family_proximity.empty:
                written += _charts.write_family_proximity(
                    self.family_proximity, outdir,
                    order=self.summary.get("family_proximity_order"),
                    clusters=self.summary.get("family_proximity_clusters"),
                    formats=chart_formats, dpi=dpi, titles=chart_titles,
                )
        return written


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)


# --------------------------------------------------------------------------- #
# small statistics helpers (no SciPy)
# --------------------------------------------------------------------------- #
def _bh_fdr(pvals):
    """Benjamini-Hochberg q-values for a 1-D array of p-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1.0)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]     # enforce monotonicity
    q = np.empty(n)
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


_LGAMMA = np.frompyfunc(math.lgamma, 1, 1)


def _binom_test(k, n, p):
    """Two-sided binomial test p-value for *k* successes in *n* trials
    (exact for n <= 2000, normal approximation beyond).  No SciPy."""
    if n <= 0:
        return 1.0
    k = int(k)
    p = min(max(float(p), 0.0), 1.0)
    if p <= 0.0 or p >= 1.0:
        return 1.0 if k == int(round(n * p)) else 0.0
    if n <= 2000:
        i = np.arange(n + 1)
        log_choose = (
            math.lgamma(n + 1)
            - _LGAMMA(i + 1).astype(float)
            - _LGAMMA(n - i + 1).astype(float)
        )
        logpmf = log_choose + i * math.log(p) + (n - i) * math.log1p(-p)
        pmf = np.exp(logpmf - logpmf.max())
        pmf /= pmf.sum()
        obs = pmf[k]
        return float(min(1.0, pmf[pmf <= obs * (1 + 1e-7)].sum()))
    mu = n * p
    sd = math.sqrt(n * p * (1 - p))
    if sd == 0:
        return 1.0
    z = abs(k - mu) / sd
    return float(min(1.0, math.erfc(z / math.sqrt(2.0))))


def _poisson_sf(k, mu):
    """Upper-tail Poisson probability  P(X >= k)  for mean *mu*."""
    if mu <= 0:
        return 1.0 if k <= 0 else 0.0
    if k <= 0:
        return 1.0
    if mu < 400 and k < 1500:
        term = math.exp(-mu)
        cdf = term
        for i in range(1, int(k)):
            term *= mu / i
            cdf += term
        return float(min(1.0, max(0.0, 1.0 - cdf)))
    z = (k - 0.5 - mu) / math.sqrt(mu)
    return float(min(1.0, 0.5 * math.erfc(z / math.sqrt(2.0))))


# --------------------------------------------------------------------------- #
def compute_analytics(
    features,
    genome,
    subtelomere_fraction=0.10,
    subtelomere_bp=None,
    cluster_gap=50000,
    proximal_window=None,
    n_permutations=1000,
    seed=0,
    colocalization=False,
    coloc_window=50000,
    ripley_scales=None,
    hotspot_window=100000,
    hotspot_step=None,
    n_proximity_clusters=None,
    alpha=0.05,
):
    """Compute the full analytics bundle.  See the module docstring for outputs."""
    all_genes = features[features["kind"] == "gene"].copy()
    if all_genes.empty:
        raise ValueError("No gene-family features to analyse.")
    if "placed" not in all_genes.columns:
        all_genes["placed"] = True

    # chromosome-level metrics are only meaningful for placed genes
    genes = all_genes[all_genes["placed"]].copy()
    if genes.empty:
        raise ValueError(
            "Every gene-family feature is on an unplaced contig; no "
            "chromosome-level analysis is possible."
        )

    proximal_window = int(proximal_window or 5 * cluster_gap)
    hotspot_step = int(hotspot_step or max(1, hotspot_window // 2))

    genes["mid"] = (genes["start"] + genes["end"]) // 2
    genes["length"] = genes["end"] - genes["start"] + 1
    genes["chrom_length"] = genes["chromosome"].map(genome).astype("int64")

    # distance to the nearer chromosome end (proxy for telomere proximity)
    left = genes["start"] - 1
    right = genes["chrom_length"] - genes["end"]
    genes["dist_to_telomere"] = np.minimum(left, right).clip(lower=0)
    genes["norm_dist_to_telomere"] = genes["dist_to_telomere"] / (genes["chrom_length"] / 2.0)

    if subtelomere_bp is not None:
        thr = pd.Series(float(subtelomere_bp), index=genes.index)
    else:
        thr = genes["chrom_length"] * float(subtelomere_fraction)
    genes["subtelomeric"] = genes["dist_to_telomere"] <= thr

    genome_mb = sum(genome.values()) / 1e6

    # ---- centromere geometry (optional) ------------------------------------
    cen_df = _centromere_frame(features, genome)
    if cen_df is not None:
        genes = _centromere_metrics(genes, cen_df)

    # ---- per-gene duplication mode --------------------------------------- #
    dup_modes, dup_map = _duplication_modes(genes, cluster_gap, proximal_window)
    genes["dup_mode"] = genes["gene_id"].map(dup_map).fillna("dispersed")

    # ---- core tables ---------------------------------------------------- #
    family_by_chrom = _family_by_chromosome(genes, genome)
    family_summary = _family_summary(genes, genome_mb)
    family_summary = _merge_placement(family_summary, all_genes)
    telomere_bias = _telomere_bias(genes, genome, n_permutations, seed, alpha)
    arrays, array_summary = _tandem_arrays(genes, cluster_gap)
    strand_bias = _strand_bias(genes, alpha)
    chrom_enrich = _chromosome_enrichment(genes, genome, alpha)
    chrom_rich = _chromosome_richness(genes, genome)
    profile = _positional_profile(genes)
    ripley = _ripley(
        genes, genome, ripley_scales or _DEFAULT_RIPLEY_SCALES,
        min(int(n_permutations), _RIPLEY_MAX_PERM), seed, alpha,
    )
    prox_mat, prox_order, prox_clusters = _family_proximity(
        genes, n_proximity_clusters
    )
    hotspots = _hotspots(genes, genome, hotspot_window, hotspot_step, alpha)

    centromere_bias = pd.DataFrame()
    arm_bias = pd.DataFrame()
    if cen_df is not None:
        centromere_bias = _centromere_bias(genes, cen_df, n_permutations, seed, alpha)
        arm_bias = _arm_bias(genes, cen_df, genome, alpha)

    coloc = pd.DataFrame()
    if colocalization:
        coloc = _colocalization(genes, genome, coloc_window, n_permutations, seed, alpha)

    summary = _build_summary(
        all_genes, genes, genome, family_summary, telomere_bias, array_summary,
        coloc, dup_modes, ripley, chrom_enrich, strand_bias, centromere_bias,
        arm_bias, hotspots, prox_order, prox_clusters, alpha,
    )
    gpf_df = pd.DataFrame(
        summary["genes_per_family"],
        columns=["gene_family", "on_chromosomes", "on_unplaced", "total"],
    )

    gene_cols = [
        "gene_id", "gene_family", "chromosome", "start", "end", "strand",
        "mid", "length", "dist_to_telomere", "norm_dist_to_telomere",
        "subtelomeric", "dup_mode",
    ]
    if cen_df is not None:
        gene_cols += ["dist_to_centromere", "norm_dist_to_centromere", "arm"]

    return AnalyticsResult(
        family_summary=family_summary,
        family_by_chromosome=family_by_chrom,
        gene_metrics=genes[gene_cols].reset_index(drop=True),
        telomere_bias=telomere_bias,
        tandem_arrays=arrays,
        array_summary=array_summary,
        genes_per_family=gpf_df,
        colocalization=coloc,
        ripley=ripley,
        duplication_modes=dup_modes,
        chromosome_enrichment=chrom_enrich,
        strand_bias=strand_bias,
        chromosome_richness=chrom_rich,
        centromere_bias=centromere_bias,
        arm_bias=arm_bias,
        positional_profile=profile,
        family_proximity=prox_mat,
        hotspots=hotspots,
        summary=summary,
    )


# --------------------------------------------------------------------------- #
def _family_by_chromosome(genes, genome):
    grp = genes.groupby(["gene_family", "chromosome"], sort=True)
    out = grp.agg(
        n_genes=("gene_id", "size"),
        first_start=("start", "min"),
        last_end=("end", "max"),
        n_subtelomeric=("subtelomeric", "sum"),
    ).reset_index()
    out["chrom_length"] = out["chromosome"].map(genome).astype("int64")
    out["density_per_mb"] = out["n_genes"] / (out["chrom_length"] / 1e6)
    out["span_bp"] = out["last_end"] - out["first_start"] + 1
    return out


def _family_summary(genes, genome_mb):
    grp = genes.groupby("gene_family", sort=True)
    out = grp.agg(
        n_genes=("gene_id", "size"),
        n_chromosomes=("chromosome", "nunique"),
        mean_gene_length=("length", "mean"),
        median_gene_length=("length", "median"),
        n_subtelomeric=("subtelomeric", "sum"),
        mean_norm_dist_to_telomere=("norm_dist_to_telomere", "mean"),
    ).reset_index()
    plus = grp["strand"].apply(lambda s: float((s == "+").mean())).reset_index(name="frac_plus_strand")
    out = out.merge(plus, on="gene_family")
    out["genes_per_mb_genome"] = out["n_genes"] / genome_mb
    out["frac_subtelomeric"] = out["n_subtelomeric"] / out["n_genes"]
    return out.sort_values("n_genes", ascending=False).reset_index(drop=True)


def _merge_placement(family_summary, all_genes):
    """Add ``n_on_chromosomes`` / ``n_on_unplaced`` (from the full gene set) to
    the per-family summary, including families seen only on unplaced contigs."""
    pl = all_genes.groupby("gene_family", sort=True)["placed"].agg(
        n_on_chromosomes="sum", n_total="count"
    ).reset_index()
    pl["n_on_chromosomes"] = pl["n_on_chromosomes"].astype(int)
    pl["n_on_unplaced"] = (pl["n_total"] - pl["n_on_chromosomes"]).astype(int)
    merged = family_summary.merge(
        pl[["gene_family", "n_on_chromosomes", "n_on_unplaced"]],
        on="gene_family", how="outer",
    )
    for col in ("n_on_chromosomes", "n_on_unplaced"):
        merged[col] = merged[col].fillna(0).astype(int)
    if "n_genes" in merged:
        merged["n_genes"] = merged["n_genes"].fillna(
            merged["n_on_chromosomes"] + merged["n_on_unplaced"]
        ).astype(int)
    return merged.sort_values("n_genes", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
def _telomere_bias(genes, genome, n_permutations, seed, alpha):
    """Permutation test: is a family closer to / further from telomeres than
    expected if its genes were placed uniformly at random on the same
    chromosomes (chromosome assignment and gene lengths preserved)?"""
    rng = np.random.default_rng(seed)
    chrom_len = {c: int(genome[c]) for c in genome}
    rows = []
    for fam, sub in genes.groupby("gene_family", sort=True):
        lengths = sub["length"].to_numpy()
        chroms = sub["chromosome"].to_numpy()
        clen = np.array([chrom_len[c] for c in chroms], dtype=float)
        obs = float(sub["norm_dist_to_telomere"].mean())

        n = len(sub)
        span = np.maximum(clen - lengths, 1.0)
        starts = rng.random((n_permutations, n)) * span  # 0-based
        ends = starts + lengths
        dist = np.minimum(starts, clen - ends)
        dist = np.clip(dist, 0, None)
        null_means = (dist / (clen / 2.0)).mean(axis=1)

        p_tel = float((null_means <= obs).mean())
        p_int = float((null_means >= obs).mean())
        direction = "telomere-proximal" if obs < null_means.mean() else "interior"
        rows.append(
            {
                "gene_family": fam,
                "n_genes": n,
                "observed_mean_norm_dist": obs,
                "null_mean_norm_dist": float(null_means.mean()),
                "null_sd": float(null_means.std(ddof=1)) if n_permutations > 1 else 0.0,
                "p_toward_telomere": p_tel,
                "p_toward_interior": p_int,
                "direction": direction,
                "significant": bool(min(p_tel, p_int) <= alpha),
                "n_permutations": n_permutations,
            }
        )
    df = pd.DataFrame(rows)
    if not df.empty:
        two = np.clip(2.0 * np.minimum(df["p_toward_telomere"], df["p_toward_interior"]), 0, 1)
        df["q_value"] = _bh_fdr(two)
        df["significant_fdr"] = df["q_value"] <= alpha
    return df.sort_values("p_toward_telomere").reset_index(drop=True)


# --------------------------------------------------------------------------- #
def _centromere_frame(features, genome):
    cen = features[features["kind"] == "centromere"]
    cen = cen[cen["chromosome"].isin(genome)]
    if cen.empty:
        return None
    g = cen.groupby("chromosome").agg(
        cen_start=("start", "min"), cen_end=("end", "max")
    )
    g["cen_mid"] = (g["cen_start"] + g["cen_end"]) / 2.0
    return g


def _centromere_metrics(genes, cen_df):
    genes = genes.copy()
    cs = genes["chromosome"].map(cen_df["cen_start"])
    ce = genes["chromosome"].map(cen_df["cen_end"])
    have = cs.notna()
    mid = genes["mid"].astype(float)
    d = np.where(mid < cs, cs - mid, np.where(mid > ce, mid - ce, 0.0))
    genes["dist_to_centromere"] = np.where(have, d, np.nan)
    genes["norm_dist_to_centromere"] = genes["dist_to_centromere"] / (genes["chrom_length"] / 2.0)
    genes["arm"] = np.where(
        ~have, "na",
        np.where(mid < cs, "p", np.where(mid > ce, "q", "centromeric")),
    )
    return genes


def _centromere_bias(genes, cen_df, n_permutations, seed, alpha):
    sub_all = genes[genes["chromosome"].isin(cen_df.index)].copy()
    sub_all = sub_all[sub_all["norm_dist_to_centromere"].notna()]
    if sub_all.empty:
        return pd.DataFrame()
    rng = np.random.default_rng(seed + 11)
    rows = []
    for fam, sub in sub_all.groupby("gene_family", sort=True):
        n = len(sub)
        if n < 2:
            continue
        L = sub["chrom_length"].to_numpy().astype(float)
        cs = sub["chromosome"].map(cen_df["cen_start"]).to_numpy().astype(float)
        ce = sub["chromosome"].map(cen_df["cen_end"]).to_numpy().astype(float)
        lengths = sub["length"].to_numpy().astype(float)
        obs = float(sub["norm_dist_to_centromere"].mean())
        span = np.maximum(L - lengths, 1.0)
        starts = rng.random((n_permutations, n)) * span
        mids = starts + lengths / 2.0
        dd = np.where(mids < cs, cs - mids, np.where(mids > ce, mids - ce, 0.0))
        null = (dd / (L / 2.0)).mean(axis=1)
        p_peri = float((null <= obs).mean())
        p_dist = float((null >= obs).mean())
        rows.append(
            {
                "gene_family": fam, "n_genes": n,
                "observed_mean_norm_dist_cen": obs,
                "null_mean_norm_dist_cen": float(null.mean()),
                "direction": "pericentromeric" if obs < null.mean() else "centromere-distal",
                "p_pericentromeric": p_peri, "p_distal": p_dist,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    two = np.clip(2.0 * np.minimum(df["p_pericentromeric"], df["p_distal"]), 0, 1)
    df["q_value"] = _bh_fdr(two)
    df["significant_fdr"] = df["q_value"] <= alpha
    return df.sort_values("p_pericentromeric").reset_index(drop=True)


def _arm_bias(genes, cen_df, genome, alpha):
    sub_all = genes[genes["arm"].isin(["p", "q"])]
    if sub_all.empty:
        return pd.DataFrame()
    rows = []
    for fam, sub in sub_all.groupby("gene_family", sort=True):
        p_len = q_len = 0.0
        for c in sub["chromosome"].unique():
            if c not in cen_df.index:
                continue
            p_len += float(cen_df.loc[c, "cen_start"])
            q_len += float(genome[c]) - float(cen_df.loc[c, "cen_end"])
        tot = p_len + q_len
        if tot <= 0:
            continue
        exp_p = p_len / tot
        n_p = int((sub["arm"] == "p").sum())
        n_q = int((sub["arm"] == "q").sum())
        n = n_p + n_q
        if n == 0:
            continue
        pval = _binom_test(n_p, n, exp_p)
        rows.append(
            {
                "gene_family": fam, "n_p_arm": n_p, "n_q_arm": n_q,
                "expected_frac_p": round(exp_p, 4),
                "observed_frac_p": round(n_p / n, 4),
                "p_value": pval,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["q_value"] = _bh_fdr(df["p_value"].to_numpy())
    df["significant_fdr"] = df["q_value"] <= alpha
    df["biased_toward_arm"] = np.where(
        df["observed_frac_p"] >= df["expected_frac_p"], "p", "q"
    )
    return df.sort_values("q_value").reset_index(drop=True)


# --------------------------------------------------------------------------- #
def _tandem_arrays(genes, cluster_gap):
    arrays = []
    per_family = {}
    for (fam, chrom), sub in genes.sort_values(["gene_family", "chromosome", "start"]).groupby(
        ["gene_family", "chromosome"], sort=True
    ):
        starts = sub["start"].to_numpy()
        ends = sub["end"].to_numpy()
        ids = sub["gene_id"].to_numpy()
        cluster_start = starts[0]
        cluster_end = ends[0]
        members = [ids[0]]
        for i in range(1, len(sub)):
            if starts[i] - cluster_end <= cluster_gap:
                cluster_end = max(cluster_end, ends[i])
                members.append(ids[i])
            else:
                _flush_array(arrays, per_family, fam, chrom, cluster_start, cluster_end, members)
                cluster_start, cluster_end, members = starts[i], ends[i], [ids[i]]
        _flush_array(arrays, per_family, fam, chrom, cluster_start, cluster_end, members)

    arrays_df = pd.DataFrame(
        arrays,
        columns=["gene_family", "chromosome", "start", "end", "n_genes", "span_bp", "gene_ids"],
    )

    rows = []
    for fam, sub in genes.groupby("gene_family", sort=True):
        info = per_family.get(fam, {"n_arrays": 0, "genes_in_arrays": 0, "largest": 0})
        total = len(sub)
        rows.append(
            {
                "gene_family": fam,
                "n_genes": total,
                "n_tandem_arrays": info["n_arrays"],
                "genes_in_arrays": info["genes_in_arrays"],
                "frac_clustered": info["genes_in_arrays"] / total if total else 0.0,
                "largest_array": info["largest"],
            }
        )
    return arrays_df, pd.DataFrame(rows).sort_values("frac_clustered", ascending=False).reset_index(drop=True)


def _flush_array(arrays, per_family, fam, chrom, start, end, members):
    info = per_family.setdefault(fam, {"n_arrays": 0, "genes_in_arrays": 0, "largest": 0})
    if len(members) >= 2:
        arrays.append(
            {
                "gene_family": fam,
                "chromosome": chrom,
                "start": int(start),
                "end": int(end),
                "n_genes": len(members),
                "span_bp": int(end - start + 1),
                "gene_ids": ";".join(map(str, members)),
            }
        )
        info["n_arrays"] += 1
        info["genes_in_arrays"] += len(members)
        info["largest"] = max(info["largest"], len(members))


# --------------------------------------------------------------------------- #
def _duplication_modes(genes, tandem_gap, proximal_window):
    """Classify each gene as tandem / proximal / dispersed from the gap to its
    nearest same-family paralogue on the same chromosome."""
    mode = {}
    ordered = genes.sort_values(["gene_family", "chromosome", "start"])
    for (_fam, _chrom), sub in ordered.groupby(["gene_family", "chromosome"], sort=False):
        s = sub["start"].to_numpy()
        e = sub["end"].to_numpy()
        ids = sub["gene_id"].to_numpy()
        n = len(sub)
        for i in range(n):
            gaps = []
            if i > 0:
                gaps.append(s[i] - e[i - 1])
            if i < n - 1:
                gaps.append(s[i + 1] - e[i])
            gap = min(gaps) if gaps else None
            if gap is not None and gap <= tandem_gap:
                mode[ids[i]] = "tandem"
            elif gap is not None and gap <= proximal_window:
                mode[ids[i]] = "proximal"
            else:
                mode[ids[i]] = "dispersed"

    tmp = genes.assign(_mode=genes["gene_id"].map(mode).fillna("dispersed"))
    counts = (
        tmp.groupby("gene_family")["_mode"].value_counts().unstack(fill_value=0)
    )
    for col in ("tandem", "proximal", "dispersed"):
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[["tandem", "proximal", "dispersed"]].reset_index()
    counts["n_genes"] = counts[["tandem", "proximal", "dispersed"]].sum(axis=1)
    for col in ("tandem", "proximal", "dispersed"):
        counts["frac_" + col] = counts[col] / counts["n_genes"]
    counts["predominant_mode"] = counts[["tandem", "proximal", "dispersed"]].idxmax(axis=1)
    counts = counts.rename(
        columns={"tandem": "n_tandem", "proximal": "n_proximal", "dispersed": "n_dispersed"}
    )
    counts = counts.sort_values("n_genes", ascending=False).reset_index(drop=True)
    return counts, mode


# --------------------------------------------------------------------------- #
_DEFAULT_RIPLEY_SCALES = [
    1_000, 2_000, 5_000, 10_000, 20_000, 50_000,
    100_000, 200_000, 500_000, 1_000_000,
]
_RIPLEY_MIN_N = 5
_RIPLEY_MAX_N = 400          # subsample larger families for the K computation
_RIPLEY_MAX_PERM = 499


def _ripley_k_1d(pos, L, scales):
    """Edge-corrected 1-D Ripley K at each *scale* for points *pos* on [0, L]."""
    n = len(pos)
    if n < 2:
        return np.zeros(len(scales))
    d = np.abs(pos[:, None] - pos[None, :])
    np.fill_diagonal(d, np.inf)
    out = np.empty(len(scales))
    for si, t in enumerate(scales):
        lo = np.maximum(pos - t, 0.0)
        hi = np.minimum(pos + t, L)
        w = np.clip((hi - lo) / (2.0 * t), 1e-6, 1.0)     # isotropic edge weight
        contrib = ((d <= t) / w[:, None]).sum()
        out[si] = L * contrib / (n * (n - 1))
    return out


def _ripley(genes, genome, scales, n_perm, seed, alpha):
    rng = np.random.default_rng(seed + 7)
    scales = np.array(sorted(scales), dtype=float)
    rows = []
    for fam, sub in genes.groupby("gene_family", sort=True):
        if len(sub) < _RIPLEY_MIN_N:
            continue
        posmap = {}
        for c, g in sub.groupby("chromosome"):
            pos = np.sort(g["mid"].to_numpy().astype(float))
            if len(pos) > _RIPLEY_MAX_N:
                pos = pos[np.linspace(0, len(pos) - 1, _RIPLEY_MAX_N).astype(int)]
            if len(pos) >= 2:
                posmap[c] = pos
        if not posmap:
            continue
        min_L = min(float(genome[c]) for c in posmap)
        use = scales[scales < 0.45 * min_L]
        if use.size == 0:
            continue

        def pooled_L(pm):
            num = np.zeros(use.size)
            den = 0.0
            for c, pos in pm.items():
                nc = len(pos)
                w = nc * (nc - 1)
                num += w * _ripley_k_1d(pos, float(genome[c]), use)
                den += w
            return (num / den) / 2.0 if den else np.zeros(use.size)

        stat_obs = pooled_L(posmap) - use
        perm = np.empty((n_perm, use.size))
        for p in range(n_perm):
            pm = {c: np.sort(rng.random(len(pos)) * float(genome[c]))
                  for c, pos in posmap.items()}
            perm[p] = pooled_L(pm) - use
        lo = np.percentile(perm, 100 * alpha / 2.0, axis=0)
        hi = np.percentile(perm, 100 * (1 - alpha / 2.0), axis=0)
        p_clus = (1 + (perm >= stat_obs).sum(axis=0)) / (n_perm + 1.0)
        p_disp = (1 + (perm <= stat_obs).sum(axis=0)) / (n_perm + 1.0)
        for k, t in enumerate(use):
            rows.append(
                {
                    "gene_family": fam,
                    "scale_bp": int(t),
                    "L_minus_t_obs": float(stat_obs[k]),
                    "L_minus_t_null_mean": float(perm[:, k].mean()),
                    "env_lo": float(lo[k]),
                    "env_hi": float(hi[k]),
                    "p_clustered": float(p_clus[k]),
                    "p_dispersed": float(p_disp[k]),
                    "pattern": (
                        "clustered" if stat_obs[k] > hi[k]
                        else "dispersed" if stat_obs[k] < lo[k]
                        else "random"
                    ),
                }
            )
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_clustered"] = _bh_fdr(df["p_clustered"].to_numpy())
        df["significant_fdr"] = df["q_clustered"] <= alpha
    return df


# --------------------------------------------------------------------------- #
def _chromosome_enrichment(genes, genome, alpha):
    """Per (family, chromosome): more/fewer genes than expected if the family
    were spread over the genome in proportion to chromosome length."""
    L_total = float(sum(genome.values()))
    counts = genes.groupby(["gene_family", "chromosome"]).size()
    fam_tot = genes.groupby("gene_family").size()
    rows = []
    for fam, n_f in fam_tot.items():
        for c in genome:
            k = int(counts.get((fam, c), 0))
            p = float(genome[c]) / L_total
            exp = float(n_f) * p
            rows.append(
                {
                    "gene_family": fam,
                    "chromosome": c,
                    "observed": k,
                    "expected": round(exp, 3),
                    "log2_obs_over_exp": float(np.log2((k + 0.5) / (exp + 0.5))),
                    "p_value": _binom_test(k, int(n_f), p),
                }
            )
    df = pd.DataFrame(rows)
    df["q_value"] = _bh_fdr(df["p_value"].to_numpy())
    df["enriched"] = (df["q_value"] <= alpha) & (df["observed"] > df["expected"])
    df["depleted"] = (df["q_value"] <= alpha) & (df["observed"] < df["expected"])
    return df.sort_values(["gene_family", "chromosome"]).reset_index(drop=True)


def _strand_bias(genes, alpha):
    rows = []
    for fam, sub in genes.groupby("gene_family", sort=True):
        npl = int((sub["strand"] == "+").sum())
        nmi = int((sub["strand"] == "-").sum())
        n = npl + nmi
        if n == 0:
            continue
        rows.append(
            {
                "gene_family": fam, "n_plus": npl, "n_minus": nmi,
                "frac_plus": round(npl / n, 4),
                "p_value": _binom_test(npl, n, 0.5),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["q_value"] = _bh_fdr(df["p_value"].to_numpy())
    df["significant_fdr"] = df["q_value"] <= alpha
    df["biased_toward"] = np.where(df["frac_plus"] >= 0.5, "+", "-")
    return df.sort_values("q_value").reset_index(drop=True)


def _chromosome_richness(genes, genome):
    rows = []
    for c in genome:
        sub = genes[genes["chromosome"] == c]
        n = len(sub)
        vc = sub["gene_family"].value_counts()
        if n:
            probs = (vc / n).to_numpy()
            H = float(-(probs * np.log(probs)).sum())
            nfam = int(len(vc))
            even = float(H / np.log(nfam)) if nfam > 1 else 0.0
            dom, domf = str(vc.index[0]), float(vc.iloc[0] / n)
        else:
            H, nfam, even, dom, domf = 0.0, 0, 0.0, None, 0.0
        rows.append(
            {
                "chromosome": c,
                "chrom_length_bp": int(genome[c]),
                "n_genes": n,
                "n_families": nfam,
                "shannon_diversity": round(H, 4),
                "evenness": round(even, 4),
                "dominant_family": dom,
                "dominant_family_frac": round(domf, 4),
                "genes_per_mb": round(n / (genome[c] / 1e6), 4),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
_PROFILE_BINS = 20


def _positional_profile(genes, nbins=_PROFILE_BINS):
    u = (genes["mid"] / genes["chrom_length"]).clip(0, 0.999999)
    b = (u * nbins).astype(int)
    tmp = genes.assign(_bin=b.values)
    rows = []

    def emit(label, series):
        counts = series.value_counts().reindex(range(nbins), fill_value=0)
        tot = int(counts.sum())
        for k in range(nbins):
            rows.append(
                {
                    "gene_family": label,
                    "bin": k,
                    "bin_center": round((k + 0.5) / nbins, 4),
                    "n_genes": int(counts[k]),
                    "frac_of_family": float(counts[k] / tot) if tot else 0.0,
                }
            )

    for fam, sub in tmp.groupby("gene_family", sort=True):
        emit(fam, sub["_bin"])
    emit("(all families)", tmp["_bin"])
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def _family_proximity(genes, k):
    fams = sorted(genes["gene_family"].unique())
    F = len(fams)
    if F < 2:
        return pd.DataFrame(), [], []
    by = {
        f: {
            c: np.sort(g["mid"].to_numpy().astype(float))
            for c, g in genes[genes["gene_family"] == f].groupby("chromosome")
        }
        for f in fams
    }

    def cross_nn(a, b):
        ds = []
        for c, apos in by[a].items():
            bpos = by[b].get(c)
            if bpos is None or len(bpos) == 0:
                continue
            idx = np.clip(np.searchsorted(bpos, apos), 0, len(bpos) - 1)
            cand = np.abs(bpos[idx] - apos)
            idxm = np.clip(idx - 1, 0, len(bpos) - 1)
            cand = np.minimum(cand, np.abs(bpos[idxm] - apos))
            ds.append(cand)
        return float(np.concatenate(ds).mean()) if ds else np.nan

    def within_nn(a):
        dd = []
        for _c, pos in by[a].items():
            if len(pos) < 2:
                continue
            g = np.diff(pos)
            nn = np.minimum(np.r_[g, np.inf], np.r_[np.inf, g])
            dd.append(nn[np.isfinite(nn)])
        return float(np.concatenate(dd).mean()) if dd else np.nan

    D = np.full((F, F), np.nan)
    for i in range(F):
        D[i, i] = within_nn(fams[i])
        for j in range(i + 1, F):
            vals = [v for v in (cross_nn(fams[i], fams[j]), cross_nn(fams[j], fams[i]))
                    if not np.isnan(v)]
            D[i, j] = D[j, i] = float(np.mean(vals)) if vals else np.nan
    mat = pd.DataFrame(D, index=fams, columns=fams)

    if k is None:
        k = min(F, max(2, round(F / 4.0)))
    k = int(min(max(k, 2), F))
    order, clusters = _cluster_families(mat, k)
    return mat, order, clusters


def _cluster_families(mat, k):
    """Average-linkage agglomerative clustering of families from a distance
    matrix.  Returns (leaf order, list of clusters) with the tree cut into *k*."""
    fams = list(mat.index)
    F = len(fams)
    D = mat.to_numpy(dtype=float).copy()
    finite = D[np.isfinite(D)]
    big = float(finite.max() * 1.5) if finite.size else 1.0
    D[~np.isfinite(D)] = big
    np.fill_diagonal(D, 0.0)

    members = {i: [i] for i in range(F)}
    dist = {(a, b): D[a, b] for a in range(F) for b in range(a + 1, F)}
    ids = list(range(F))
    nxt = F
    snapshot = [members[c][:] for c in ids]
    while len(ids) > 1:
        (a, b), _h = min(
            ((pair, d) for pair, d in dist.items()
             if pair[0] in ids and pair[1] in ids),
            key=lambda kv: kv[1],
        )
        na, nb = len(members[a]), len(members[b])
        new = nxt
        nxt += 1
        members[new] = members[a] + members[b]
        for c in ids:
            if c in (a, b):
                continue
            da = dist[(min(a, c), max(a, c))]
            db = dist[(min(b, c), max(b, c))]
            dist[(min(new, c), max(new, c))] = (na * da + nb * db) / (na + nb)
        ids = [c for c in ids if c not in (a, b)] + [new]
        if len(ids) == k:
            snapshot = [members[c][:] for c in ids]
    order = [fams[i] for i in members[ids[0]]]
    clusters = [[fams[i] for i in grp] for grp in snapshot]
    return order, clusters


# --------------------------------------------------------------------------- #
def _hotspots(genes, genome, window, step, alpha, min_genes=3):
    """Sliding-window scan for regions unusually dense in multigene-family genes
    (all families pooled), Poisson-tested against a uniform null and merged."""
    win_rows = []
    for c in genome:
        L = int(genome[c])
        sub = genes[genes["chromosome"] == c]
        ntot = len(sub)
        if ntot == 0 or L <= window:
            continue
        mids = sub["mid"].to_numpy()
        fams = sub["gene_family"].to_numpy()
        lam = ntot * window / L
        for ws in range(0, L - window + 1, step):
            we = ws + window
            m = (mids >= ws) & (mids < we)
            kk = int(m.sum())
            if kk < min_genes or kk <= lam:
                continue
            win_rows.append(
                (c, ws + 1, we, kk, sorted(set(fams[m])), _poisson_sf(kk, lam))
            )
    if not win_rows:
        return pd.DataFrame()
    wdf = pd.DataFrame(
        win_rows, columns=["chromosome", "w_start", "w_end", "n_genes", "fam_list", "p_value"]
    )
    wdf["q_value"] = _bh_fdr(wdf["p_value"].to_numpy())
    sig = wdf[wdf["q_value"] <= alpha]
    if sig.empty:
        return pd.DataFrame()

    merged = []
    for c, g in sig.sort_values(["chromosome", "w_start"]).groupby("chromosome"):
        cur = None
        for r in g.itertuples(index=False):
            if cur and r.w_start <= cur["end"]:
                cur["end"] = max(cur["end"], r.w_end)
                cur["min_q"] = min(cur["min_q"], r.q_value)
            else:
                if cur:
                    merged.append(cur)
                cur = {"chromosome": c, "start": r.w_start, "end": r.w_end, "min_q": r.q_value}
        if cur:
            merged.append(cur)

    out = []
    for i, m in enumerate(merged, 1):
        reg = genes[
            (genes["chromosome"] == m["chromosome"])
            & (genes["mid"] >= m["start"] - 1)
            & (genes["mid"] <= m["end"])
        ]
        span = m["end"] - m["start"] + 1
        out.append(
            {
                "hotspot_id": "hotspot_{}".format(i),
                "chromosome": m["chromosome"],
                "start": int(m["start"]),
                "end": int(m["end"]),
                "length_bp": int(span),
                "n_genes": int(len(reg)),
                "n_families": int(reg["gene_family"].nunique()),
                "families": ",".join(sorted(reg["gene_family"].unique())),
                "genes_per_mb": round(len(reg) / (span / 1e6), 2),
                "min_q_value": float(m["min_q"]),
            }
        )
    return pd.DataFrame(out)


def _write_hotspot_bed(hotspots_df, path):
    with open(path, "w") as fh:
        fh.write("#chrom\tstart\tend\tname\tscore\tstrand\n")
        for r in hotspots_df.itertuples(index=False):
            q = max(float(getattr(r, "min_q_value", 1.0) or 1.0), 1e-300)
            score = min(1000, int(round(-10.0 * math.log10(q))))
            name = "{};n={};fam={}".format(
                r.hotspot_id, int(r.n_genes), str(r.families).replace(" ", "")
            )
            fh.write(
                "{}\t{}\t{}\t{}\t{}\t.\n".format(
                    r.chromosome, int(r.start) - 1, int(r.end), name, score
                )
            )


# --------------------------------------------------------------------------- #
def _colocalization(genes, genome, window, n_permutations, seed, alpha):
    rng = np.random.default_rng(seed + 1)
    fams = sorted(genes["gene_family"].unique())
    by_fam_chrom = {
        (f, c): sub[["start", "end"]].to_numpy()
        for (f, c), sub in genes.groupby(["gene_family", "chromosome"])
    }
    chrom_len = {c: int(genome[c]) for c in genome}
    rows = []
    for i in range(len(fams)):
        for j in range(i + 1, len(fams)):
            fa, fb = fams[i], fams[j]
            obs = _count_close(genes, fa, fb, window)
            if obs == 0:
                rows.append(
                    {
                        "family_a": fa, "family_b": fb, "observed_close_pairs": 0,
                        "null_mean": 0.0, "p_value": 1.0, "significant": False,
                    }
                )
                continue
            null = np.empty(n_permutations)
            a_chroms = genes.loc[genes["gene_family"] == fa, "chromosome"].to_numpy()
            a_len = (
                genes.loc[genes["gene_family"] == fa, "end"].to_numpy()
                - genes.loc[genes["gene_family"] == fa, "start"].to_numpy()
                + 1
            )
            b_arrays = {c: by_fam_chrom.get((fb, c), np.empty((0, 2))) for c in set(a_chroms)}
            for p in range(n_permutations):
                cnt = 0
                for c in set(a_chroms):
                    mask = a_chroms == c
                    L = a_len[mask]
                    span = max(chrom_len[c] - L.max(), 1) if len(L) else 1
                    rs = rng.random(mask.sum()) * span
                    re_ = rs + L
                    bs = b_arrays[c]
                    if len(bs) == 0:
                        continue
                    for s, e in zip(rs, re_):
                        near = (bs[:, 0] - window <= e) & (bs[:, 1] + window >= s)
                        cnt += int(near.any())
                null[p] = cnt
            p_val = float((null >= obs).mean())
            rows.append(
                {
                    "family_a": fa, "family_b": fb, "observed_close_pairs": int(obs),
                    "null_mean": float(null.mean()), "p_value": p_val,
                    "significant": bool(p_val <= alpha),
                }
            )
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_value"] = _bh_fdr(df["p_value"].to_numpy())
        df["significant_fdr"] = df["q_value"] <= alpha
    return df.sort_values("p_value").reset_index(drop=True)


def _count_close(genes, fa, fb, window):
    a = genes[genes["gene_family"] == fa]
    b = genes[genes["gene_family"] == fb]
    total = 0
    for c, asub in a.groupby("chromosome"):
        bsub = b[b["chromosome"] == c]
        if bsub.empty:
            continue
        bs = bsub["start"].to_numpy()
        be = bsub["end"].to_numpy()
        for s, e in zip(asub["start"].to_numpy(), asub["end"].to_numpy()):
            near = (bs - window <= e) & (be + window >= s)
            total += int(near.any())
    return total


# --------------------------------------------------------------------------- #
def _build_summary(all_genes, genes, genome, family_summary, telomere_bias,
                   array_summary, coloc, dup_modes, ripley, chrom_enrich,
                   strand_bias, centromere_bias, arm_bias, hotspots,
                   prox_order, prox_clusters, alpha):
    round4 = lambda v: round(float(v), 4) if isinstance(v, (int, float, np.floating)) else v

    tel_table = [
        {
            "gene_family": r["gene_family"],
            "n_genes": int(r["n_genes"]),
            "observed_mean_norm_dist": round4(r["observed_mean_norm_dist"]),
            "null_mean_norm_dist": round4(r["null_mean_norm_dist"]),
            "direction": r["direction"],
            "p_toward_telomere": round4(r["p_toward_telomere"]),
            "p_toward_interior": round4(r["p_toward_interior"]),
            "q_value": round4(r.get("q_value", np.nan)),
            "significant": bool(r["significant"]),
            "significant_fdr": bool(r.get("significant_fdr", False)),
        }
        for _, r in telomere_bias.iterrows()
    ]

    gpf = all_genes.groupby("gene_family", sort=False)["placed"].agg(
        on_chromosomes="sum", total="count"
    ).reset_index()
    gpf["on_unplaced"] = gpf["total"] - gpf["on_chromosomes"]
    gpf = gpf.sort_values("total", ascending=False)
    genes_per_family = [
        {
            "gene_family": r.gene_family,
            "on_chromosomes": int(r.on_chromosomes),
            "on_unplaced": int(r.on_unplaced),
            "total": int(r.total),
        }
        for r in gpf.itertuples(index=False)
    ]
    n_unplaced = int((~all_genes["placed"]).sum())

    dup_table = [
        {
            "gene_family": r["gene_family"],
            "n_tandem": int(r["n_tandem"]),
            "n_proximal": int(r["n_proximal"]),
            "n_dispersed": int(r["n_dispersed"]),
            "predominant_mode": r["predominant_mode"],
        }
        for _, r in dup_modes.iterrows()
    ]

    ripley_summary = []
    if not ripley.empty:
        for fam, g in ripley.groupby("gene_family"):
            gg = g.sort_values("scale_bp")
            peak = gg.loc[gg["L_minus_t_obs"].idxmax()]
            sig_scales = [int(s) for s in gg.loc[gg["significant_fdr"], "scale_bp"]]
            ripley_summary.append(
                {
                    "gene_family": fam,
                    "peak_scale_bp": int(peak["scale_bp"]),
                    "peak_L_minus_t": round4(peak["L_minus_t_obs"]),
                    "significant_clustered_scales_bp": sig_scales,
                }
            )

    enrich_hits = []
    if not chrom_enrich.empty:
        for _, r in chrom_enrich[chrom_enrich["enriched"] | chrom_enrich["depleted"]].iterrows():
            enrich_hits.append(
                {
                    "gene_family": r["gene_family"],
                    "chromosome": r["chromosome"],
                    "observed": int(r["observed"]),
                    "expected": round4(r["expected"]),
                    "log2_obs_over_exp": round4(r["log2_obs_over_exp"]),
                    "q_value": round4(r["q_value"]),
                    "call": "enriched" if r["enriched"] else "depleted",
                }
            )

    strand_hits = []
    if not strand_bias.empty:
        for _, r in strand_bias[strand_bias["significant_fdr"]].iterrows():
            strand_hits.append(
                {
                    "gene_family": r["gene_family"],
                    "frac_plus": round4(r["frac_plus"]),
                    "biased_toward": r["biased_toward"],
                    "q_value": round4(r["q_value"]),
                }
            )

    cen_table = [
        {
            "gene_family": r["gene_family"],
            "n_genes": int(r["n_genes"]),
            "direction": r["direction"],
            "p_pericentromeric": round4(r["p_pericentromeric"]),
            "p_distal": round4(r["p_distal"]),
            "q_value": round4(r["q_value"]),
            "significant_fdr": bool(r["significant_fdr"]),
        }
        for _, r in centromere_bias.iterrows()
    ] if not centromere_bias.empty else []

    arm_table = [
        {
            "gene_family": r["gene_family"],
            "n_p_arm": int(r["n_p_arm"]),
            "n_q_arm": int(r["n_q_arm"]),
            "expected_frac_p": round4(r["expected_frac_p"]),
            "observed_frac_p": round4(r["observed_frac_p"]),
            "q_value": round4(r["q_value"]),
            "biased_toward_arm": r["biased_toward_arm"],
            "significant_fdr": bool(r["significant_fdr"]),
        }
        for _, r in arm_bias.iterrows()
    ] if not arm_bias.empty else []

    hotspot_regions = []
    if not hotspots.empty:
        for _, r in hotspots.head(20).iterrows():
            hotspot_regions.append(
                {
                    "hotspot_id": r["hotspot_id"],
                    "chromosome": r["chromosome"],
                    "start": int(r["start"]),
                    "end": int(r["end"]),
                    "n_genes": int(r["n_genes"]),
                    "n_families": int(r["n_families"]),
                    "families": r["families"],
                    "genes_per_mb": round4(r["genes_per_mb"]),
                    "min_q_value": round4(r["min_q_value"]),
                }
            )

    prox_clusters_json = [
        {"cluster": i + 1, "families": grp} for i, grp in enumerate(prox_clusters or [])
    ]

    summary = {
        "n_genes": int(len(all_genes)),
        "n_genes_on_chromosomes": int(all_genes["placed"].sum()),
        "n_genes_on_unplaced": n_unplaced,
        "n_unplaced_contigs": int(all_genes.loc[~all_genes["placed"], "chromosome"].nunique()),
        "n_families": int(all_genes["gene_family"].nunique()),
        "n_chromosomes": int(len(genome)),
        "genome_size_bp": int(sum(genome.values())),
        "genes_subtelomeric": int(genes["subtelomeric"].sum()),
        "frac_subtelomeric": float(genes["subtelomeric"].mean()),
        "genes_per_family": genes_per_family,
        "largest_family": _top(family_summary, "n_genes"),
        "densest_family_per_mb": _top(family_summary, "genes_per_mb_genome"),
        "most_clustered_family": _top(
            array_summary[array_summary["n_genes"] >= 3], "frac_clustered"
        ),
        "telomere_bias_table": tel_table,
        "families_with_telomere_bias": [
            {
                "gene_family": r["gene_family"],
                "direction": r["direction"],
                "p_toward_telomere": r["p_toward_telomere"],
                "p_toward_interior": r["p_toward_interior"],
            }
            for _, r in telomere_bias[telomere_bias["significant"]].iterrows()
        ],
        "duplication_modes": dup_table,
        "ripley_clustering": ripley_summary,
        "chromosome_enrichment_hits": enrich_hits,
        "strand_biased_families": strand_hits,
        "centromere_bias_table": cen_table,
        "arm_bias_table": arm_table,
        "n_hotspots": int(len(hotspots)),
        "hotspot_regions": hotspot_regions,
        "family_proximity_order": list(prox_order or []),
        "family_proximity_clusters": prox_clusters_json,
        "alpha": alpha,
        "fdr_method": "benjamini-hochberg",
    }
    if coloc is not None and not coloc.empty:
        summary["colocalization_table"] = [
            {
                "family_a": r["family_a"],
                "family_b": r["family_b"],
                "observed_close_pairs": int(r["observed_close_pairs"]),
                "null_mean": round4(r["null_mean"]),
                "p_value": round4(r["p_value"]),
                "q_value": round4(r.get("q_value", np.nan)),
                "significant": bool(r["significant"]),
                "significant_fdr": bool(r.get("significant_fdr", False)),
            }
            for _, r in coloc.iterrows()
        ]
        summary["significant_colocalizations"] = [
            {"family_a": r["family_a"], "family_b": r["family_b"], "p_value": r["p_value"]}
            for _, r in coloc[coloc["significant"]].iterrows()
        ]
    return summary


def _top(df, col):
    if df is None or df.empty or col not in df:
        return None
    row = df.loc[df[col].idxmax()]
    return {"gene_family": row["gene_family"], col: float(row[col])}
