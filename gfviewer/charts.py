"""Stand-alone analytics charts.

Kept separate from :mod:`gfviewer.render` (the chromosome ideogram) so they can
be produced independently by the CLI (``--analytics``) and by the web portal's
on-demand download endpoint, in any of the usual figure formats.

Figures: genes-per-family bar chart, positional density profile
("metachromosome" plot), Ripley's L multi-scale clustering curve, and the
family-proximity heat map.
"""

import os

import matplotlib
import numpy as np

if os.environ.get("GFVIEWER_MPL_BACKEND"):
    matplotlib.use(os.environ["GFVIEWER_MPL_BACKEND"])
else:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

CM = 1 / 2.54
_RASTER = {"png", "jpg", "jpeg", "tif", "tiff"}
_UNPLACED_HATCH = "////"
_UNPLACED_FACE = "#d8d8d8"
_UNPLACED_EDGE = "#9a9a9a"
_DEFAULT_BAR = (0.48, 0.12, 0.17)


def genes_per_family_figure(rows, color_map=None, title="Genes per family"):
    """Return a matplotlib ``Figure`` with a horizontal stacked bar per family.

    *rows* is the ``summary["genes_per_family"]`` list -- dicts with keys
    ``gene_family``, ``on_chromosomes``, ``on_unplaced``, ``total`` -- ordered
    however the caller wants them stacked top-to-bottom.
    """
    rows = list(rows)
    if not rows:
        raise ValueError("No per-family gene counts to chart.")
    color_map = color_map or {}

    labels = [r["gene_family"] for r in rows]
    on_chrom = [int(r.get("on_chromosomes", 0)) for r in rows]
    on_unpl = [int(r.get("on_unplaced", 0)) for r in rows]
    any_unpl = any(on_unpl)

    height = max(4.0, 0.42 * len(rows) + 1.6) + (0.9 if any_unpl else 0.0)
    fig, ax = plt.subplots(figsize=(18 * CM, height * CM), dpi=200)
    y = list(range(len(rows)))[::-1]  # first row at the top

    for yi, fam, nc, nu in zip(y, labels, on_chrom, on_unpl):
        col = color_map.get(fam, _DEFAULT_BAR)
        ax.barh(yi, nc, color=col, edgecolor="none", height=0.62)
        if nu:
            ax.barh(yi, nu, left=nc, height=0.62, facecolor=_UNPLACED_FACE,
                    edgecolor=_UNPLACED_EDGE, hatch=_UNPLACED_HATCH, linewidth=0.5)
        ax.text(nc + nu, yi, "  " + str(nc + nu), va="center", ha="left", fontsize=8)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Number of genes", fontsize=9)
    ax.margins(x=0.08, y=0.01)
    ax.spines[["top", "right"]].set_visible(False)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)

    if any_unpl:
        reserve = min(0.22, 1.4 / height)
        fig.tight_layout(rect=(0, reserve, 1, 1))
        fig.legend(
            handles=[
                Patch(facecolor=_DEFAULT_BAR, label="on chromosomes"),
                Patch(facecolor=_UNPLACED_FACE, edgecolor=_UNPLACED_EDGE,
                      hatch=_UNPLACED_HATCH, label="on unplaced / stray contigs"),
            ],
            loc="lower center", ncol=2, fontsize=8, frameon=False,
        )
    else:
        fig.tight_layout()
    return fig


def save_figure(fig, path, dpi=200):
    """Write *fig* to *path*; format inferred from the extension."""
    fmt = os.path.splitext(path)[1].lstrip(".").lower()
    kw = {"facecolor": "white", "bbox_inches": "tight"}
    if fmt in _RASTER:
        kw["dpi"] = dpi
    if fmt in ("jpg", "jpeg"):
        kw["pil_kwargs"] = {"quality": 92}
    fig.savefig(path, **kw)
    plt.close(fig)
    return path


def _write(make_fig, outdir, basename, formats, dpi):
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for fmt in formats:
        ext = "jpg" if fmt.lower() == "jpeg" else fmt.lower()
        paths.append(save_figure(make_fig(), os.path.join(outdir, "{}.{}".format(basename, ext)), dpi))
    return paths


def write_genes_per_family(rows, outdir, basename="analytics_genes_per_family",
                           formats=("pdf",), color_map=None, dpi=200, titles=True):
    """Render + save the chart in every requested format.  Returns the paths.

    *titles* keeps (True) or drops (False) the built-in chart heading in the
    saved image files.
    """
    ttl = "Genes per family" if titles else None
    return _write(lambda: genes_per_family_figure(rows, color_map=color_map, title=ttl),
                  outdir, basename, formats, dpi)


# --------------------------------------------------------------------------- #
# positional density profile ("metachromosome" plot)
# --------------------------------------------------------------------------- #
def positional_profile_figure(df, color_map=None, title="Positional density profile"):
    """Line plot of each family's gene density along the normalised chromosome
    (0 = start telomere, 1 = end telomere).  *df* is the ``positional_profile``
    table (columns ``gene_family``, ``bin_center``, ``frac_of_family``)."""
    color_map = color_map or {}
    fams = [f for f in dict.fromkeys(df["gene_family"]) if f != "(all families)"]
    nbins = int(df["bin"].max()) + 1 if "bin" in df.columns else 20

    fig, ax = plt.subplots(figsize=(20 * CM, 11 * CM), dpi=200)
    for fam in fams:
        g = df[df["gene_family"] == fam].sort_values("bin_center")
        ax.plot(g["bin_center"], g["frac_of_family"], lw=1.3,
                color=color_map.get(fam), label=fam)
    allg = df[df["gene_family"] == "(all families)"].sort_values("bin_center")
    if not allg.empty:
        ax.plot(allg["bin_center"], allg["frac_of_family"], lw=2.4, ls="--",
                color="#333333", label="all families", zorder=6)
    ax.axhline(1.0 / nbins, color="#999999", lw=0.8, ls=":")
    ax.axvline(0.5, color="#cccccc", lw=0.8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("position along chromosome  (0 = start telomere, 0.5 = centre, 1 = end telomere)",
                  fontsize=9)
    ax.set_ylabel("fraction of the family's genes", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    ax.legend(fontsize=7, ncol=max(1, min(6, (len(fams) + 1) // 2)),
              frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    return fig


def write_positional_profile(df, outdir, basename="analytics_positional_profile",
                             formats=("pdf",), color_map=None, dpi=200, titles=True):
    ttl = "Positional density profile" if titles else None
    return _write(lambda: positional_profile_figure(df, color_map=color_map, title=ttl),
                  outdir, basename, formats, dpi)


# --------------------------------------------------------------------------- #
# Ripley's L multi-scale clustering
# --------------------------------------------------------------------------- #
def ripley_figure(df, color_map=None, title="Multi-scale clustering (1-D Ripley's L)"):
    """L(t) - t versus scale t, one line per family.  Filled markers flag scales
    that are significant after FDR correction.  *df* is the ``ripley`` table."""
    color_map = color_map or {}
    fig, ax = plt.subplots(figsize=(20 * CM, 11 * CM), dpi=200)
    for fam, g in df.groupby("gene_family"):
        g = g.sort_values("scale_bp")
        col = color_map.get(fam)
        ax.plot(g["scale_bp"], g["L_minus_t_obs"], lw=1.2, marker="o", ms=3,
                color=col, label=fam)
        if "significant_fdr" in g.columns:
            sig = g[g["significant_fdr"]]
            if not sig.empty:
                ax.scatter(sig["scale_bp"], sig["L_minus_t_obs"], s=34, color=col,
                           edgecolor="black", linewidths=0.4, zorder=6)
    ax.axhline(0.0, color="#666666", lw=0.9)
    ax.set_xscale("log")
    ax.set_xlabel("spatial scale  t  (bp, log axis)", fontsize=9)
    ax.set_ylabel("L(t) − t   (bp;  > 0 ⇒ clustered,  < 0 ⇒ over-dispersed)",
                  fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    n = df["gene_family"].nunique()
    ax.legend(fontsize=7, ncol=max(1, min(6, (n + 1) // 2)), frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    return fig


def write_ripley(df, outdir, basename="analytics_ripley",
                 formats=("pdf",), color_map=None, dpi=200, titles=True):
    ttl = "Multi-scale clustering (1-D Ripley's L)" if titles else None
    return _write(lambda: ripley_figure(df, color_map=color_map, title=ttl),
                  outdir, basename, formats, dpi)


# --------------------------------------------------------------------------- #
# family-proximity heat map
# --------------------------------------------------------------------------- #
def family_proximity_figure(mat, order=None, clusters=None,
                            title="Family proximity (mean cross-nearest-neighbour distance)"):
    """Heat map of the family x family distance matrix, families reordered by the
    average-linkage clustering; darker = physically closer.  *mat* is a square
    DataFrame indexed and columned by family name."""
    fams = [f for f in (order or list(mat.index)) if f in mat.index]
    if len(fams) < 2:
        fig, ax = plt.subplots(figsize=(8 * CM, 6 * CM))
        ax.text(0.5, 0.5, "needs ≥ 2 families", ha="center", va="center")
        ax.axis("off")
        return fig
    M = mat.loc[fams, fams].to_numpy(dtype=float)
    finite = M[np.isfinite(M)]
    vmax = float(np.percentile(finite, 95)) if finite.size else 1.0
    disp = np.where(np.isfinite(M), M, vmax) / 1e3          # -> kb

    n = len(fams)
    side = max(9.0, 0.46 * n + 4.0)
    fig, ax = plt.subplots(figsize=(side * CM, side * CM), dpi=200)
    im = ax.imshow(disp, cmap="viridis_r", vmin=0, vmax=vmax / 1e3, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(fams, rotation=90, fontsize=7)
    ax.set_yticklabels(fams, fontsize=7)

    if clusters:
        pos = {f: i for i, f in enumerate(fams)}
        for grp in clusters:
            idx = sorted(pos[f] for f in grp.get("families", grp) if f in pos)
            if not idx:
                continue
            a, b = idx[0] - 0.5, idx[-1] + 0.5
            ax.add_patch(Rectangle((a, a), b - a, b - a, fill=False,
                                   edgecolor="crimson", lw=1.6))
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("mean cross-nearest-neighbour distance (kb)", fontsize=8)
    if title:
        ax.set_title(title, fontsize=11, fontweight="bold", pad=10)
    fig.tight_layout()
    return fig


def write_family_proximity(mat, outdir, basename="analytics_family_proximity",
                           order=None, clusters=None, formats=("pdf",), dpi=200,
                           titles=True):
    ttl = ("Family proximity (mean cross-nearest-neighbour distance)"
           if titles else None)
    return _write(lambda: family_proximity_figure(mat, order=order,
                                                  clusters=clusters, title=ttl),
                  outdir, basename, formats, dpi)
