"""High-level orchestration shared by the CLI and the web app.

``run`` performs the whole flow -- read genome, load + validate features, build
the palette, optionally compute analytics, render, and save -- and returns a
:class:`PipelineResult` describing everything produced.
"""

import os
from dataclasses import dataclass, field

from gfviewer import analytics as _analytics
from gfviewer import render as _render
from gfviewer.genome import load_genome
from gfviewer.io import load_features
from gfviewer.palette import build_palette
from gfviewer.style import StyleConfig


@dataclass
class PipelineResult:
    genome: dict
    features: object                 # pandas DataFrame
    color_map: dict
    style: StyleConfig
    warnings: list = field(default_factory=list)
    outputs: dict = field(default_factory=dict)     # {fmt: [paths]}
    analytics: object = None                         # AnalyticsResult or None
    analytics_files: list = field(default_factory=list)
    figures: list = field(default_factory=list)
    collapsed_families: set = field(default_factory=set)

    @property
    def families(self):
        return [f for f in self.color_map if f != "Other"]


def run(
    data_paths,
    genome_path,
    outdir,
    style=None,
    color_file=None,
    mapping_file=None,
    family_attr=None,
    id_attr=None,
    gff_types=None,
    collapse_rare=False,
    keep_families=None,
    formats=None,
    basename="gfviewer",
    with_analytics=False,
    analytics_kwargs=None,
    on_unknown_chrom="error",
    coord_bounds="clip",
    keep_figures=False,
    render_single_page=False,
):
    """Run the full GFViewer pipeline.  See :class:`PipelineResult`."""
    style = style or StyleConfig()
    if formats:
        style.formats = list(formats)
    style.validate()
    os.makedirs(outdir, exist_ok=True)
    warnings = []

    genome = load_genome(genome_path)

    features, io_warn = load_features(
        data_paths, genome,
        family_attr=family_attr, id_attr=id_attr, gff_types=gff_types,
        mapping_file=mapping_file, on_unknown_chrom=on_unknown_chrom,
        coord_bounds=coord_bounds,
    )
    warnings += io_warn

    genes = features[features["kind"] == "gene"]
    fam_order = list(dict.fromkeys(genes["gene_family"].tolist()))
    fam_counts = genes["gene_family"].value_counts().to_dict()

    color_map, pal_warn, collapsed = build_palette(
        fam_order, color_file=color_file, collapse_rare=collapse_rare,
        family_counts=fam_counts, keep=keep_families,
    )
    warnings += pal_warn

    if collapsed:
        features = features.copy()
        mask = features["gene_family"].isin(collapsed)
        features.loc[mask, "gene_family"] = "Other"

    # analytics operate on the (possibly collapsed) feature set
    result = PipelineResult(
        genome=genome, features=features, color_map=color_map, style=style,
        warnings=warnings, collapsed_families=collapsed,
    )

    if with_analytics:
        akw = dict(analytics_kwargs or {})
        result.analytics = _analytics.compute_analytics(features, genome, **akw)
        result.analytics_files = result.analytics.write(
            outdir, chart_formats=style.formats, color_map=color_map, dpi=style.dpi,
            chart_titles=style.show_titles,
        )

    if render_single_page:
        style.rows_per_page = max(len(genome), 1)

    figs = _render.render(features, genome, style, color_map)
    result.outputs = _render.save_figures(figs, outdir, basename, style.formats, dpi=style.dpi)
    if keep_figures:
        result.figures = figs
    else:
        _render.close_all(figs)

    return result
