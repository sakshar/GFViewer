"""Command-line interface for GFViewer.

Backward compatible with the original flags (``-d -g -o -c -l -t -p -cen``)
while adding: multiple input files, BED/GFF3/GTF, unlimited-ish families with a
soft cap, full styling, multi-format export, style files, and analytics.
"""

import argparse
import os
import sys

from gfviewer import __version__
from gfviewer.errors import GFViewerError
from gfviewer.palette import color_guide
from gfviewer.pipeline import run
from gfviewer.style import StyleConfig


def build_parser():
    p = argparse.ArgumentParser(
        prog="gfviewer",
        description="Visualize the localization of multigene families across chromosomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--version", action="version", version="GFViewer " + __version__)
    p.add_argument(
        "--color-guide", action="store_true",
        help="print the built-in colour palette and exit",
    )

    io = p.add_argument_group("input / output")
    io.add_argument(
        "-d", "--data", "--data_file", dest="data", nargs="+", metavar="FILE",
        help="one or more annotation files (.xlsx/.csv/.tsv/.bed/.gff3/.gtf). "
        "With several files, each file's base name is used as the gene family "
        "unless the file already carries family information.",
    )
    io.add_argument(
        "-g", "--genome", "--genome_file", dest="genome", metavar="FILE",
        help="genome FASTA, .fai, chrom.sizes, or '<chrom>,<length>' table",
    )
    io.add_argument(
        "-o", "--output", "--output_directory", dest="output", metavar="DIR",
        help="output directory",
    )
    io.add_argument("-c", "--colors", "--color_map_file", dest="colors", metavar="FILE",
                    help="colour-map file (<family>,<code|name|#hex> or <family>,r,g,b)")
    io.add_argument("-m", "--mapping", metavar="FILE",
                    help="gene_id,gene_family mapping (single annotation file only)")
    io.add_argument("--family-attr", metavar="KEY",
                    help="GFF3/GTF attribute holding the gene family")
    io.add_argument("--id-attr", metavar="KEY",
                    help="GFF3/GTF attribute holding the gene id")
    io.add_argument("--gff-types", metavar="T", nargs="+",
                    help="GFF3/GTF feature types to keep (default: gene-like)")
    io.add_argument("--basename", default="gfviewer", help="output file base name")
    io.add_argument(
        "-f", "--format", dest="formats", nargs="+", default=None,
        metavar="FMT", help="export formats: pdf svg png jpg tif eps",
    )
    io.add_argument("--on-unknown-chrom", choices=["error", "drop", "unplaced"],
                    default="error",
                    help="features on a sequence missing from the genome file: "
                    "error out, drop them, or keep them as unplaced contigs")
    io.add_argument("--coord-bounds", choices=["clip", "error", "keep"], default="clip")

    fam = p.add_argument_group("gene families")
    fam.add_argument("--collapse-rare", action="store_true",
                     help="merge infrequent families into 'Other' when above the cap")
    fam.add_argument("--keep", nargs="+", metavar="FAMILY", default=None,
                     help="families to keep out of the 'Other' bucket")

    an = p.add_argument_group("analytics")
    an.add_argument("--analytics", action="store_true",
                    help="also compute per-family statistics and write CSV/JSON")
    an.add_argument("--analytics-only", action="store_true",
                    help="compute analytics and skip rendering")
    an.add_argument("--subtelomere-fraction", type=float, default=0.10)
    an.add_argument("--subtelomere-bp", type=int, default=None)
    an.add_argument("--cluster-gap", type=int, default=50000,
                    help="max gap (bp) between paralogues in a tandem array")
    an.add_argument("--proximal-window", type=int, default=None,
                    help="max gap (bp) for the 'proximal' duplication class "
                    "(default: 5x --cluster-gap)")
    an.add_argument("--ripley-scales", type=int, nargs="+", default=None,
                    metavar="BP", help="scales for the 1-D Ripley K/L test")
    an.add_argument("--hotspot-window", type=int, default=100000,
                    help="window (bp) for the multigene-family hotspot scan")
    an.add_argument("--hotspot-step", type=int, default=None,
                    help="step (bp) for the hotspot scan (default: window / 2)")
    an.add_argument("--proximity-clusters", type=int, default=None,
                    help="number of family clusters to cut the proximity "
                    "dendrogram into (default: auto)")
    an.add_argument("--permutations", type=int, default=1000)
    an.add_argument("--seed", type=int, default=0)
    an.add_argument("--colocalization", action="store_true")

    st = p.add_argument_group("style (overrides --style file / defaults)")
    st.add_argument("--style", metavar="FILE", help="YAML/JSON style file")
    st.add_argument("--save-style", metavar="FILE", help="write the resolved style and exit")
    st.add_argument("--orientation", choices=["horizontal", "vertical"],
                    help="draw chromosomes horizontally (default) or vertically")
    st.add_argument("-t", "--telomere", "--telomere_length", dest="telomere_length", type=int)
    st.add_argument("-p", "--per-page", "--number_of_chromosomes_per_page",
                    dest="per_page", type=int,
                    help="chromosomes per page (rows x columns)")
    st.add_argument("--columns", type=int,
                    help="chromosomes side by side (0 = auto)")
    st.add_argument("--row-height", type=float, dest="row_height",
                    help="cm of space per chromosome across its thickness")
    st.add_argument("--length-cm", type=float, dest="length_cm",
                    help="length-axis budget in cm (vertical orientation)")
    st.add_argument("--single-page", action="store_true",
                    help="put every chromosome on one page")
    st.add_argument("--show-unplaced", action="store_true",
                    help="also draw unplaced / stray contigs, not just chromosomes")
    st.add_argument("--title")
    st.add_argument("--subtitle")
    st.add_argument("--no-titles", action="store_true",
                    help="omit titles from the exported chromosome figure and "
                    "the analytics-chart image files")
    st.add_argument("--dpi", type=int)
    st.add_argument("--page-size", nargs=2, type=float, metavar=("W_CM", "H_CM"))
    st.add_argument("--background")
    st.add_argument("--tick-style", choices=["line", "lollipop", "triangle", "box", "arrow"])
    st.add_argument("--no-split-strand", action="store_true")
    st.add_argument("-cen", "--centromeres", dest="centromeres", action="store_true",
                    help="draw centromeres")
    st.add_argument("--label-mode", choices=["none", "family", "gene_id"])
    st.add_argument("--label-size", type=float)
    st.add_argument("--font", help="font family for labels/legend/chrom names")
    st.add_argument("-l", "--legend-location", "--legend_location", dest="legend_location",
                    choices=["upper left", "upper center", "upper right",
                             "center left", "center", "center right",
                             "lower left", "lower center", "lower right",
                             "outside right", "outside bottom", "none"])
    st.add_argument("--legend-columns", type=int)
    st.add_argument("--legend-size", type=float)
    st.add_argument("--legend-title")
    st.add_argument("--legend-frame", action="store_true")
    st.add_argument("--legend-separate-page", action="store_true")
    st.add_argument("--no-legend", action="store_true")
    return p


def _resolve_style(args):
    style = StyleConfig.load(args.style) if args.style else StyleConfig()

    if args.orientation:
        style.orientation = args.orientation
    if args.telomere_length is not None:
        style.telomere_length = args.telomere_length
    if args.columns is not None:
        style.chromosomes_per_row = args.columns
    if args.row_height is not None:
        style.row_height_cm = args.row_height
    if args.length_cm is not None:
        style.length_cm = args.length_cm
    if args.per_page is not None:
        cols = style.chromosomes_per_row or 1
        style.rows_per_page = max(1, int(round(args.per_page / cols)))
    if args.title is not None:
        style.title = args.title
    if args.subtitle is not None:
        style.subtitle = args.subtitle
    if args.no_titles:
        style.show_titles = False
    if args.dpi is not None:
        style.dpi = args.dpi
    if args.page_size:
        style.page_width_cm, style.page_height_cm = args.page_size
    if args.background:
        style.background = args.background
    if args.tick_style:
        style.tick_style = args.tick_style
    if args.no_split_strand:
        style.split_by_strand = False
    if args.centromeres:
        style.show_centromeres = True
    if args.show_unplaced:
        style.show_unplaced = True
    if args.label_mode:
        style.label_mode = args.label_mode
    if args.label_size is not None:
        style.label_font_size = args.label_size
    if args.font:
        style.label_font_family = args.font
        style.legend_font_family = args.font
        style.chrom_label_font_family = args.font
        style.title_font_family = args.font
    if args.legend_location:
        style.legend_location = args.legend_location
        if args.legend_location == "none":
            style.legend_show = False
    if args.legend_columns is not None:
        style.legend_columns = args.legend_columns
    if args.legend_size is not None:
        style.legend_font_size = args.legend_size
    if args.legend_title is not None:
        style.legend_title = args.legend_title
    if args.legend_frame:
        style.legend_frame = True
    if args.legend_separate_page:
        style.legend_separate_page = True
    if args.no_legend:
        style.legend_show = False
    if args.formats:
        style.formats = list(args.formats)
    return style.validate()


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.color_guide:
        print("GFViewer colour palette (index: name):\n")
        for idx, name, rgb in color_guide():
            print("  {:>2}: {:<9} rgb{}".format(idx, name, tuple(round(v, 2) for v in rgb)))
        print("\nBeyond 20, colours are generated automatically. "
              "Grey is reserved for centromeres.")
        return 0

    missing = [n for n in ("data", "genome", "output") if not getattr(args, n)]
    if missing:
        print("error: missing required argument(s): --" + ", --".join(missing), file=sys.stderr)
        return 2

    try:
        style = _resolve_style(args)
        if args.save_style:
            style.save(args.save_style)
            print("Wrote style to", args.save_style)
            return 0

        akw = dict(
            subtelomere_fraction=args.subtelomere_fraction,
            subtelomere_bp=args.subtelomere_bp,
            cluster_gap=args.cluster_gap,
            proximal_window=args.proximal_window,
            ripley_scales=args.ripley_scales,
            hotspot_window=args.hotspot_window,
            hotspot_step=args.hotspot_step,
            n_proximity_clusters=args.proximity_clusters,
            n_permutations=args.permutations,
            seed=args.seed,
            colocalization=args.colocalization,
        )
        result = run(
            data_paths=args.data,
            genome_path=args.genome,
            outdir=args.output,
            style=style,
            color_file=args.colors,
            mapping_file=args.mapping,
            family_attr=args.family_attr,
            id_attr=args.id_attr,
            gff_types=args.gff_types,
            collapse_rare=args.collapse_rare,
            keep_families=args.keep,
            basename=args.basename,
            with_analytics=args.analytics or args.analytics_only,
            analytics_kwargs=akw,
            on_unknown_chrom=args.on_unknown_chrom,
            coord_bounds=args.coord_bounds,
            render_single_page=args.single_page,
        )
    except GFViewerError as exc:
        print("error:", exc, file=sys.stderr)
        return 1

    for w in result.warnings:
        print("warning:", w, file=sys.stderr)
    print("Gene families ({}):".format(len(result.families)), ", ".join(result.families))
    if not args.analytics_only:
        for fmt, paths in result.outputs.items():
            for path in paths:
                print("  wrote", os.path.relpath(path))
    for path in result.analytics_files:
        print("  wrote", os.path.relpath(path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
