"""matplotlib ideogram renderer for GFViewer.

Chromosomes are drawn as round-capped bars -- **horizontal** (stacked in rows)
or **vertical** (placed side by side) depending on ``style.orientation``.  Genes
appear as coloured marks along each bar (``+`` strand on one side of the axis,
``-`` on the other, when ``split_by_strand`` is on).  When centromeres are shown
with ``centromere_style="constriction"`` the bar is drawn as a proper
cytogenetic ideogram: a p-arm and a q-arm separated by a dark primary
constriction.

Everything is vector, so one render can be written to PDF / SVG / PNG / JPG /
TIFF / EPS, and the SVG carries stable ``id`` attributes (``gfv-fam-<family>``,
``gfv-legend``, ``gfv-chrom-<name>``, ``gfv-label-<gene>``, ``gfv-title``) for
the interactive web editor.

All geometry is in **centimetres** with ``set_aspect('equal')``.

Public API
----------
``render(features, genome, style, color_map)``          -> list[Figure]
``save_figures(figs, outdir, basename, formats, dpi)``  -> {fmt: [paths]}
``figure_to_svg(fig)``                                  -> str
"""

import io as _io
import math
import os
import re

import matplotlib

if os.environ.get("GFVIEWER_MPL_BACKEND"):
    matplotlib.use(os.environ["GFVIEWER_MPL_BACKEND"])
else:  # this tool only ever writes files
    matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

CM = 1 / 2.54
PT_TO_CM = 2.54 / 72.0
_RASTER = {"png", "jpg", "jpeg", "tif", "tiff"}
_BOTTOM_LEGEND = {"lower left", "lower center", "lower right", "outside bottom"}
_SHELF_GAP = 1.2


def _slug(text):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(text))


def _label_gutter_cm(style, max_chars):
    """Perpendicular distance (cm) a gene label extends *beyond the tick tip*,
    on one side of a chromosome bar, for the current label mode.

    This is what neighbouring chromosomes have to be spaced apart by (on top of
    the bar thickness and the tick) so labels coming off adjacent chromosomes do
    not overlap.  ``max_chars`` is the longest label string that will be drawn.
    """
    if style.label_mode == "none" or max_chars <= 0:
        return 0.0
    char_cm = style.label_font_size * PT_TO_CM
    horiz = style.orientation == "horizontal"
    rotated = style.label_mode == "gene_id" and horiz
    if rotated:
        # text is turned 90 deg, so the whole string sticks out perpendicular
        reach_max = 0.32 + 2 * char_cm            # n_bands = 3, band_step = char_cm
        text_cm = max_chars * char_cm * 0.62
        clear = 0.15
    elif horiz:
        # upright family labels: roughly one text line tall
        reach_max = 0.32 + (char_cm * 1.5)        # n_bands = 2, band_step = 1.5*char_cm
        text_cm = char_cm * 1.3
        clear = 0.15
    else:
        # vertical orientation: upright text placed left/right of the bar, so the
        # whole string width sticks out perpendicular towards the next column.
        # DejaVu Sans alphanumerics average ~0.66 em; leave real clearance so the
        # two facing label stacks of adjacent columns never touch.
        reach_max = 0.32 + (char_cm * 1.5)
        text_cm = max_chars * char_cm * 0.68
        clear = 0.4
    return reach_max + text_cm + clear


# --------------------------------------------------------------------------- #
# page geometry
# --------------------------------------------------------------------------- #
class _Layout:
    """Orientation-aware page geometry, in centimetres."""

    def __init__(self, style, n_on_page, has_title=None, label_gutter=0.0):
        self.orient = style.orientation
        want_title = style.title if has_title is None else has_title

        cpr = int(style.chromosomes_per_row)
        if cpr <= 0:
            cpr = 1 if self.orient == "horizontal" else min(max(n_on_page, 1), 12)
        self.cpr = max(1, cpr)
        self.lines = max(1, math.ceil(n_on_page / self.cpr))

        self.gap = min(0.85, max(0.0, style.chromosome_gap))

        # ---- chromosome pitch, auto-widened for labels --------------------- #
        # ``across`` is the centre-to-centre spacing of neighbouring
        # chromosomes (row pitch when horizontal, column pitch when vertical).
        # When gene / family labels are drawn they stick out perpendicular to
        # the bar, so the pitch must cover, on *each* facing side, half the bar
        # plus the tick plus the label gutter -- otherwise labels from adjacent
        # chromosomes overlap.  Because the drawing block is inset by half a
        # pitch at each end, sizing the pitch this way also keeps the outermost
        # labels on the canvas.
        base_across = max(1.2, float(style.row_height_cm))
        self.across = base_across
        want_labels = style.label_mode != "none" and label_gutter > 0.0
        if want_labels:
            prov_bh = min(min(style.body_width, 0.35) * base_across, 1.1)
            prov_tick = max(0.15, style.tick_length) * prov_bh
            need = 2.0 * (prov_bh / 2.0 + prov_tick + label_gutter) + 0.3
            # cap the automatic growth so a few very long IDs can't explode the page
            self.across = min(max(base_across, need), base_across + 12.0)

        self.body_h = min(min(style.body_width, 0.35) * self.across, 1.1)
        self.tick_len = max(0.15, style.tick_length) * self.body_h

        # whatever the widened pitch still can't absorb spills into the outer
        # page margin (below) rather than off the edge of the figure
        gutter_room = self.across / 2.0 - (self.body_h / 2.0 + self.tick_len)
        self._label_overflow = max(0.0, label_gutter - gutter_room) if want_labels else 0.0

        self.m_top = 2.6 if want_title else 1.0
        bottom_legend = (
            style.legend_show and style.legend_location in _BOTTOM_LEGEND
            and not style.legend_separate_page
        )

        # a bottom/centred legend needs a minimum canvas width to fit
        legend_min_w = 22.0 if (
            style.legend_show and not style.legend_separate_page
            and (bottom_legend or style.legend_location in ("upper center", "center"))
        ) else 0.0

        if self.orient == "horizontal":
            self.m_left = 2.6
            self.m_right = 0.5
            self.m_bottom = (3.4 if bottom_legend else 0.7)
            # labels stick out above the top row / below the bottom row
            self.m_top += self._label_overflow
            self.m_bottom += self._label_overflow
            self.W = float(style.page_width_cm)
            block = self.across * self.lines
            if style.auto_page_height:
                self.H = self.m_top + block + self.m_bottom
            else:
                self.H = float(style.page_height_cm)
            self.along_budget = (self.W - self.m_left - self.m_right) / self.cpr
            self.x_pad = self.m_left
        else:  # vertical
            self.m_left = 1.0
            self.m_right = 0.5
            # labels stick out left of the first column / right of the last one
            self.m_left += self._label_overflow
            self.m_right += self._label_overflow
            self.m_bottom = (3.4 if bottom_legend else 1.4)   # room for chrom names
            self.along_budget = max(4.0, float(style.length_cm))
            self.H = (
                self.m_top + self.lines * self.along_budget
                + (self.lines - 1) * _SHELF_GAP + self.m_bottom
            )
            block = self.across * self.cpr
            content_w = self.m_left + block + self.m_right
            self.W = max(content_w, legend_min_w)
            self.x_pad = (self.W - block) / 2.0   # centre the chromosome columns

    # ---- per-chromosome coordinate frame -------------------------------- #
    def placer(self, index):
        line, col = divmod(index, self.cpr)
        if self.orient == "horizontal":
            usable_w = self.W - self.m_left - self.m_right
            col_w = usable_w / self.cpr
            x0 = self.m_left + col * col_w
            y_top = self.H - self.m_top
            yc = y_top - (line + 0.5) * self.across
            return _Placer("horizontal", (x0, yc), (1.0, 0.0), (0.0, 1.0),
                           col_w * (1.0 - self.gap))
        # vertical
        xc = self.x_pad + (col + 0.5) * self.across
        shelf_top = self.H - self.m_top - line * (self.along_budget + _SHELF_GAP)
        return _Placer("vertical", (xc, shelf_top), (0.0, -1.0), (1.0, 0.0),
                       self.along_budget * (1.0 - self.gap))


class _Placer:
    """Maps chromosome-local ``(l, t)`` -> page ``(x, y)`` in cm.

    ``l`` runs from 0 at the chromosome's first base along the length axis;
    ``t`` is the signed offset across the bar (``+`` = the ``+``-strand side).
    """

    def __init__(self, orient, origin, u_l, u_t, along_max):
        self.orient = orient
        self.ox, self.oy = origin
        self.ulx, self.uly = u_l
        self.utx, self.uty = u_t
        self.along_max = along_max

    def pt(self, l, t):
        return (self.ox + l * self.ulx + t * self.utx,
                self.oy + l * self.uly + t * self.uty)

    def seg(self, l0, t0, l1, t1):
        return (self.pt(l0, t0), self.pt(l1, t1))

    def poly(self, pts_lt):
        return [self.pt(l, t) for l, t in pts_lt]


# --------------------------------------------------------------------------- #
def render(features, genome, style, color_map):
    """Render *features* to a list of figures (one per page)."""
    style.validate()

    # optional family subset -- drives marks, labels *and* the legend, since they
    # are all keyed off ``color_map``; ignore unknown ids, never blank the figure
    if style.only_families:
        wanted = set(style.only_families)
        sub = {k: v for k, v in color_map.items() if k in wanted}
        if sub:
            color_map = sub

    if style.show_unplaced:
        names = list(genome)
    else:
        from gfviewer.genome import chromosome_id_list

        names = chromosome_id_list(genome) or list(genome)

    # optional user subset -- keep genome order; ignore unknown / stale ids and
    # never fall through to an empty figure
    if style.only_chromosomes:
        wanted = set(style.only_chromosomes)
        subset = [n for n in names if n in wanted]
        if subset:
            names = subset

    cpr = int(style.chromosomes_per_row)
    if cpr <= 0:
        cpr = 1 if style.orientation == "horizontal" else min(len(names), 12)
    per_page = max(1, max(1, cpr) * style.rows_per_page)
    n_pages = math.ceil(len(names) / per_page)

    genes = features[features["kind"] == "gene"]
    cents = features[features["kind"] == "centromere"]
    ref_bp = (max(genome.values()) + 2 * style.telomere_length) if style.shared_scale else None

    # widen the spacing between chromosomes so labels from neighbouring
    # chromosomes cannot overlap (and cannot run off the page)
    label_gutter = 0.0
    if style.label_mode == "family":
        chars = max((len(str(f)) for f in color_map), default=0)
        label_gutter = _label_gutter_cm(style, min(chars, 28))
    elif style.label_mode == "gene_id":
        labelled = genes[genes["gene_family"].isin(color_map)]
        chars = max((len(str(x)) for x in labelled["gene_id"]), default=0)
        label_gutter = _label_gutter_cm(style, min(chars, 28))

    figs = []
    for page in range(n_pages):
        page_names = names[page * per_page:(page + 1) * per_page]
        has_title = bool(style.title) and page == 0 and getattr(style, "show_titles", True)
        fig, ax, lay = _new_page(style, len(page_names), has_title, label_gutter)
        if has_title:
            _title(ax, lay, style)

        fam_segs, fam_tips, labels = {}, {}, []
        for i, chrom in enumerate(page_names):
            _draw_chromosome(
                ax, lay, i, chrom, int(genome[chrom]), style, ref_bp,
                genes[genes["chromosome"] == chrom],
                cents[cents["chromosome"] == chrom],
                color_map, fam_segs, fam_tips, labels,
            )
        _emit_family_marks(ax, style, color_map, fam_segs, fam_tips)
        if style.label_mode != "none":
            _emit_labels(ax, style, labels, lay)
        if (
            style.legend_show and not style.legend_separate_page
            and (style.legend_per_page or page == 0)
        ):
            _draw_legend(ax, lay, style, color_map)
        figs.append(fig)

    if style.legend_separate_page and style.legend_show:
        figs.append(_legend_page(style, color_map))
    return figs


# --------------------------------------------------------------------------- #
def _new_page(style, n_on_page, has_title=None, label_gutter=0.0):
    lay = _Layout(style, n_on_page, has_title, label_gutter)
    fig = plt.figure(figsize=(lay.W * CM, lay.H * CM), dpi=style.dpi)
    fig.patch.set_facecolor(style.background)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, lay.W)
    ax.set_ylim(0, lay.H)
    ax.set_aspect("equal")
    ax.set_facecolor(style.background)
    ax.axis("off")
    return fig, ax, lay


def _title(ax, lay, style):
    t = ax.text(
        lay.W / 2.0, lay.H - 1.1, style.title, ha="center", va="center",
        fontsize=style.title_font_size, fontfamily=style.title_font_family,
        color=style.title_color, weight="bold",
    )
    t.set_gid("gfv-title")
    if style.subtitle:
        ax.text(
            lay.W / 2.0, lay.H - 1.9, style.subtitle, ha="center", va="center",
            fontsize=style.title_font_size * 0.68,
            fontfamily=style.title_font_family, color=style.title_color,
        )


# --------------------------------------------------------------------------- #
# chromosome + centromere
# --------------------------------------------------------------------------- #
def _arm_points(l_a, l_b, half_h, round_a, round_b, n=12):
    """Closed polygon (in local ``(l, t)``) for a bar segment that is rounded on
    the *a* end and/or the *b* end (radius = ``half_h``, i.e. a capsule end)."""
    r = half_h
    la = l_a + (r if round_a else 0.0)
    lb = l_b - (r if round_b else 0.0)
    pts = [(la if round_a else l_a, -r), (lb if round_b else l_b, -r)]
    if round_b:
        for k in range(1, n + 1):
            ang = -math.pi / 2 + math.pi * k / n
            pts.append((lb + r * math.cos(ang), r * math.sin(ang)))
    else:
        pts.append((l_b, r))
    pts.append((la if round_a else l_a, r))
    if round_a:
        for k in range(1, n + 1):
            ang = math.pi / 2 + math.pi * k / n
            pts.append((la + r * math.cos(ang), r * math.sin(ang)))
    return pts


def _draw_chromosome(
    ax, lay, index, name, length, style, ref_bp,
    chrom_genes, chrom_cents, color_map, fam_segs, fam_tips, labels,
):
    tel = style.telomere_length
    total_bp = length + 2 * tel
    ref = ref_bp if ref_bp is not None else total_bp
    pl = lay.placer(index)
    bar_len = max(pl.along_max * (total_bp / ref), lay.body_h + 0.05)
    bh = lay.body_h
    r = bh / 2.0

    def map_pos(pos):
        return r + (pos / length) * (bar_len - 2 * r) if length else r

    # ---- centromere geometry -------------------------------------------
    cen_span = None
    draw_constriction = (
        style.show_centromeres and not chrom_cents.empty
        and style.centromere_style == "constriction"
    )
    if not chrom_cents.empty:
        cen_span = (int(chrom_cents["start"].min()), int(chrom_cents["end"].max()))

    edge_kw = dict(linewidth=style.body_edge_width,
                   edgecolor=style.body_edge_color, joinstyle="round")

    if draw_constriction and cen_span:
        p_split = min(max(map_pos(cen_span[0]), r * 1.3), bar_len - r * 1.3)
        q_split = min(max(map_pos(cen_span[1]), p_split + r * 0.6), bar_len - r * 1.1)
        # p-arm (rounded outer / flat inner)
        p = Polygon(pl.poly(_arm_points(0.0, p_split, r, True, False)),
                    closed=True, facecolor=style.body_color, zorder=2, **edge_kw)
        p.set_gid("gfv-chrom-" + _slug(name))
        ax.add_patch(p)
        # q-arm (flat inner / rounded outer)
        ax.add_patch(Polygon(
            pl.poly(_arm_points(q_split, bar_len, r, False, True)),
            closed=True, facecolor=style.body_color, zorder=2, **edge_kw,
        ))
        # primary constriction (dark bowtie)
        lm = (p_split + q_split) / 2.0
        ax.add_patch(Polygon(
            pl.poly([(p_split, r), (lm, r * 0.22), (q_split, r),
                     (q_split, -r), (lm, -r * 0.22), (p_split, -r)]),
            closed=True, facecolor=style.centromere_color,
            edgecolor=style.body_edge_color, linewidth=style.body_edge_width * 0.8,
            zorder=4,
        ))
    else:
        body = Polygon(
            pl.poly(_arm_points(0.0, bar_len, r, True, True)),
            closed=True, facecolor=style.body_color, zorder=2, **edge_kw,
        )
        body.set_gid("gfv-chrom-" + _slug(name))
        ax.add_patch(body)
        if (
            style.show_centromeres and cen_span
            and style.centromere_style == "marker"
        ):
            lm = min(max(map_pos(sum(cen_span) / 2.0), r), bar_len - r)
            ax.add_patch(Polygon(
                pl.poly([(lm, r * 1.4), (lm - r * 0.7, 0), (lm, -r * 1.4),
                         (lm + r * 0.7, 0)]),
                closed=True, facecolor=style.centromere_color,
                edgecolor="none", zorder=5,
            ))

    # ---- gene marks --------------------------------------------------
    # along-axis page-space extent of this chromosome bar -- labels for this
    # chromosome are decluttered only within this window (plus a small pad) so a
    # dense cluster on one chromosome never pushes labels onto / off another.
    e0, e1 = pl.pt(0.0, 0.0), pl.pt(bar_len, 0.0)
    if lay.orient == "horizontal":
        anchor = (name, min(e0[0], e1[0]), max(e0[0], e1[0]))
    else:
        anchor = (name, min(e0[1], e1[1]), max(e0[1], e1[1]))

    fam_here = {}
    for _, g in chrom_genes.iterrows():
        fam = g["gene_family"]
        if fam not in color_map:
            continue
        gl = map_pos((g["start"] + g["end"]) / 2.0)
        pos_side = (g["strand"] != "-") if style.split_by_strand else True
        sgn = 1.0 if pos_side else -1.0
        edge_t = sgn * r
        tip_t = sgn * (r + lay.tick_len)
        fam_segs.setdefault(fam, []).append(pl.seg(gl, edge_t, gl, tip_t))
        tip_xy = pl.pt(gl, tip_t)
        if style.tick_style in ("lollipop", "triangle", "arrow"):
            fam_tips.setdefault(fam, []).append((tip_xy[0], tip_xy[1], pos_side))
        if style.label_mode == "gene_id":
            labels.append((tip_xy[0], tip_xy[1], sgn, str(g["gene_id"]),
                           color_map[fam], anchor))
        elif style.label_mode == "family":
            fam_here.setdefault((fam, pos_side), []).append(gl)

    if style.label_mode == "family":
        for (fam, pos_side), gls in fam_here.items():
            gls.sort()
            sgn = 1.0 if pos_side else -1.0
            tip = pl.pt(gls[len(gls) // 2], sgn * (r + lay.tick_len))
            labels.append((tip[0], tip[1], sgn, fam, color_map[fam], anchor))

    # ---- chromosome name --------------------------------------------
    if lay.orient == "horizontal":
        nx, ny = pl.pt(-0.3, 0.0)
        ha, va = "right", "center"
    else:
        nx, ny = pl.pt(bar_len + 0.4, 0.0)
        ha, va = "center", "top"
    ax.text(nx, ny, str(name), ha=ha, va=va, clip_on=False,
            fontsize=style.chrom_label_font_size,
            fontfamily=style.chrom_label_font_family, color=style.chrom_label_color)


# --------------------------------------------------------------------------- #
def _emit_family_marks(ax, style, color_map, fam_segs, fam_tips):
    for fam, segs in fam_segs.items():
        lc = LineCollection(
            segs, colors=[color_map[fam]], linewidths=style.tick_width,
            alpha=style.tick_alpha, zorder=6, capstyle="round",
        )
        lc.set_gid("gfv-fam-" + _slug(fam))
        ax.add_collection(lc)
    horiz = style.orientation == "horizontal"
    for fam, tips in fam_tips.items():
        pos_pts = [(x, y) for x, y, p in tips if p]
        neg_pts = [(x, y) for x, y, p in tips if not p]
        size = max(10, (style.tick_width * 7) ** 1.05)
        if style.tick_style == "lollipop":
            m_pos = m_neg = "o"
        elif horiz:
            m_pos, m_neg = "^", "v"
        else:
            m_pos, m_neg = ">", "<"
        for pts, marker in ((pos_pts, m_pos), (neg_pts, m_neg)):
            if not pts:
                continue
            pc = ax.scatter(
                [p[0] for p in pts], [p[1] for p in pts], s=size,
                c=[color_map[fam]], marker=marker, edgecolors="none",
                zorder=7, alpha=style.tick_alpha,
            )
            pc.set_gid("gfv-fam-" + _slug(fam) + "-tips")


# --------------------------------------------------------------------------- #
def _emit_labels(ax, style, labels, lay):
    if not labels:
        return
    horiz = lay.orient == "horizontal"
    rotated = style.label_mode == "gene_id" and horiz
    char_cm = style.label_font_size * PT_TO_CM
    gap = char_cm * (1.15 if rotated else (1.4 if not horiz else 4.2))
    band_step = char_cm * (1.0 if rotated else 1.5)
    n_bands = 3 if rotated else 2

    # hard page bounds -- labels are decluttered along the length axis and must
    # never cross these, whatever the chromosome asks for
    if horiz:
        page_lo, page_hi = lay.m_left - 0.5, lay.W - lay.m_right
    else:
        page_lo, page_hi = max(0.4, lay.m_bottom - 0.4), lay.H - lay.m_top - 0.2

    # group by chromosome: declutter each chromosome's labels only within its
    # own along-axis span (+ a little), so a dense cluster on one chromosome
    # can't shove labels across a neighbour or off the canvas
    groups = {}
    for it in labels:
        groups.setdefault(it[5][0], []).append(it)

    for glabels in groups.values():
        _name, c_lo, c_hi = glabels[0][5]
        pad = 0.5 if horiz else 1.0
        lo, hi = c_lo - pad, c_hi + pad
        pos_n = sum(1 for it in glabels if it[2] > 0)
        need = (max(pos_n, len(glabels) - pos_n, 1) - 1) * gap
        if hi - lo < need:                      # borrow room for a busy chromosome
            grow = (need - (hi - lo)) / 2.0
            lo, hi = lo - grow, hi + grow
        lo, hi = max(lo, page_lo), min(hi, page_hi)
        if horiz and not rotated:               # centre-anchored family labels
            longest = max((len(str(it[3])) for it in glabels), default=0)
            inset = min(0.5 * longest * char_cm * 0.6, max(0.0, (hi - lo) / 2.0 - 1.0))
            lo, hi = lo + inset, hi - inset

        for want_sgn in (1.0, -1.0):
            items = [it for it in glabels if it[2] == want_sgn]
            if not items:
                continue
            cap = max(1, style.label_max_per_chrom)
            items.sort(key=lambda t: (t[0] if horiz else t[1]))
            if len(items) > cap:
                items = items[:: math.ceil(len(items) / cap)]
            spread = _declutter([(it[0] if horiz else it[1]) for it in items],
                                gap, lo, hi)
            for i, (it, sp) in enumerate(zip(items, spread)):
                tx, ty, sgn, text, color = it[0], it[1], it[2], it[3], it[4]
                reach = 0.32 + (i % n_bands) * band_step
                if horiz:
                    lx, ly = sp, ty + sgn * reach
                    ha, va = "center", ("bottom" if sgn > 0 else "top")
                else:
                    lx, ly = tx + sgn * reach, sp
                    ha, va = ("left" if sgn > 0 else "right"), "center"
                if style.label_leader_lines and (abs(lx - tx) > 0.05 or abs(ly - ty) > 0.05):
                    ax.add_line(Line2D([tx, lx], [ty, ly], color=color,
                                       linewidth=0.4, alpha=0.55, zorder=5))
                t = ax.text(lx, ly, text, rotation=90 if rotated else 0,
                            ha=ha, va=va, fontsize=style.label_font_size,
                            fontfamily=style.label_font_family,
                            color=style.label_color, zorder=8, clip_on=False)
                t.set_gid("gfv-label-" + _slug(text))


def _declutter(vals, min_gap, lo, hi):
    """Spread *vals* so neighbours are at least *min_gap* apart, keeping the
    whole set inside ``[lo, hi]``.  If they cannot all fit at *min_gap*, the gap
    is compressed so nothing spills past the bounds."""
    n = len(vals)
    if n == 0:
        return []
    span = max(hi - lo, 0.0)
    if n > 1 and (n - 1) * min_gap > span:
        min_gap = span / (n - 1)
    out = list(vals)
    for i in range(1, n):
        if out[i] - out[i - 1] < min_gap:
            out[i] = out[i - 1] + min_gap
    if out[-1] > hi:
        shift = out[-1] - hi
        out = [x - shift for x in out]
        for i in range(n - 2, -1, -1):
            if out[i + 1] - out[i] < min_gap:
                out[i] = out[i + 1] - min_gap
    if out[0] < lo:
        out[0] = lo
        for i in range(1, n):
            if out[i] - out[i - 1] < min_gap:
                out[i] = out[i - 1] + min_gap
    return out


# --------------------------------------------------------------------------- #
# legend
# --------------------------------------------------------------------------- #
def _legend_handles(style, color_map):
    handles = []
    for fam, rgb in color_map.items():
        if style.legend_marker == "line":
            handles.append(Line2D([0], [0], color=rgb, lw=3, label=fam))
        else:
            mk = "o" if style.legend_marker == "circle" else "s"
            handles.append(Line2D(
                [0], [0], marker=mk, color="none", markerfacecolor=rgb,
                markeredgecolor="none", markersize=9, label=fam,
            ))
    if style.show_centromeres:
        handles.append(Line2D(
            [0], [0], marker="D", color="none",
            markerfacecolor=style.centromere_color, markeredgecolor="none",
            markersize=9, label="Centromere",
        ))
    return handles


def _legend_loc(style):
    return {
        "outside right": ("center left", (1.005, 0.5)),
        "outside bottom": ("upper center", (0.5, -0.01)),
    }.get(style.legend_location, (style.legend_location, None))


def _auto_ncol(n):
    if n <= 6:
        return n
    if n <= 12:
        return int(math.ceil(n / 2))
    if n <= 24:
        return int(math.ceil(n / 3))
    return 6


def _draw_legend(ax, lay, style, color_map):
    handles = _legend_handles(style, color_map)
    ncol = style.legend_columns or _auto_ncol(len(handles))
    loc, anchor = _legend_loc(style)
    leg = ax.legend(
        handles=handles, loc=loc, bbox_to_anchor=anchor, ncol=ncol,
        frameon=style.legend_frame, fontsize=style.legend_font_size,
        title=style.legend_title or None,
        prop={"family": style.legend_font_family},
        borderaxespad=0.6, handletextpad=0.5, columnspacing=1.1,
    )
    leg.set_gid("gfv-legend")
    if leg.get_title() and leg.get_title().get_text():
        leg.get_title().set_fontsize(style.legend_font_size * 1.05)
    return leg


def _legend_page(style, color_map):
    fig, ax, lay = _new_page(style, 1)
    handles = _legend_handles(style, color_map)
    ncol = style.legend_columns or max(1, math.ceil(len(handles) / 25))
    leg = ax.legend(
        handles=handles, loc="center", ncol=ncol, frameon=style.legend_frame,
        fontsize=style.legend_font_size,
        title=style.legend_title or "Gene family",
        prop={"family": style.legend_font_family},
    )
    leg.set_gid("gfv-legend")
    return fig


# --------------------------------------------------------------------------- #
# saving
# --------------------------------------------------------------------------- #
def save_figures(figs, outdir, basename, formats, dpi=None):
    """Write *figs* in each requested format.  Returns ``{fmt: [paths]}``."""
    os.makedirs(outdir, exist_ok=True)
    formats = [f.lower().lstrip(".") for f in formats]
    out = {}
    for fmt in formats:
        ext = "jpg" if fmt == "jpeg" else fmt
        paths = []
        if fmt == "pdf" and len(figs) > 1:
            path = os.path.join(outdir, basename + ".pdf")
            with PdfPages(path) as pp:
                for fig in figs:
                    pp.savefig(fig, facecolor=fig.get_facecolor())
            paths.append(path)
        else:
            for i, fig in enumerate(figs, start=1):
                suffix = "" if len(figs) == 1 else ".p{}".format(i)
                path = os.path.join(outdir, "{}{}.{}".format(basename, suffix, ext))
                kw = {"facecolor": fig.get_facecolor()}
                if fmt in _RASTER:
                    kw["dpi"] = dpi or 200
                if fmt in ("jpg", "jpeg"):
                    kw["pil_kwargs"] = {"quality": 92}
                fig.savefig(path, **kw)
                paths.append(path)
        out[fmt] = paths
    return out


def figure_to_svg(fig):
    buf = _io.BytesIO()
    fig.savefig(buf, format="svg", facecolor=fig.get_facecolor())
    return buf.getvalue().decode("utf-8")


def close_all(figs):
    for fig in figs:
        plt.close(fig)
