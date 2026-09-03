"""Style configuration for the GFViewer renderer.

:class:`StyleConfig` is a flat dataclass holding every visual knob the renderer
understands.  It round-trips to a YAML or JSON file so a user can save a look
and reuse it (``--style my_style.yaml``), and it is also the payload the
interactive web editor sends back when the user tweaks the figure.
"""

import json
import os
from dataclasses import asdict, dataclass, field, fields

try:  # PyYAML is a hard dependency but degrade gracefully if absent
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

from gfviewer.errors import InputValidationError

_LEGEND_LOCATIONS = {
    "upper left", "upper center", "upper right",
    "center left", "center", "center right",
    "lower left", "lower center", "lower right",
    "right", "best", "outside right", "outside bottom", "none",
}
_TICK_STYLES = {"line", "lollipop", "triangle", "box", "arrow"}
_LABEL_MODES = {"none", "family", "gene_id"}
_ORIENTATIONS = {"horizontal", "vertical"}
_CENTROMERE_STYLES = {"constriction", "marker", "none"}


@dataclass
class StyleConfig:
    # ---- canvas ---------------------------------------------------------------
    orientation: str = "horizontal"     # horizontal | vertical
    page_width_cm: float = 26.0
    page_height_cm: float = 20.0        # used as-is only when auto_page_height is off
    auto_page_height: bool = True       # size the page to the content
    row_height_cm: float = 3.0          # space per chromosome across its thickness
    length_cm: float = 16.0             # length-axis budget (vertical orientation)
    dpi: int = 200
    background: str = "white"
    title: str = ""
    subtitle: str = ""
    show_titles: bool = True            # bake titles into exported figure /
    #                                    analytics-chart image files (the
    #                                    chromosome-figure title/subtitle and the
    #                                    built-in chart headings). Off = clean
    #                                    images for a manuscript with its own captions.

    # ---- layout -------------------------------------------------------------
    chromosomes_per_row: int = 0        # chromosomes packed side by side (0 = auto:
    #                                    1 for horizontal, all-across for vertical)
    only_chromosomes: list = field(default_factory=list)  # [] = every chromosome;
    #                                    otherwise restrict the figure to these ids
    only_families: list = field(default_factory=list)     # [] = every family;
    #                                    otherwise draw (and list in the legend)
    #                                    only these families
    rows_per_page: int = 24
    shared_scale: bool = True          # same bp/― scale across all chromosomes
    chromosome_gap: float = 0.28       # fraction of a slot left as gutter
    telomere_length: int = 10000       # bp drawn as rounded caps
    show_unplaced: bool = False         # also draw unplaced / stray contigs

    # ---- chromosome body ---------------------------------------------------
    body_width: float = 0.28           # in slot-height units
    body_color: str = "#f0f0f0"
    body_edge_color: str = "#333333"
    body_edge_width: float = 1.0
    centromere_style: str = "constriction"   # constriction | marker | none
    centromere_color: str = "#222222"        # dark primary-constriction patch
    show_centromeres: bool = False

    # ---- gene marks ------------------------------------------------------
    tick_style: str = "line"           # see _TICK_STYLES
    tick_length: float = 0.9           # fraction of half body-width it overshoots
    tick_width: float = 1.2
    tick_alpha: float = 0.95
    split_by_strand: bool = True       # + above axis, - below

    # ---- labels ----------------------------------------------------------
    label_mode: str = "none"           # none | family | gene_id
    label_font_size: float = 6.0
    label_font_family: str = "DejaVu Sans"
    label_color: str = "#222222"
    label_max_per_chrom: int = 40
    label_leader_lines: bool = True
    label_alternate_sides: bool = True

    # ---- legend --------------------------------------------------------
    legend_show: bool = True
    legend_location: str = "lower center"
    legend_columns: int = 0            # 0 = auto
    legend_font_size: float = 8.0
    legend_font_family: str = "DejaVu Sans"
    legend_title: str = "Gene family"
    legend_frame: bool = False
    legend_marker: str = "square"      # square | circle | line
    legend_per_page: bool = True       # repeat legend on every page
    legend_separate_page: bool = False  # dedicated legend page instead

    # ---- typography (chromosome names + title) -------------------------
    chrom_label_font_size: float = 9.0
    chrom_label_font_family: str = "DejaVu Sans"
    chrom_label_color: str = "#000000"
    title_font_size: float = 14.0
    title_font_family: str = "DejaVu Sans"
    title_color: str = "#000000"

    # ---- export ---------------------------------------------------------
    formats: list = field(default_factory=lambda: ["pdf"])

    # --------------------------------------------------------------------- #
    def validate(self):
        errs = []
        if self.legend_location not in _LEGEND_LOCATIONS:
            errs.append(
                "legend_location {!r} is not one of {}".format(
                    self.legend_location, sorted(_LEGEND_LOCATIONS)
                )
            )
        if self.tick_style not in _TICK_STYLES:
            errs.append("tick_style {!r} not in {}".format(self.tick_style, sorted(_TICK_STYLES)))
        if self.label_mode not in _LABEL_MODES:
            errs.append("label_mode {!r} not in {}".format(self.label_mode, sorted(_LABEL_MODES)))
        if self.orientation not in _ORIENTATIONS:
            errs.append("orientation {!r} not in {}".format(self.orientation, sorted(_ORIENTATIONS)))
        if self.centromere_style not in _CENTROMERE_STYLES:
            errs.append("centromere_style {!r} invalid".format(self.centromere_style))
        for name in ("page_width_cm", "page_height_cm", "body_width"):
            if getattr(self, name) <= 0:
                errs.append("{} must be > 0".format(name))
        if self.dpi < 30 or self.dpi > 1200:
            errs.append("dpi must be between 30 and 1200")
        if self.chromosomes_per_row < 0:
            errs.append("chromosomes_per_row must be >= 0 (0 = auto)")
        if not isinstance(self.only_chromosomes, (list, tuple)):
            errs.append("only_chromosomes must be a list of sequence names")
        if not isinstance(self.only_families, (list, tuple)):
            errs.append("only_families must be a list of family names")
        if not isinstance(self.show_titles, bool):
            errs.append("show_titles must be true or false")
        if self.rows_per_page < 1:
            errs.append("rows_per_page must be >= 1")
        if self.telomere_length < 0:
            errs.append("telomere_length must be >= 0")
        valid_fmt = {"pdf", "svg", "png", "jpg", "jpeg", "tif", "tiff", "eps", "ps"}
        bad = [f for f in self.formats if f.lower() not in valid_fmt]
        if bad:
            errs.append("unsupported export format(s): {}".format(", ".join(bad)))
        if errs:
            raise InputValidationError("Invalid style settings:", hints=errs)
        return self

    # --------------------------------------------------------------------- #
    def to_dict(self):
        return asdict(self)

    def to_json(self, **kw):
        return json.dumps(self.to_dict(), indent=2, **kw)

    def save(self, path):
        data = self.to_dict()
        if path.lower().endswith((".yaml", ".yml")):
            if yaml is None:  # pragma: no cover
                raise InputValidationError("PyYAML is not installed; save as .json instead.")
            with open(path, "w") as fh:
                yaml.safe_dump(data, fh, sort_keys=False)
        else:
            with open(path, "w") as fh:
                json.dump(data, fh, indent=2)

    # --------------------------------------------------------------------- #
    @classmethod
    def _field_names(cls):
        return {f.name for f in fields(cls)}

    @classmethod
    def from_dict(cls, data, strict=False):
        if data is None:
            return cls()
        if not isinstance(data, dict):
            raise InputValidationError("Style must be a mapping of option -> value.")
        known = cls._field_names()
        unknown = [k for k in data if k not in known]
        if unknown and strict:
            raise InputValidationError("Unknown style option(s): " + ", ".join(unknown))
        clean = {k: v for k, v in data.items() if k in known}
        cfg = cls(**clean)
        return cfg

    @classmethod
    def load(cls, path):
        if not os.path.isfile(path):
            raise InputValidationError("Style file not found: {!r}".format(path))
        with open(path, "r") as fh:
            text = fh.read()
        try:
            if path.lower().endswith((".yaml", ".yml")):
                if yaml is None:  # pragma: no cover
                    raise InputValidationError("PyYAML is not installed; use a .json style file.")
                data = yaml.safe_load(text)
            else:
                data = json.loads(text)
        except Exception as exc:  # noqa: BLE001
            raise InputValidationError("Could not parse style file {!r}: {}".format(path, exc))
        return cls.from_dict(data).validate()


DEFAULT_STYLE = StyleConfig()
