"""GFViewer: visualize the localization of multigene families across chromosomes.

The package is organized as:

* :mod:`gfviewer.genome`     -- read chromosome lengths (FASTA / .fai / chrom.sizes / table)
* :mod:`gfviewer.io`         -- unified feature loader (table / BED / GFF3 / GTF) + validation
* :mod:`gfviewer.palette`    -- programmatic colour-palette generation
* :mod:`gfviewer.style`      -- :class:`StyleConfig` and YAML/JSON style files
* :mod:`gfviewer.render`     -- matplotlib ideogram renderer (PDF/SVG/PNG/JPG/...)
* :mod:`gfviewer.analytics`  -- per-family statistics, clustering, localization bias
* :mod:`gfviewer.cli`        -- command-line interface
"""

__version__ = "2.0.0"

from gfviewer.errors import GFViewerError, InputValidationError  # noqa: E402

__all__ = ["__version__", "GFViewerError", "InputValidationError"]
