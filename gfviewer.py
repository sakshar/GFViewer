#!/usr/bin/env python
"""Runnable shim so ``python gfviewer.py ...`` keeps working in a source checkout.

The real implementation lives in the :mod:`gfviewer` package
(:mod:`gfviewer.cli`).  The previous BioPython/BasicChromosome engine is kept,
unmodified, at :mod:`gfviewer.legacy` for one release.
"""

import sys

from gfviewer.cli import main

if __name__ == "__main__":
    sys.exit(main())
