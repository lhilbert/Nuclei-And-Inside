"""
nucleus3d -- 3D nucleus segmentation and per-nucleus quantification for
Nikon .nd2 z-stacks, plus the analysis of the resulting tables.

The package is in two halves that meet only at the measurement table:

    nucleus3d.core       images -> nuclei_measurements.csv  (+ optional
                         per-nucleus substacks)
    nucleus3d.analysis   nuclei_measurements.csv -> figures

Run them as two steps, via the templates in `scripts/`:

    python scripts/run_segmentation.py     # step 1, hours
    python scripts/run_analysis.py         # step 2, seconds

Only the handful of names below are re-exported here -- the entry point of
each half, its parameter object, and the figure functions. Everything else
is reachable through the subpackage it lives in
(`from nucleus3d.core.quantify import midplane_metrics`), which keeps the
top-level namespace readable and makes the import say which half a function
belongs to.
"""

from .core import SegParams, run
from .core.io import load_field, describe_file
from .analysis import (feature_pca, midplane_scatter, zclip_diagnostics,
                       nucleus_gallery, nucleus_mosaic, pca_summary)

__all__ = [
    "SegParams", "run", "load_field", "describe_file",
    "feature_pca", "midplane_scatter", "zclip_diagnostics",
    "nucleus_gallery", "nucleus_mosaic", "pca_summary",
]
