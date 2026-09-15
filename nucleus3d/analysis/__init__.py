"""
Step 2: tables to figures.

Everything here takes `nuclei_measurements.csv` (or the in-memory
DataFrame) and needs nothing else -- except the mosaics, which re-read
mid-plane crops from the original .nd2 files and therefore need the image
data reachable.

    features    all-feature PCA: feature_pca, feature_matrix
    plots       point figures: midplane_scatter, zclip_diagnostics,
                pca_summary
    mosaic      image figures: nucleus_gallery, nucleus_mosaic
    style       shared palette and font sizes

Two filters are applied by default by every figure function here, and both
matter more than they sound: nuclei whose widest plane is a slab face
(`mid_at_z_border`) and nuclei touching an image edge in xy
(`touches_xy_border`). See the README section on truncation.
"""

from .features import feature_pca, feature_matrix
from .plots import midplane_scatter, zclip_diagnostics, pca_summary
from .mosaic import nucleus_gallery, nucleus_mosaic, auto_window_um

__all__ = [
    "feature_pca", "feature_matrix",
    "midplane_scatter", "zclip_diagnostics", "pca_summary",
    "nucleus_gallery", "nucleus_mosaic", "auto_window_um",
]
