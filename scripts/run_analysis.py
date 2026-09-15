#!/usr/bin/env python
"""
=============================================================================
 STEP 2 of 2 -- tables to figures.  Edit the SETTINGS block and run.
=============================================================================

    python scripts/run_analysis.py

Reads the table written by step 1 and needs nothing else -- no .nd2 files,
no segmentation -- EXCEPT for the mosaics, which re-read mid-plane crops
from the original images. Set MAKE_MOSAICS = False if the image data is not
reachable from here.

Writes into <FIGURE_DIR>:

    zclip_diagnostics.png    slab truncation -- read this one FIRST
    midplane_scatter.png     one panel per condition
    nucleus_gallery.png      a dozen captioned example nuclei
    nucleus_mosaic.png       dense image mosaic on the same plane
    pca_summary.png          all features at once: scree, loadings, scores
    pca_mosaic.png           the mosaic laid out on PC1/PC2
    pca_scores.csv           PC1..PC6 per nucleus, plus the input features
=============================================================================
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")        # headless: the library does not set this for you
import pandas as pd          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nucleus3d.analysis import (feature_pca, midplane_scatter,   # noqa: E402
                                nucleus_gallery, nucleus_mosaic,
                                pca_summary, zclip_diagnostics)

# =============================================================================
# SETTINGS -- edit below
# =============================================================================

RESULTS_DIR = "results"
"""Where step 1 wrote nuclei_measurements.csv."""

FIGURE_DIR = "results/figures"

N_Z = None
"""Planes per stack, for the slab diagnostic's first panel. None omits the
slab-face markers."""

CONDITION_ORDER = None
"""Panel order, e.g. ["Culture/Control", "Culture/Treated"]. None sorts
alphabetically -- which puts treatments before controls, so set it."""

CONDITION_LABELS = None
"""Short names for the panels, e.g. {"Hoechst_..._11Aug2021/Control":
"Control"}. None uses the last path element of each condition folder."""

MAKE_MOSAICS = True
"""The mosaics re-read crops from the .nd2 files; turn off if the image
data is not reachable from this machine."""

PSEUDOTIME_COLUMN = "mid_dna_cv_corr"
"""Column for the one-dimensional mosaic. "PC1" is the multivariate
alternative once you have looked at the PCA."""

# =============================================================================
# END OF SETTINGS
# =============================================================================


def main():
    table = pd.read_csv(os.path.join(RESULTS_DIR, "nuclei_measurements.csv"))
    os.makedirs(FIGURE_DIR, exist_ok=True)
    fig_path = lambda name: os.path.join(FIGURE_DIR, name)      # noqa: E731

    # FIRST, not last. If the clipping rate differs between conditions,
    # dropping clipped nuclei samples the conditions differently and every
    # figure below compares unlike populations.
    zc = zclip_diagnostics(table, fig_path("zclip_diagnostics.png"), n_z=N_Z,
                           condition_labels=CONDITION_LABELS)
    print("slab diagnostics ->", zc["path"])
    print("  clipped in z: " + ", ".join(f"{k} {v:.0f}%"
                                         for k, v in sorted(zc["clipped_pct"].items())))
    print("  clipped/unclipped median ratio: " + ", ".join(
        f"{k.replace('mid_', '')} {v:.2f}" for k, v in zc["clipped_over_inside"].items()))

    sc = midplane_scatter(table, fig_path("midplane_scatter.png"),
                          conditions=CONDITION_ORDER,
                          condition_labels=CONDITION_LABELS)
    # the two drop counts are sequential, not independent: the xy filter runs
    # on what the z filter left, so the three numbers sum to len(table)
    print(f"\nscatter -> {sc['path']}")
    print(f"  nuclei drawn: {sum(sc['kept'].values())} of {len(table)}"
          f" -- dropped {sum(sc['dropped_z_border'].values())} z-clipped,"
          f" then {sum(sc['dropped_xy_border'].values())} of the rest touching"
          f" the xy edge")

    # All features at once. Read res["dropped"] to see what was pruned, and
    # the specimen-split caveat in the README before reading a PC difference
    # as biology.
    res = feature_pca(table)
    pca_summary(res, fig_path("pca_summary.png"), conditions=CONDITION_ORDER,
                condition_labels=CONDITION_LABELS)
    res["scores"].to_csv(fig_path("pca_scores.csv"), index=False)
    ex = res["explained"]["variance_ratio"]
    print(f"\nPCA -> {fig_path('pca_summary.png')}")
    print(f"  {len(res['features'])} features from {res['n_samples']} nuclei, "
          f"PC1 {100 * ex.iloc[0]:.0f}%, PC2 {100 * ex.iloc[1]:.0f}%, "
          f"PC3 {100 * ex.iloc[2]:.0f}%")

    if MAKE_MOSAICS:
        g = nucleus_gallery(table, fig_path("nucleus_gallery.png"))
        m = nucleus_mosaic(table, fig_path("nucleus_mosaic.png"))
        p = nucleus_mosaic(table, fig_path("nucleus_pseudotime.png"),
                           grid=(8, 12), sort_by=PSEUDOTIME_COLUMN)
        nucleus_mosaic(res["scores"], fig_path("pca_mosaic.png"),
                       x="PC1", y="PC2")
        print(f"\nimage figures -> {FIGURE_DIR}/")
        print(f"  gallery {g['n_filled']}/{g['n_cells']} cells, "
              f"mosaic {m['n_filled']}/{m['n_cells']}, "
              f"pseudotime {p['n_filled']}/{p['n_cells']}, "
              f"tile {m['window_um']:.0f} um")


if __name__ == "__main__":
    main()
