#!/usr/bin/env python
"""
=============================================================================
 ANALYSIS TEMPLATE -- copy this file per experiment, edit the block below,
 and run it. Nothing outside the SETTINGS block needs to change.
=============================================================================

    python run_analysis.py

Before the first run on a NEW dataset, check the channel names:

    python -c "from nucleus3d import describe_file; print(describe_file('yourfile.nd2'))"

Then look at the QC figures in <OUTPUT_DIR>/qc/ before trusting the table.
=============================================================================
"""

import os
import sys

# Work whether or not the package has been pip-installed: put this file's
# own folder on the import path, so a copied repository runs as-is.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nucleus3d import SegParams, run          # noqa: E402

# =============================================================================
# SETTINGS -- edit below
# =============================================================================

# --- data ---------------------------------------------------------------
INPUT_DIR = "/path/to/nd2/folder"
"""Folder of .nd2 files, searched recursively. A single .nd2 also works.
Files in subfolders get a `condition` column from the subfolder name, so
one folder per treatment gives you grouped output for free."""

OUTPUT_DIR = "results"

DNA_CHANNEL = "DAPI"
"""
Channel to segment -- name or index. CHECK THIS on any new dataset.

Do not trust the folder name. One dataset here is called "...JF646Hoechst..."
but has no Hoechst channel at all: the DNA stain is on 'Cy5'. Another
dataset has no DNA channel whatsoever, and the lab MATLAB script segments
nuclei from the Pol II channel instead. Both give confident, meaningless
results if the wrong channel is set.
"""

POSITIONS = None
"""Stage positions to process: None for all, or e.g. [0, 5, 10] to test
parameters quickly before committing to a full run."""

# --- outputs ------------------------------------------------------------
SAVE_QC_FIGURES = True
"""One five-panel validation figure per field. Keep on."""

SAVE_NUCLEUS_BOXES = False
"""Write one 3D OME-TIFF per nucleus (all channels + mask). Off by default:
it is the slow, disk-hungry part -- roughly 1-3 MB per nucleus."""

BOX_PAD_UM = 1.0             # lateral margin around each nucleus
BOX_INCLUDE_MASK = True      # append the segmentation mask as a channel

# --- segmentation parameters -------------------------------------------
# Defaults are validated on vt-iSIM embryo (JF646-Hoechst) and cultured-cell
# (DAPI) slab data. Each line notes the direction to move it; see README.
PARAMS = SegParams(
    # bandpass
    sigma_small_um=1.0,        # fine scale; raise if the data is noisy
    sigma_large_um=10.0,       # ~ one nuclear diameter; scale with nucleus size

    # threshold
    thresh_factor=1.0,         # LOWER (0.4-0.6) to catch dim nuclei

    # splitting touching nuclei
    seed_depth_um=1.2,         # RAISE to split less, LOWER to split more
    seed_min_distance_um=4.0,
    ignore_z_border=True,      # keep True for thin slabs / z-clipped nuclei

    # boundary refinement
    refine_boundaries=True,    # keep True: DoG alone cuts inside the edge
    refine_factor=1.0,         # <1 grows boundaries, >1 shrinks them

    # size and shape filters
    min_volume_um3=15.0,       # SLAB volume, not nuclear volume, when clipped
    min_area_um2=8.0,          # largest cross-section; portable across slabs
    min_solidity=0.7,          # mitotic figures score ~0.84 and pass
    clear_border=False,        # True drops nuclei touching the image edge
)

# =============================================================================
# END OF SETTINGS
# =============================================================================


def main():
    out = run(
        input_dir=INPUT_DIR,
        outdir=OUTPUT_DIR,
        dna_channel=DNA_CHANNEL,
        params=PARAMS,
        positions=POSITIONS,
        save_qc=SAVE_QC_FIGURES,
        save_boxes=SAVE_NUCLEUS_BOXES,
        box_pad_um=BOX_PAD_UM,
        box_include_mask=BOX_INCLUDE_MASK,
    )

    table = out["measurements"]
    if len(table) and "condition" in table.columns:
        print("\nper condition:")
        print(table.groupby("condition").agg(
            fields=("field", "nunique"),
            nuclei=("label", "size"),
            median_area_um2=("max_area_um2", "median"),
            median_volume_um3=("volume_um3", "median"),
        ).round(2).to_string())


if __name__ == "__main__":
    main()
