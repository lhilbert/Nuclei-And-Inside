#!/usr/bin/env python
"""
=============================================================================
 STEP 1 of 2 -- images to tables.  Copy this file per experiment, edit the
 SETTINGS block, and run it. Nothing outside that block needs to change.
=============================================================================

    python scripts/run_segmentation.py        # hours; then step 2:
    python scripts/run_analysis.py            # seconds

Segments every field, quantifies every nucleus, and writes

    <OUTPUT_DIR>/nuclei_measurements.csv      one row per nucleus
    <OUTPUT_DIR>/field_summary.csv            one row per field
    <OUTPUT_DIR>/run_parameters.json          settings + worker plan
    <OUTPUT_DIR>/qc/*.png                     per-field validation figures
    <OUTPUT_DIR>/nucleus_boxes/*.ome.tif      one substack per nucleus (OPT-IN)
    <OUTPUT_DIR>/nucleus_boxes_index.csv      substack index, by nucleus_uid

Before the first run on a NEW dataset, check the channel names -- the folder
name is not evidence of what is in the file:

    python -c "from nucleus3d import describe_file; print(describe_file('yourfile.nd2'))"

Then look at the QC figures in <OUTPUT_DIR>/qc/ before trusting the table.
=============================================================================
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")        # headless: the library does not set this for you

# Work whether or not the package has been pip-installed: put the repository
# root on the import path, so a fresh clone runs as-is.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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

N_WORKERS = "auto"
"""
Fields to process at once: 1 for serial, an integer, or "auto".

"auto" is memory-driven, not core-driven, and the difference is large.
Segmenting one 31 x 716 x 794 field peaks near 3 GB, so what limits the
pool is how much RAM is FREE, not how many cores exist. "auto" measures the
real peak on the first field, in an isolated subprocess, and sizes the pool
from that measurement.
"""

MAX_WORKERS = None
"""Hard cap on the pool, e.g. 4 to leave the machine usable. None = no cap."""

MEMORY_FRACTION = 0.75
"""Share of currently-available RAM the pool may claim. Lower it if you
want to keep working while a run is going."""

# --- outputs ------------------------------------------------------------
SAVE_QC_FIGURES = True
"""One five-panel validation figure per field. Keep on."""

SAVE_NUCLEUS_BOXES = False
"""
Write one 3D OME-TIFF substack per nucleus -- all channels, plus the
segmentation mask as an extra channel. OPT-IN, because it is the
disk-hungry part: roughly 1-3 MB per nucleus, so ~1-2 GB for a dataset of
800.

They land in <OUTPUT_DIR>/nucleus_boxes/, named
<file>_p<position>_nuc<label>.ome.tif (e.g.
SetC_Control_004_crop_p00_nuc001.ome.tif, ~3 MB each on the example data),
and <OUTPUT_DIR>/nucleus_boxes_index.csv indexes them by `nucleus_uid` --
the same key as in nuclei_measurements.csv, so the two join directly. Each file carries its own provenance (source file,
position, label, bounding box in source coordinates) in the OME
description; read it with nucleus3d.core.export.read_box_provenance.

On a re-run, a substack whose provenance and mask still match what would
be written is left alone rather than rewritten (reuse_boxes=True).
"""

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
        n_workers=N_WORKERS,
        max_workers=MAX_WORKERS,
        memory_fraction=MEMORY_FRACTION,
    )

    table = out["measurements"]
    n_fields = len(out["field_summary"])
    print(f"\n{len(table)} nuclei from {n_fields} fields "
          f"-> {os.path.join(OUTPUT_DIR, 'nuclei_measurements.csv')}")
    if SAVE_NUCLEUS_BOXES:
        print(f"substacks -> {os.path.join(OUTPUT_DIR, 'nucleus_boxes')}/ "
              f"(indexed by nucleus_boxes.csv)")

    if len(table) and "condition" in table.columns:
        print("\nper condition:")
        print(table.groupby("condition").agg(
            fields=("field", "nunique"),
            nuclei=("label", "size"),
            median_area_um2=("max_area_um2", "median"),
            median_volume_um3=("volume_um3", "median"),
        ).round(2).to_string())

    print("\nNext: check <OUTPUT_DIR>/qc/ , then run scripts/run_analysis.py")


if __name__ == "__main__":
    main()
