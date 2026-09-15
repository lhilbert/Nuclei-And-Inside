#!/usr/bin/env python
"""TEMPLATE. Copy this, edit the SETTINGS block, run it. Nothing below SETTINGS needs changing.

    python run_analysis.py

Before your first run on a new dataset:

    python -c "from antenna3d import describe_file; print(describe_file('yourfile.nd2'))"

and set CHANNEL_DNA / CHANNEL_ACTIN to the OME names it prints, and DXY_UM / DZ_UM to the voxel
size. Both matter more than anything else here: a wrong channel produces a full, confident,
meaningless table, and a wrong voxel size puts every length and the whole filter scale at the
wrong physical size.
"""
from pathlib import Path

from antenna3d import (Channels, DetectParams, GraphParams, Optics, Params, SegParams,
                       TraceParams, E1_100X_FROZEN, run_folder)

# ============================================================================== SETTINGS =====

INPUT_DIR = "~/data/HAKactin_spheres/crop"
OUTPUT_DIR = "results/antenna3d"

# Put files in per-treatment subfolders of INPUT_DIR and the `condition` column appears,
# named after the subfolder. Files at the top level get condition = None.

# Files matching any of these are ignored. A deconvolved twin of a field is the classic one:
# counting both puts one field twice into every density, with nothing looking wrong.
EXCLUDE = ("*Deconvolved*",)

# --- channels. THE NAME IS NOT THE STAIN; run describe_file first. ---------------------------
CHANNEL_DNA = "DAPI"
CHANNEL_ACTIN = "GFP"

# --- optics. Everything scale-like derives from these. ---------------------------------------
DXY_UM = 0.04599853515625     # lateral voxel size. A file that disagrees is REFUSED, not run.
DZ_UM = 0.2                   # axial step
LATERAL_FWHM_UM = 0.2302      # MEASURED PSF. Use antenna3d.theoretical(525.0, NA) if you have
AXIAL_FWHM_UM = 0.5923        # no measurement - but measure it; it is worth an afternoon.
NA = 1.49
STEM_MARKER = "100x"          # a file whose name lacks this is refused. None to disable.

# The ten PSF-derived scales. None derives them from the PSF above, which is what NEW DATA
# should do. `E1_100X_FROZEN` reproduces the published widefield-100x numbers exactly; it is
# 0.78x-1.35x what derivation gives, because it was built from the theoretical PSF and then
# hand-rounded. It is all-or-nothing: you cannot freeze some of them.
FROZEN_SCALES = E1_100X_FROZEN

# --- nuclei ----------------------------------------------------------------------------------
# A folder of label volumes, one per field, named <stem>_p<NN>.tif or <stem>.tif. THIS IS THE
# INTENDED PATH: use the masks you already have, from nucleus3d or anywhere else, and neither
# cellpose nor torch is imported. Set to None to segment with cellpose instead.
LABELS_DIR = None

MIN_VOLUME_UM3 = 100.0        # below this it is debris
MAX_VOLUME_UM3 = 4000.0       # above this it is a merge
DROP_BORDER_NUCLEI = True     # nuclei touching an xy edge are cut laterally; drop them

# --- detection. Read "No threshold is re-tuned per condition" in the README first. ------------
HIGH_K = 3.0                  # hysteresis seed, in units of the PURE-NOISE response
LOW_K = 1.0                   # hysteresis grow
MIN_BRANCH_UM = 0.5           # shorter than this is not an antenna. Swept by run_acceptance.py.

MEASURE_WIDTH = True          # per-segment apparent FWHM. The slowest part; False to skip.

# --- what to run and what to keep -------------------------------------------------------------
FIELDS = None                 # None = all. Or a list of stems, to test on a few first.
LIMIT_NUCLEI = None           # None = all. Or an int, for a quick look.
MULTI_SERIES = False          # True only if one .nd2 really holds many stage positions
KEEP_WORK = "small"           # "small" (~30 kB/nucleus, needed by run_acceptance.py)
                              # "full"  (~9 MB/nucleus, adds the branch-point columns there)
                              # False   (keeps nothing; run_acceptance.py will not work)
QC = True                     # the figures. Leave this on.

# ========================================================================== END SETTINGS =====

PARAMS = Params(
    optics=Optics(lateral_fwhm_um=LATERAL_FWHM_UM, axial_fwhm_um=AXIAL_FWHM_UM,
                  dxy_um=DXY_UM, dz_um=DZ_UM, na=NA, stem_marker=STEM_MARKER,
                  frozen_scales=FROZEN_SCALES),
    channels=Channels(dna=CHANNEL_DNA, actin=CHANNEL_ACTIN),
    segment=SegParams(min_volume_um3=MIN_VOLUME_UM3, max_volume_um3=MAX_VOLUME_UM3,
                      drop_xy_border_touching=DROP_BORDER_NUCLEI),
    detect=DetectParams(low_k=LOW_K, high_k=HIGH_K),
    trace=TraceParams(min_branch_um=MIN_BRANCH_UM),
    graph=GraphParams(measure_width=MEASURE_WIDTH),
)

if __name__ == "__main__":
    run_folder(INPUT_DIR, OUTPUT_DIR, PARAMS, labels_dir=LABELS_DIR, fields=FIELDS,
               limit_nuclei=LIMIT_NUCLEI, keep_work=KEEP_WORK, qc=QC,
               multi_series=MULTI_SERIES, exclude=EXCLUDE)
    print(f"\nNow look at {Path(OUTPUT_DIR).expanduser() / 'qc'} before you use the table.")
    print("If you have a probe-absent control, run run_acceptance.py next.")
