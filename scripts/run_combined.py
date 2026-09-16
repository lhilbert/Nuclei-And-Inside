#!/usr/bin/env python
"""
=============================================================================
 ANALYSIS TEMPLATE -- nuclei, and the antennas inside them, in one run.
 Copy this file per experiment, edit the SETTINGS block, and run it.
=============================================================================

    python scripts/run_combined.py

It runs nucleus3d over your .nd2 files, hands its label volumes to antenna3d,
and joins the two tables on `nucleus_uid`.

Before the first run on a NEW dataset, check the channel names and the voxel
size. They are the two most expensive things to get wrong, and both are
silent:

    python -c "from nucleus3d import describe_file; print(describe_file('yourfile.nd2'))"
    python -c "from antenna3d import describe_file; print(describe_file('yourfile.nd2'))"

(Yes, both -- they report different things. antenna3d also prints the voxel
size, the slab thickness and the optical regime it would demand.)

THE TWO PIPELINES DO NOT ACCEPT THE SAME DATA. antenna3d declares an optical
regime and refuses a file whose voxel size disagrees, because everything
scale-like in it derives from the PSF. The shipped defaults are widefield
100x. On iSIM sampling it raises SamplingError rather than running a filter
whose kernel is narrower than the sample spacing -- so nucleus3d will happily
process data that antenna3d correctly declines. Run scripts/run_segmentation.py alone
for those.

AND READ THE CAVEAT. On the data antenna3d was built from, a
reporter-NEGATIVE control reported antennas in 66 of 66 nuclei. Run
scripts/run_acceptance.py on your own probe-absent control before quoting any
antenna number. See docs/antenna3d.md.
=============================================================================
"""

import os
import sys

# The repository root: this template lives in scripts/, the packages and
# combine.py do not.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Qualified imports, deliberately. BOTH packages export a `SegParams` and a
# `describe_file`, with different fields and different return shapes. A
# `from ... import *` here silently binds one name to the other package's
# object.
import antenna3d as a3                                        # noqa: E402
import nucleus3d as n3                                        # noqa: E402
from combine import join_nuclei_and_antennas                  # noqa: E402

# =============================================================================
# SETTINGS -- edit below
# =============================================================================

# --- data ---------------------------------------------------------------
INPUT_DIR = "/path/to/nd2/folder"
"""Folder of .nd2 files, searched recursively. One subfolder per treatment
gives both tables a `condition` column named after the subfolder -- and
run_acceptance.py keys its arms off exactly that column, so both pipelines
must be pointed at the SAME root for it to work."""

OUTPUT_DIR = "results"
"""Both pipelines write a `field_summary.csv`, a `run_parameters.json` and a
`failures.csv`. They get one subtree each, below, so they do not overwrite
each other."""

NUCLEUS_DIR = os.path.join(OUTPUT_DIR, "nuclei")
ANTENNA_DIR = os.path.join(OUTPUT_DIR, "antennas")
JOINED_CSV = os.path.join(OUTPUT_DIR, "nuclei_and_antennas.csv")

# --- channels. THE NAME IS NOT THE STAIN. Run describe_file first. -------
DNA_CHANNEL = "DAPI"
ACTIN_CHANNEL = "GFP"

# --- optics. Everything scale-like in antenna3d derives from these. ------
DXY_UM = 0.04599853515625     # a file that disagrees is REFUSED, not run
DZ_UM = 0.2
LATERAL_FWHM_UM = 0.2302      # MEASURED PSF. a3.theoretical(525.0, NA) if you
AXIAL_FWHM_UM = 0.5923        # have no measurement -- but measure it.
NA = 1.49
STEM_MARKER = "100x"          # a file whose name lacks this is refused; None to disable
FROZEN_SCALES = a3.E1_100X_FROZEN   # None derives all ten from the PSF above

# --- what to run --------------------------------------------------------
POSITIONS = None              # None = all; [0, 5, 10] to test first
LIMIT_NUCLEI = None           # None = all; an int for a quick look
MULTI_SERIES = False          # True only if one .nd2 really holds many positions
EXCLUDE = ("*Deconvolved*",)  # a deconvolved twin counts the same field twice

STOP_AFTER_SEGMENTATION = False
"""Stop once the nuclei are segmented, so you can open the QC figures before
committing to the long half of the run. Recommended on a new dataset."""

# --- nucleus segmentation (nucleus3d) -----------------------------------
NUCLEUS_PARAMS = n3.SegParams(
    sigma_small_um=1.0,
    sigma_large_um=10.0,
    thresh_factor=1.0,         # LOWER (0.4-0.6) to catch dim nuclei
    seed_depth_um=1.2,         # RAISE to split less, LOWER to split more
    seed_min_distance_um=4.0,
    ignore_z_border=True,
    refine_boundaries=True,
    refine_factor=1.0,
    min_volume_um3=15.0,       # SLAB volume, not nuclear volume, when clipped
    min_area_um2=8.0,
    min_solidity=0.7,
    clear_border=False,
)

# --- antenna detection (antenna3d) --------------------------------------
# Read "no threshold is ever re-tuned per condition" in docs/antenna3d.md.
HIGH_K = 3.0                  # hysteresis seed, in units of the PURE-NOISE response
LOW_K = 1.0
MIN_BRANCH_UM = 0.5           # shorter than this is not an antenna
MEASURE_WIDTH = True          # per-segment apparent FWHM; the slowest part
KEEP_WORK = "small"           # "small" is what run_acceptance.py re-traces from

# antenna3d's own nucleus gate, applied to the masks nucleus3d hands over.
# It is STRICTER than nucleus3d's on purpose -- a nucleus with no usable
# denominator has no usable density -- so expect fewer antenna rows than
# nucleus rows. The join reports the difference.
ANTENNA_MIN_VOLUME_UM3 = 100.0
ANTENNA_MAX_VOLUME_UM3 = 4000.0
DROP_BORDER_NUCLEI = True

# =============================================================================
# END OF SETTINGS
# =============================================================================

ANTENNA_PARAMS = a3.Params(
    optics=a3.Optics(lateral_fwhm_um=LATERAL_FWHM_UM, axial_fwhm_um=AXIAL_FWHM_UM,
                     dxy_um=DXY_UM, dz_um=DZ_UM, na=NA, stem_marker=STEM_MARKER,
                     frozen_scales=FROZEN_SCALES),
    channels=a3.Channels(dna=DNA_CHANNEL, actin=ACTIN_CHANNEL),
    segment=a3.SegParams(min_volume_um3=ANTENNA_MIN_VOLUME_UM3,
                         max_volume_um3=ANTENNA_MAX_VOLUME_UM3,
                         drop_xy_border_touching=DROP_BORDER_NUCLEI),
    detect=a3.DetectParams(low_k=LOW_K, high_k=HIGH_K),
    trace=a3.TraceParams(min_branch_um=MIN_BRANCH_UM),
    graph=a3.GraphParams(measure_width=MEASURE_WIDTH),
)


def main():
    print("=" * 74)
    print("STAGE 1/2  nuclei")
    print("=" * 74)
    nuclei = n3.run(
        input_dir=INPUT_DIR,
        outdir=NUCLEUS_DIR,
        dna_channel=DNA_CHANNEL,
        params=NUCLEUS_PARAMS,
        positions=POSITIONS,
        save_qc=True,
        save_labels=True,          # the handoff. Without it stage 2 has no masks.
    )
    labels_dir = nuclei["labels_dir"]

    print(f"\nLook at {os.path.join(NUCLEUS_DIR, 'qc')} before you trust stage 2.")
    print("Panel 3 is the missed-nucleus check; a whole cyan nucleus is a real miss.")
    if STOP_AFTER_SEGMENTATION:
        print("\nSTOP_AFTER_SEGMENTATION is on. Set it to False and re-run "
              "to continue -- stage 1 runs again, so on a large dataset check "
              "the figures from a short POSITIONS list first.")
        return

    print("\n" + "=" * 74)
    print("STAGE 2/2  antennas, inside those nuclei")
    print("=" * 74)
    a3.run_folder(
        INPUT_DIR, ANTENNA_DIR, ANTENNA_PARAMS,
        labels_dir=labels_dir,
        positions=POSITIONS,
        limit_nuclei=LIMIT_NUCLEI,
        keep_work=KEEP_WORK,
        qc=True,
        multi_series=MULTI_SERIES,
        exclude=EXCLUDE,
    )

    print("\n" + "=" * 74)
    print("JOIN")
    print("=" * 74)
    join_nuclei_and_antennas(NUCLEUS_DIR, ANTENNA_DIR, out_csv=JOINED_CSV)

    print("\nNow, in order:")
    print(f"  1. {os.path.join(NUCLEUS_DIR, 'qc')}       segmentation")
    print(f"  2. {os.path.join(ANTENNA_DIR, 'qc', 'traces')}  is it a filament, "
          f"or the rim of a punctum?")
    print("  3. scripts/run_acceptance.py, if you have a probe-absent control.")
    print("     Without one, nothing here can tell you the antennas are real.")


if __name__ == "__main__":
    main()
