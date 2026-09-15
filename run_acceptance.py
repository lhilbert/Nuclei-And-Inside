#!/usr/bin/env python
"""TEMPLATE. The probe-absent control test. Run it on an output tree run_analysis.py produced.

    python run_acceptance.py

**What this is.** If you have a control with NO PROBE in it - not a vehicle control, an actual
absent reporter - the correct length density there is zero. This measures whether the detector
agrees. It is the only test here that can distinguish a working detector from a confidently
broken one, and it is worth acquiring such a control in order to be able to run it.

It needs `keep_work="small"` or `"full"` in run_analysis.py, because it RE-TRACES from the
stored binary at each `min_branch_um` rather than filtering edges - pruning a spur changes the
topology, so a shorter minimum branch is not the same graph with more edges kept.

**On the data this package was built from, it FAILS.** 66 of 66 probe-absent nuclei reported
antennas. See the README.
"""
from antenna3d import run_acceptance
from run_analysis import OUTPUT_DIR, PARAMS

# ============================================================================== SETTINGS =====

# Values of the `condition` column, i.e. the names of your input subfolders.
POSITIVE = "HAKactin"         # the probe-PRESENT condition
NEGATIVE = "no-HAKactin"      # the probe-ABSENT control. The right answer here is ZERO.

# The curve. Reporting one value of min_branch_um is how a detector looks better than it is.
MIN_BRANCH_GRID = (0.25, 0.5, 1.0, 2.0)

# ========================================================================== END SETTINGS =====

if __name__ == "__main__":
    run_acceptance(OUTPUT_DIR, POSITIVE, NEGATIVE, PARAMS, grid_um=MIN_BRANCH_GRID)
    print(f"\nWritten: {OUTPUT_DIR}/acceptance.csv, acceptance_per_field.csv "
          f"and acceptance_per_nucleus.csv")
