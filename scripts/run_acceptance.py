#!/usr/bin/env python
"""TEMPLATE. The probe-absent control test. Run it on an output tree run_antennas.py produced.

    python run_acceptance.py

**What this is.** If you have a control with NO PROBE in it - not a vehicle control, an actual
absent reporter - the correct length density there is zero. This measures whether the detector
agrees. It is the only test here that can distinguish a working detector from a confidently
broken one, and it is worth acquiring such a control in order to be able to run it.

It needs `keep_work="small"` or `"full"` in run_antennas.py, because it RE-TRACES from the
stored binary at each `min_branch_um` rather than filtering edges - pruning a spur changes the
topology, so a shorter minimum branch is not the same graph with more edges kept.

**On the data this package was built from, it FAILS.** 66 of 66 probe-absent nuclei reported
antennas. See the README.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))          # scripts/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from antenna3d import run_acceptance                                    # noqa: E402

# ============================================================================== SETTINGS =====

# Where the output tree and the parameters come from. `run_antennas` if you ran the antenna
# pipeline on its own, `run_combined` if you ran both. The parameters MUST be the ones that
# produced the tree: this test RE-TRACES from the stored binary, and re-tracing at other
# parameters measures nothing.
from run_antennas import OUTPUT_DIR, PARAMS                             # noqa: E402
# from run_combined import ANTENNA_DIR as OUTPUT_DIR, ANTENNA_PARAMS as PARAMS

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
