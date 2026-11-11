#!/usr/bin/env python
"""TEMPLATE. Copy this, edit the SETTINGS block, run it.

    python scripts/run_spots.py

Oligopaint spot detection inside nucleus3d's nuclei -- the MATLAB
`Segmentation_FeatureExtraction.m` OP block, ported. Run
scripts/run_segmentation.py first: this reads its labels/ folder and never
segments nuclei itself. A field without a label volume is listed in
failures.csv and skipped.

Check the channel indices on every new dataset:

    python -c "from nucleus3d import describe_file; print(describe_file('yourfile.nd2'))"

Channels are taken by INDEX (Fiji C<n> = index n-1) and the name is only
checked: several channels may share a name ("Blank").
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from spots3d import ChannelSpec, SpotParams, run_folder      # noqa: E402

# ============================================================================== SETTINGS =====

INPUT_PATH = "/path/to/nd2/folder"      # folder (searched recursively) or one .nd2
LABELS_DIR = "results/labels"            # nucleus3d output; read only
OUTPUT_DIR = "results_spots"             # never point this at the nucleus3d output folder

FILES = None             # None = every .nd2 with a label volume; or a list of file stems

# --- channels. Defaults per channel are the MATLAB OP values:
#     sigma_small_um=0.07, sigma_large_um=0.3, n_sd=5, dbscan_eps_um=0.2,
#     min_vol_um3=0.05, quant_blur_um=0.07. Override per channel if needed.
CHANNELS = (
    ChannelSpec("OP1_Cy3", index=3, expected_nd2_name="Mcherry"),
    ChannelSpec("OP2_A488", index=4, expected_nd2_name="GFP"),
    ChannelSpec("OP3_Cy5", index=1, expected_nd2_name="Cy5"),
)

EROSION_UM = 0.5         # MATLAB nuc_segErosion; 0 = use the nucleus3d mask as is
QC_MAX_SPOTS = 4         # nuclei with more spots are FLAGGED, never removed
QC = True                # figures in OUTPUT_DIR/qc. Look at them.

# ========================================================================== END SETTINGS =====

PARAMS = SpotParams(channels=CHANNELS, erosion_um=EROSION_UM, qc_max_spots=QC_MAX_SPOTS)

if __name__ == "__main__":
    run_folder(INPUT_PATH, LABELS_DIR, OUTPUT_DIR, PARAMS, files=FILES, qc=QC)
    print(f"\nNow look at {Path(OUTPUT_DIR) / 'qc'} before you use the tables.")
