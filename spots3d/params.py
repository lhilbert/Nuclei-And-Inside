"""Parameters, all lengths in MICRONS.

Every default is the value in the Hilbert-lab MATLAB script
`Segmentation_FeatureExtraction.m` (VisitorGene_PseudoTime@10fc172); the line it
comes from is cited next to it. Change them per experiment in
`scripts/run_spots.py`, not here.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, replace


@dataclass(frozen=True)
class ChannelSpec:
    """One Oligopaint channel: where it is in the .nd2, and how to detect in it.

    `index` is the position in the .nd2 channel list and is what is read.
    `expected_nd2_name` is checked against the file, never used to look the
    channel up: names are not unique (the TripleOP files carry four channels
    all called "Blank"), and a lookup by name silently takes the first one.
    """

    name: str
    """Label used in every output column and file name, e.g. 'OP1_Cy3'."""

    index: int
    """0-based channel index in the .nd2 (Fiji C<n> = index n-1)."""

    expected_nd2_name: str | None = None
    """Channel name the .nd2 must carry at `index`; None skips the check."""

    sigma_small_um: float = 0.07       # OP_segBlurSigma_small, l.100
    sigma_large_um: float = 0.3        # OP_segBlurSigma_large, l.101
    n_sd: float = 5.0                  # OP_seg_numStdDev, l.102
    dbscan_eps_um: float = 0.2         # OP_DBSCAN_epsilon, l.104; 0 = no merging
    min_vol_um3: float = 0.05          # OP_minVol, l.110; strict '>' gate, l.797-799

    quant_blur_um: float = 0.07
    """Gaussian blur before this channel is QUANTIFIED (not detected). MATLAB
    blurs only the Oligopaint quantification channel: quantBlurSigma =
    [0,0,0.07] against quantChannels = [3,1,2] with OP on channel 2, l.65-66."""


#: The 06Oct2026_TripleOPTesting layout, verified on all 15 files:
#: ['Blank', 'Cy5', 'Blank', 'Mcherry', 'GFP', 'Blank', 'Dapi', 'Blank'].
TRIPLE_OP_CHANNELS = (
    ChannelSpec("OP1_Cy3", index=3, expected_nd2_name="Mcherry"),
    ChannelSpec("OP2_A488", index=4, expected_nd2_name="GFP"),
    ChannelSpec("OP3_Cy5", index=1, expected_nd2_name="Cy5"),
)


@dataclass(frozen=True)
class SpotParams:
    channels: tuple = field(default_factory=lambda: TRIPLE_OP_CHANNELS)
    """Channels detected in. Every spot is also quantified in all of them,
    as MATLAB quantifies every OP spot in every quantChannel (l.1380-1406)."""

    erosion_um: float = 0.5
    """nuc_segErosion, l.70. The nuclear mask is eroded in-plane by a disk of
    radius round(erosion_um / pixel) before detection (l.455-459). Only the
    mask USED here is eroded; the label volume on disk is never changed.
    0 disables it."""

    small_sigma_cutoff_px: float = 0.3
    """Below this the small Gaussian is skipped and the DoG is I - G(large),
    l.725-732."""

    qc_max_spots: int = 4
    """A nucleus with MORE spots than this is flagged (MATLAB:
    fractionAboveFour, l.1993-2003). Flagged only; nothing is excluded."""

    image_dtype: str = "float32"
    """MATLAB imageClass = 'single', l.46. Means and SDs accumulate in float64
    regardless, as in MATLAB (l.744-745)."""

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def with_(self, **kw) -> "SpotParams":
        return replace(self, **kw)
