"""
spots3d -- Oligopaint spot detection inside nuclei, ported from the Hilbert-lab
MATLAB pipeline (`Segmentation_FeatureExtraction.m`, VisitorGene_PseudoTime).

A sibling of `antenna3d`: it reads nucleus3d's label volumes
(`labels/<stem>_pNN.tif`) instead of segmenting nuclei, and its rows carry the
same `nucleus_uid`, so they join to `nuclei_measurements.csv`.

    params    ChannelSpec, SpotParams -- MATLAB defaults, line-cited
    io        lazy per-channel .nd2 reads; the label volume, grid-checked
    filters   MATLAB imgaussfilt (2D per plane) and the DoG on it
    detect    erosion, mean + k*SD threshold, 18-conn components, volume gate, DBSCAN
    features  volume, solidity, centroid, mid-plane elongation, intensity ratios
    pipeline  detect_field (arrays), run_field / run_folder (files)
    qc        overview and per-nucleus figures

See docs/spots3d.md for the MATLAB mapping and every known deviation.
"""

from .params import ChannelSpec, SpotParams, TRIPLE_OP_CHANNELS
from .pipeline import detect_field, run_field, run_folder

__all__ = ["ChannelSpec", "SpotParams", "TRIPLE_OP_CHANNELS",
           "detect_field", "run_field", "run_folder"]
