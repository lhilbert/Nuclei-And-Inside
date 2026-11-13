"""MATLAB `imgaussfilt` on a z-stack, and the background-correcting DoG built on it.

`imgaussfilt` applied to a 3D array filters each z-plane independently with a
2D Gaussian; there is no smoothing along z (that would be `imgaussfilt3`).
Its defaults are matched exactly:

    kernel size  2*ceil(2*sigma) + 1   ->  scipy radius = ceil(2*sigma)
    padding      'replicate'           ->  scipy mode = 'nearest'

scipy's own defaults (truncate at 4 sigma, mode 'reflect') differ in both,
and the difference is largest exactly where it matters here: on a cropped
nucleus, whose bounding box edge is often within a kernel width of a spot.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import ndimage as ndi


def imgaussfilt2d(vol: np.ndarray, sigma_px: float) -> np.ndarray:
    """Per-plane 2D Gaussian of a (Z, Y, X) array, MATLAB `imgaussfilt` defaults."""
    r = int(math.ceil(2.0 * sigma_px))
    return ndi.gaussian_filter(vol, sigma=(0.0, sigma_px, sigma_px), mode="nearest",
                               radius=(0, r, r))


def dog_matlab(vol: np.ndarray, sigma_small_um: float, sigma_large_um: float,
               pixel_um: float, small_cutoff_px: float = 0.3) -> np.ndarray:
    """Background-corrected image, `Segmentation_FeatureExtraction.m` l.725-732.

        G(small) - G(large)    if sigma_small / pixel > cutoff
        I        - G(large)    otherwise (a sub-0.3 px Gaussian is the identity)

    Sigmas are converted with the IN-PLANE pixel size only, as in MATLAB.
    """
    s_small = sigma_small_um / pixel_um
    s_large = sigma_large_um / pixel_um
    low = imgaussfilt2d(vol, s_large)
    if s_small > small_cutoff_px:
        return imgaussfilt2d(vol, s_small) - low
    return vol - low
