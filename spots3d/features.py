"""Per-spot features, `Segmentation_FeatureExtraction.m` l.1278-1406.

    volume       n_voxels * pixel^2 * dz                          l.1282-1283
    solidity     volume / convex volume                           l.1284
    centroid     unweighted voxel mean                            l.1285-1286
    elongation   major / minor axis of the CENTRAL z-plane of the
                 spot's own bounding box, 2D ellipse moments      l.1307-1359
    intensity    median(spot) / median(nucleus), per quantified
                 channel, after that channel's optional blur      l.1380-1406

No feature is used to filter: MATLAB measures solidity and elongation but
gates Oligopaint spots on volume only.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import regionprops

from .filters import imgaussfilt2d


def ellipse_axes_matlab(mask2d: np.ndarray):
    """MajorAxisLength, MinorAxisLength exactly as MATLAB regionprops computes them.

    Second central moments of the pixel coordinates plus 1/12 (the moment of
    a unit pixel), then the axes of the ellipse with the same moments. skimage's
    axis lengths omit the 1/12 and so differ for small objects.

    MATLAB is called on uint8(mask), i.e. a label image in which every
    foreground pixel is region 1 -- disconnected pieces count as one region,
    and so they do here.
    """
    y, x = np.nonzero(mask2d)
    n = x.size
    xc = x - x.mean()
    yc = -(y - y.mean())
    uxx = np.sum(xc ** 2) / n + 1.0 / 12.0
    uyy = np.sum(yc ** 2) / n + 1.0 / 12.0
    uxy = np.sum(xc * yc) / n
    common = math.sqrt((uxx - uyy) ** 2 + 4.0 * uxy ** 2)
    major = 2.0 * math.sqrt(2.0) * math.sqrt(uxx + uyy + common)
    minor = 2.0 * math.sqrt(2.0) * math.sqrt(max(uxx + uyy - common, 0.0))
    return major, minor


def elongation_matlab(obj_image: np.ndarray) -> float:
    """Major/minor axis ratio on the plane ceil(nz/2) (1-based) of the object's bbox."""
    nz = obj_image.shape[0]
    plane = obj_image[int(math.ceil(nz / 2.0)) - 1]
    if not plane.any():
        return float("nan")
    major, minor = ellipse_axes_matlab(plane)
    return major / minor if minor > 0 else float("nan")


def spot_features(lab: np.ndarray, voxel_um, origin_zyx) -> list[dict]:
    """One dict per spot (labels 1..K of `lab`, in that order).

    `origin_zyx` is the voxel index of the crop's corner in the field, so the
    centroids come out in field coordinates. The MATLAB-convention centroid
    (crop-local, 1-based, [x, y, z], in microns -- what MATLAB stores in
    OP_centCell) is kept too, for comparing against a MATLAB run.
    """
    dz, dy, dx = voxel_um
    oz, oy, ox = origin_zyx
    rows = []
    for r in regionprops(lab):
        cz, cy, cx = r.centroid
        try:
            solidity = float(r.solidity)
        except Exception:              # degenerate hull
            solidity = float("nan")
        rows.append(dict(
            spot_id=int(r.label),
            n_voxels=int(r.area),
            volume_um3=float(r.area) * dx * dx * dz,
            solidity=solidity,
            elongation=elongation_matlab(r.image),
            centroid_z_px=oz + cz, centroid_y_px=oy + cy, centroid_x_px=ox + cx,
            centroid_z_um=(oz + cz) * dz, centroid_y_um=(oy + cy) * dy,
            centroid_x_um=(ox + cx) * dx,
            matlab_cent_x_um=(cx + 1) * dx, matlab_cent_y_um=(cy + 1) * dx,
            matlab_cent_z_um=(cz + 1) * dz,
            bbox_z=f"{oz + r.bbox[0]}-{oz + r.bbox[3] - 1}",
        ))
    return rows


def spot_intensities(lab: np.ndarray, quant_crop: np.ndarray, mask: np.ndarray,
                     blur_um: float, pixel_um: float):
    """Median(spot) / median(nucleus) of one quantified channel, l.1380-1401.

    The blur is applied to the crop, as in MATLAB. Also returns the nucleus
    median used as the denominator.
    """
    img = quant_crop
    if blur_um > 0:
        img = imgaussfilt2d(img, blur_um / pixel_um)
    nuc_median = float(np.median(img[mask]))
    n = int(lab.max())
    if n == 0:
        return np.zeros(0), nuc_median
    spot_median = np.asarray(ndi.median(img, lab, np.arange(1, n + 1)), dtype=float)
    return spot_median / nuc_median, nuc_median
