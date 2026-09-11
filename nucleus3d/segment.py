"""
3D nucleus segmentation from a DNA-stain channel.

Python port of the Hilbert-lab MATLAB pipeline (`PseudoTimeCourse_CroppedImages.m`
+ `otsuLimit.m`), with two additions that the MATLAB version does not have:
a watershed split for touching nuclei, and a boundary-refinement step.

Pipeline
    1. difference-of-Gaussians bandpass     suppress noise and broad background
    2. Otsu threshold on positive DoG       (otsu_limit, ported from otsuLimit.m)
    3. watershed split of touching nuclei   h-maxima seeds on a distance map
    4. boundary refinement on raw intensity edges from the image, not the DoG
    5. filter by volume, cross-section, solidity

All length parameters are in MICRONS. They are converted to pixels using the
voxel size read from the file, so the same numbers transfer between datasets
with different sampling.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass, asdict
from scipy import ndimage as ndi
from skimage.measure import label, regionprops
from skimage.morphology import h_maxima
from skimage.segmentation import watershed


# ======================================================================
# parameters
# ======================================================================

@dataclass
class SegParams:
    """
    Segmentation parameters, all lengths in microns.

    The defaults were validated on two datasets: vt-iSIM fixed zebrafish
    embryos (JF646-Hoechst, 3.4 um slabs) and cultured cells (DAPI, 3.1 um
    slabs). See README for which knob to turn when.
    """

    # --- bandpass -----------------------------------------------------
    sigma_small_um: float = 1.0
    """Fine DoG scale. Below the noise scale, above nothing of interest."""

    sigma_large_um: float = 10.0
    """Coarse DoG scale, roughly a nuclear diameter. Scale with nucleus size."""

    isotropic_filter: bool = True
    """True: true 3D filtering accounting for dz. False: 2D slice-wise (MATLAB)."""

    # --- threshold ----------------------------------------------------
    thresh_factor: float = 1.0
    """
    Multiplier on the Otsu threshold. LOWER to catch dim nuclei.

    Otsu assumes a bimodal histogram; a field containing both bright
    mitotic figures and dim interphase nuclei is effectively trimodal, and
    the threshold lands too high. On such fields 0.4-0.6 recovers the dim
    population without inflating the bright one.
    """

    # --- splitting touching nuclei ------------------------------------
    watershed_split: bool = True
    seed_depth_um: float = 1.2
    """
    h-maxima depth for watershed seeding. THE knob for split behaviour:
    raise it to split less (fixes fragmented nuclei), lower it to split
    more (fixes merged neighbours).
    """

    seed_min_distance_um: float = 4.0
    """Minimum spacing between seeds; merges seeds on a single broad plateau."""

    ignore_z_border: bool = True
    """
    Treat the top and bottom stack faces as interior when computing the
    distance map. Correct whenever nuclei are cut off in z -- see
    `_distance_map` for why it matters.
    """

    # --- boundary refinement ------------------------------------------
    refine_boundaries: bool = True
    """Re-cut boundaries on raw intensity after DoG detection. See `_refine`."""

    refine_sigma_um: float = 0.4
    refine_factor: float = 1.0
    """Multiplier on the refinement threshold; <1 grows nuclei, >1 shrinks them."""

    # --- size and shape filters ---------------------------------------
    min_volume_um3: float = 15.0
    """
    Minimum segmented volume. NOTE this is a SLAB volume when nuclei are
    z-clipped, not a nuclear volume -- the MATLAB Nuc_min_vol = 40 assumed
    whole nuclei from a full stack and does not transfer.
    """

    min_area_um2: float = 8.0
    """
    Minimum largest cross-section over z. Independent of slab thickness,
    so this is the more portable size filter of the two.
    """

    min_solidity: float = 0.7
    """
    Minimum solidity (volume / convex-hull volume). Mitotic chromosome
    masses score ~0.84 and pass; set to ~0.3 only if you see condensed
    figures being rejected.
    """

    fill_holes: bool = True
    clear_border: bool = False
    """Drop nuclei touching the XY image edge. Off by default: edge nuclei
    are still valid for counting, just not for total-intensity work."""

    def to_dict(self):
        return asdict(self)


# ======================================================================
# thresholding -- port of otsuLimit.m
# ======================================================================

def otsu_limit(values, nbins=1000, limits=(0.0, np.inf)):
    """
    Otsu threshold on a histogram of `values`, searching only within `limits`.

    The histogram spans the full data range (matching MATLAB's
    `hist(x, 1000)` bin width); the search for the optimal split is then
    restricted to bins inside `limits`. Restricting to [0, inf) on a DoG
    image excludes the large negative background lobe, which would
    otherwise drag the threshold down.
    """
    values = np.asarray(values).ravel()
    counts, edges = np.histogram(values, bins=nbins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    keep = (centers >= limits[0]) & (centers <= limits[1])
    centers, counts = centers[keep], counts[keep].astype(float)
    if counts.sum() == 0:
        raise ValueError("no histogram counts inside the requested limits")

    p = counts / counts.sum()
    omega = np.cumsum(p)                    # cumulative class weight
    mu = np.cumsum(p * centers)             # cumulative class mean
    mu_t = mu[-1]

    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b = (mu_t * omega - mu) ** 2 / (omega * (1.0 - omega))
    sigma_b = np.nan_to_num(sigma_b, nan=-np.inf, posinf=-np.inf)

    return float(centers[int(np.argmax(sigma_b))])


# ======================================================================
# internals
# ======================================================================

def _dog(vol, voxel_um, sigma_small_um, sigma_large_um, isotropic=True):
    """Difference-of-Gaussians bandpass, sigmas given in microns."""
    vol = vol.astype(np.float32)
    dz, dy, dx = voxel_um

    if isotropic:
        s_small = (sigma_small_um / dz, sigma_small_um / dy, sigma_small_um / dx)
        s_large = (sigma_large_um / dz, sigma_large_um / dy, sigma_large_um / dx)
    else:
        # match MATLAB imgaussfilt: in-plane only, each z-slice independently
        s_small = (0.0, sigma_small_um / dy, sigma_small_um / dx)
        s_large = (0.0, sigma_large_um / dy, sigma_large_um / dx)

    return ndi.gaussian_filter(vol, s_small) - ndi.gaussian_filter(vol, s_large)


def _distance_map(mask, voxel_um, ignore_z_border=True):
    """
    Euclidean distance transform of `mask`, in microns.

    Why the z-border handling matters: in a thin slab the nuclei are cut off
    at the top and bottom of the stack. A plain EDT counts those cut faces as
    background, so no voxel is ever further than half the slab thickness from
    "outside". The distance map collapses into a z-dominated plateau whose
    local maxima have nothing to do with nuclear centres, and the watershed
    shreds single nuclei into pieces.

    With `ignore_z_border`, the first and last slices are replicated before
    the transform so the cut faces behave as interior, and the distance then
    reflects genuine lateral geometry.
    """
    if not ignore_z_border:
        return ndi.distance_transform_edt(mask, sampling=voxel_um)

    nz = mask.shape[0]
    # Pad past any distance we care about, but cap it: distances saturate
    # beyond a nuclear radius while padding cost grows linearly, and an
    # uncapped pad on a wide field makes this the memory bottleneck.
    pad = int(np.clip(round(12.0 / voxel_um[0]), 1, 64))
    padded = np.concatenate(
        [np.repeat(mask[:1], pad, axis=0), mask, np.repeat(mask[-1:], pad, axis=0)],
        axis=0,
    )
    return ndi.distance_transform_edt(padded, sampling=voxel_um)[pad:pad + nz]


def _split_touching(mask, voxel_um, seed_min_distance_um, seed_depth_um,
                    ignore_z_border=True):
    """
    Distance-transform watershed with h-maxima seeding.

    Seeds are regional maxima of the distance map that stand at least
    `seed_depth_um` above the saddle connecting them to any higher maximum.
    That criterion describes the *neck* between two nuclei rather than raw
    peak height, which makes it far less prone to fragmenting a single
    nucleus than plain peak picking.
    """
    dist = _distance_map(mask, voxel_um, ignore_z_border)

    seeds = (h_maxima(dist, seed_depth_um) > 0) & mask
    markers = label(seeds, connectivity=3)

    # collapse seeds that sit closer together than one nuclear radius
    if markers.max() > 1:
        keep, coords = [], ndi.center_of_mass(seeds, markers,
                                              range(1, markers.max() + 1))
        for c in np.atleast_2d(coords):
            cu = np.asarray(c) * np.asarray(voxel_um)
            if all(np.linalg.norm(cu - k) >= seed_min_distance_um for k in keep):
                keep.append(cu)
        if len(keep) < markers.max():
            merged = np.zeros_like(markers)
            for i, k in enumerate(keep, start=1):
                idx = tuple(int(round(v / s)) for v, s in zip(k, voxel_um))
                idx = tuple(min(max(v, 0), n - 1) for v, n in zip(idx, mask.shape))
                merged[idx] = i
            markers = merged

    if markers.max() == 0:
        return label(mask, connectivity=2)

    return watershed(-dist, markers, mask=mask)


def _refine(labels, vol, voxel_um, p):
    """
    Re-cut nucleus boundaries on the raw intensity.

    The DoG bandpass finds and separates nuclei well, but its zero-crossing
    sits *inside* the true edge: subtracting the coarse Gaussian removes some
    of the nucleus's own low-frequency content, so the thresholded region is
    systematically eroded (cross-sections come out 35-50% low).

    Here the DoG labels serve only as seeds, and the final boundary comes from
    an Otsu threshold on the lightly smoothed raw image, which follows the
    visible chromatin edge. Seeded watershed keeps neighbours apart while each
    nucleus grows out to its own intensity edge.
    """
    sm = ndi.gaussian_filter(vol.astype(np.float32),
                             tuple(p.refine_sigma_um / s for s in voxel_um))

    fg = sm > p.refine_factor * otsu_limit(sm, limits=(-np.inf, np.inf))
    fg |= labels > 0                      # never shrink below the detected cores
    fg = np.stack([ndi.binary_fill_holes(sl) for sl in fg])

    grown = watershed(-sm, labels, mask=fg)
    grown[~fg] = 0                        # drop foreground containing no seed
    return grown


def _drop_xy_border(labels):
    """Remove objects touching the XY edge (nuclei legitimately touch z)."""
    edge = np.zeros(labels.shape, dtype=bool)
    edge[:, 0, :] = edge[:, -1, :] = edge[:, :, 0] = edge[:, :, -1] = True
    for lb in np.unique(labels[edge]):
        if lb:
            labels[labels == lb] = 0
    return labels


# ======================================================================
# main entry point
# ======================================================================

def segment_nuclei(vol, voxel_um, params=None, return_intermediates=False):
    """
    Segment nuclei in a 3D DNA-stain volume.

    Parameters
    ----------
    vol : (Z, Y, X) array
        The DNA-stain channel.
    voxel_um : (dz, dy, dx)
        Voxel size in microns.
    params : SegParams or None
    return_intermediates : bool
        Also return the DoG image, threshold and raw mask (for QC figures).

    Returns
    -------
    labels : (Z, Y, X) int array, 0 = background
    props : DataFrame, one row per retained nucleus
    intermediates : dict, only if `return_intermediates`
    """
    p = params or SegParams()
    voxel_volume = float(np.prod(voxel_um))
    pixel_area = voxel_um[1] * voxel_um[2]

    # 1. bandpass
    dog = _dog(vol, voxel_um, p.sigma_small_um, p.sigma_large_um,
               p.isotropic_filter)

    # 2. threshold
    thresh = otsu_limit(dog, limits=(0.0, np.inf))
    mask = dog > p.thresh_factor * thresh
    if p.fill_holes:
        mask = np.stack([ndi.binary_fill_holes(sl) for sl in mask])

    # 3. split touching nuclei
    labels = (_split_touching(mask, voxel_um, p.seed_min_distance_um,
                              p.seed_depth_um, p.ignore_z_border)
              if p.watershed_split else label(mask, connectivity=2))

    # 4. refine boundaries against the raw image
    if p.refine_boundaries:
        labels = _refine(labels, vol, voxel_um, p)

    if p.clear_border:
        labels = _drop_xy_border(labels)

    # 5. measure and filter
    min_voxels = p.min_volume_um3 / voxel_volume
    rows, keep = [], np.zeros(int(labels.max()) + 1, dtype=bool)

    for r in regionprops(labels, intensity_image=vol.astype(np.float32)):
        if r.area < min_voxels:
            continue

        # largest cross-section over z -- a size criterion independent of
        # how much of the nucleus the slab happens to contain
        max_area_um2 = float(r.image.sum(axis=(1, 2)).max()) * pixel_area
        if max_area_um2 < p.min_area_um2:
            continue

        try:
            solidity = r.area / r.image_convex.sum()
        except Exception:            # degenerate hull (flat or tiny object)
            solidity = 0.0
        if solidity < p.min_solidity:
            continue

        keep[r.label] = True
        zc, yc, xc = r.centroid
        rows.append(dict(
            label=int(r.label),
            volume_um3=r.area * voxel_volume,
            max_area_um2=max_area_um2,
            equiv_diam_xy_um=2.0 * np.sqrt(max_area_um2 / np.pi),
            solidity=solidity,
            n_voxels=int(r.area),
            centroid_z_um=zc * voxel_um[0],
            centroid_y_um=yc * voxel_um[1],
            centroid_x_um=xc * voxel_um[2],
            touches_xy_border=bool(r.bbox[1] == 0 or r.bbox[2] == 0
                                   or r.bbox[4] == vol.shape[1]
                                   or r.bbox[5] == vol.shape[2]),
        ))

    labels = np.where(keep[labels], labels, 0)
    props = pd.DataFrame(rows)
    if not props.empty:
        props = props.sort_values("volume_um3", ascending=False).reset_index(drop=True)

    if return_intermediates:
        return labels, props, dict(dog=dog, threshold=thresh, mask=mask)
    return labels, props
