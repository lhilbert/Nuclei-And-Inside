"""Spot detection inside one nucleus: a line-by-line port of the MATLAB OP block.

Order, as in `Segmentation_FeatureExtraction.m`:

    nuclear mask eroded in-plane by 0.5 um            l.455-459   erode_mask_union
    crop to the bounding box of the (eroded) nucleus  l.686-701   (pipeline)
    DoG background correction on the crop             l.725-732   filters.dog_matlab
    mean + k*SD of the DoG over the nucleus mask      l.755-759   threshold_in_nucleus
    18-connected components                           l.795
    strict volume gate, n_voxels > min_vol / voxel    l.796-800   components_gated
    DBSCAN (minPts = 1) on component centroids, merge l.1223-1243 dbscan_merge

Component and spot ORDER also follows MATLAB (bwconncomp scans column-major
over [Y X Z]), so that spot k here is meant to be spot k in a MATLAB run on the
same mask (not yet verified against MATLAB; see docs/spots3d.md).
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from skimage.measure import label as cc_label

from .filters import dog_matlab


def matlab_round(x: float) -> int:
    """MATLAB round(): half away from zero (numpy rounds half to even)."""
    return int(np.sign(x) * np.floor(abs(x) + 0.5))


def erosion_radius_px(erosion_um: float, pixel_um: float) -> int:
    """round(nuc_segErosion / pixelSize), l.457."""
    return matlab_round(erosion_um / pixel_um)


def erode_mask_union(mask: np.ndarray, radius_px: int) -> np.ndarray:
    """In-plane erosion of a binary (Z, Y, X) mask by a disk of `radius_px`.

    MATLAB erodes the binary mask of ALL nuclei together (l.458), so faces
    where two nuclei touch are not eroded; this does the same and is meant to
    be applied to `labels > 0`, then intersected with each nucleus.

    Computed as `EDT > r` per plane, which is exactly erosion by the
    Euclidean disk {x^2 + y^2 <= r^2} with the image border counting as
    foreground -- MATLAB's padding for imerode. The disk is exact; MATLAB's
    default strel('disk', r) is a polygonal (N=4) approximation of it, which
    is a documented deviation (docs/spots3d.md).
    """
    if radius_px <= 0:
        return mask.astype(bool, copy=True)
    out = np.zeros(mask.shape, dtype=bool)
    for z in range(mask.shape[0]):
        plane = mask[z].astype(bool)
        if plane.all():
            out[z] = True
        elif plane.any():
            out[z] = ndi.distance_transform_edt(plane) > radius_px
    return out


def threshold_in_nucleus(dog: np.ndarray, mask: np.ndarray, n_sd: float):
    """`DoG .* mask > mean + n_sd * std` over the mask voxels, l.755-759.

    Mean and SD accumulate in float64 and the SD is the N-1 one, as MATLAB's
    std. Voxels outside the mask enter the comparison as 0, exactly like the
    MATLAB product -- they pass only if the threshold is negative.
    """
    vals = dog[mask].astype(np.float64)
    mu = float(vals.mean())
    sd = float(vals.std(ddof=1)) if vals.size > 1 else 0.0
    thr = mu + n_sd * sd
    binary = np.where(mask, dog, 0) > thr
    return binary, mu, sd, thr


def _matlab_order(lab: np.ndarray, n: int) -> np.ndarray:
    """Labels 1..n reordered as bwconncomp would number them.

    bwconncomp numbers components by their first voxel in MATLAB's linear
    (column-major, [Y X Z]) order: y fastest, then x, then z.
    """
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    nz, ny, nx = lab.shape
    z, y, x = np.indices(lab.shape, sparse=True)
    key = (z * (nx * ny) + x * ny + y).astype(np.int64)
    first = ndi.minimum(np.broadcast_to(key, lab.shape), labels=lab,
                        index=np.arange(1, n + 1))
    return np.argsort(np.asarray(first), kind="stable") + 1


def components_gated(binary: np.ndarray, voxel_um, min_vol_um3: float):
    """18-connected components, kept if n_voxels > min_vol / voxel volume (strict).

    Returns a label image numbered 1..K in MATLAB order, and the number of
    components before the gate.
    """
    dz, dy, dx = voxel_um
    lab, n = cc_label(binary, connectivity=2, return_num=True)   # 2 hops = 18-conn in 3D
    order = _matlab_order(lab, n)
    counts = np.bincount(lab.ravel(), minlength=n + 1)
    # MATLAB: minPixels = OP_minVol./(pixelSize.^2)./zStepSize, with pixelSize the x size
    min_px = min_vol_um3 / (dx ** 2) / dz
    kept = [old for old in order if counts[old] > min_px]
    remap = np.zeros(n + 1, dtype=np.int32)
    remap[kept] = np.arange(1, len(kept) + 1, dtype=np.int32)
    return remap[lab], n


def centroids_um(lab: np.ndarray, n: int, voxel_um) -> np.ndarray:
    """Unweighted voxel centroids (z, y, x) in microns of labels 1..n."""
    if n == 0:
        return np.zeros((0, 3))
    idx = np.arange(1, n + 1)
    c = np.array(ndi.center_of_mass(np.ones(lab.shape, np.float32), lab, idx))
    return c.reshape(-1, 3) * np.asarray(voxel_um, dtype=float)


def dbscan_merge(lab: np.ndarray, voxel_um, eps_um: float):
    """MATLAB `dbscan(centroids_um, eps, 1)`, then union of each cluster's voxels.

    With minPts = 1 every point is a core point, so a DBSCAN cluster is exactly
    a connected component of the graph joining centroids at distance <= eps
    (rangesearch's criterion). Chains merge: A-B and B-C within eps put A, B
    and C into one spot even if A-C is not. Clusters are numbered by their
    first member, as MATLAB's dbscan and unique() leave them.

    Returns the merged label image and, per merged spot, how many components
    it was built from.
    """
    n = int(lab.max())
    if n == 0 or eps_um <= 0:
        return lab, np.ones(n, dtype=np.int32)
    pts = centroids_um(lab, n, voxel_um)
    pairs = np.array(sorted(cKDTree(pts).query_pairs(eps_um)), dtype=np.int64).reshape(-1, 2)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, comp = connected_components(graph, directed=False)
    # renumber clusters by their lowest member, i.e. in first-seen order
    first_seen = {}
    for c in comp:
        first_seen.setdefault(int(c), len(first_seen) + 1)
    new_id = np.array([first_seen[int(c)] for c in comp], dtype=np.int32)
    remap = np.concatenate([[0], new_id]).astype(np.int32)
    n_merged = np.bincount(new_id, minlength=len(first_seen) + 1)[1:].astype(np.int32)
    return remap[lab], n_merged


def detect_in_nucleus(crop: np.ndarray, mask: np.ndarray, voxel_um, spec,
                      small_cutoff_px: float = 0.3) -> dict:
    """Run the whole MATLAB OP detection on one nucleus crop.

    `crop` is the raw channel cut to the bounding box of `mask` (both
    (Z, Y, X)); `spec` a ChannelSpec. Returns the spot label image on the crop
    (numbered in MATLAB order) and the numbers behind it.
    """
    dog = dog_matlab(crop, spec.sigma_small_um, spec.sigma_large_um, voxel_um[2],
                     small_cutoff_px)
    binary, mu, sd, thr = threshold_in_nucleus(dog, mask, spec.n_sd)
    lab, n_raw = components_gated(binary, voxel_um, spec.min_vol_um3)
    n_gated = int(lab.max())
    lab, n_merged = dbscan_merge(lab, voxel_um, spec.dbscan_eps_um)
    return dict(labels=lab, n_merged=n_merged, dog=dog, thr_mean=mu, thr_sd=sd,
                threshold=thr, n_components_raw=int(n_raw),
                n_components_gated=n_gated, n_above_threshold=int(binary.sum()))
