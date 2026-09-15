"""Destriping, background flattening and the noise estimate. All of it in microns.

The noise estimate is the load-bearing part. Every threshold downstream is a multiple of the
response of pure noise at the same parameters, so what `noise_sd` means decides what all of
them mean - which is why the choice of estimator is a declared parameter
(`DetectParams.noise_estimator`) and not a line inside a function.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d, zoom

from .grid import Grid


def destripe_plane(plane: np.ndarray, axis: int = 0, smooth_px: float = 25.0) -> np.ndarray:
    """Remove a fixed-pattern column offset estimated over the whole frame.

    Measured on this data: the column-median high-pass has sd 27-73 ADU against a per-pixel
    noise sd of 412-2122, i.e. 3-6% of noise. Small, but coherent along the whole 2280-row
    frame, so it is worth removing and costs almost nothing.

    **Estimated on the FULL frame**, because a column median over a nucleus-sized crop is
    dominated by real structure rather than by the stripe.
    """
    med = np.median(plane, axis=axis)
    offset = med - gaussian_filter1d(med, smooth_px, mode="nearest")
    return plane - (offset[None, :] if axis == 0 else offset[:, None])


def flatten_in_mask(vol: np.ndarray, mask: np.ndarray, sigma_vox) -> np.ndarray:
    """Background flattening by NORMALISED CONVOLUTION over the nuclear mask.

    ``background = G_s * (I . M) / (G_s * M)``, so the in-mask estimate uses only in-mask data
    and is independent of the cytoplasm by construction.

    This is not decoration. A plain Gaussian high-pass on a padded crop has infinite support,
    so it averages the cytoplasm into the background estimate within ~sigma of the boundary.
    The cytoplasm here is BRIGHTER than the nucleus (measured in/out median ratio 0.98), so a
    plain high-pass systematically over-subtracts a shell just inside the envelope - which is
    exactly the shell the antenna signal lives in.
    """
    m = mask.astype(np.float32)
    num = gaussian_filter(vol.astype(np.float32) * m, sigma_vox, mode="nearest")
    den = gaussian_filter(m, sigma_vox, mode="nearest")
    bg = np.divide(num, den, out=np.zeros_like(num), where=den > 1e-3)
    out = vol.astype(np.float32) - bg
    out[den <= 1e-3] = 0.0
    return out


#: MAD -> sd for a normal.
MAD_TO_SD = 1.4826
#: A difference of two INDEPENDENT samples has sqrt(2) x the sd of one. Verified rather than
#: assumed: the estimate is flat in lag on a white control to four decimal places, which bounds
#: any noise-correlation bias at ~3%.
DIFF_TO_SINGLE = 1.4142135623730951

NOISE_ESTIMATORS = ("first_difference_x", "in_mask_mad")


def robust_sd(a: np.ndarray, mask: np.ndarray | None = None) -> float:
    """MAD-based sd of the values themselves. **CONTAMINATED BY STRUCTURE, differentially.**

    Kept because it is what the earlier published artefacts were produced with, and because
    being able to reproduce them matters. Its defect is not an arithmetic bug - it is that a MAD
    over an unordered voxel set cannot distinguish structure from noise even in principle, and
    it is taken INSIDE the nuclear mask, over exactly the voxels that hold the signal.

    Measured against a first-difference estimate it is inflated 1.62x in structure-bearing
    nuclei against 1.06x in non-bearing ones. Since ``gamma = 5 * noise_sd`` and the frangi
    response goes as ``(S/gamma)^2`` below saturation, the response FALLS with the structure it
    is meant to detect, and falls most where there is most to detect.
    """
    v = a[mask] if mask is not None else a.ravel()
    if v.size == 0:
        return 0.0
    return float(MAD_TO_SD * np.median(np.abs(v - np.median(v))))


def noise_sd_first_difference(a: np.ndarray, mask: np.ndarray | None = None,
                              axis: int = -1, lag: int = 1) -> float:
    """Noise sd from the MAD of a lag-1 difference. Insensitive to SMOOTH structure.

    Differencing cancels anything that varies slowly compared with the lag, which is what a
    diffuse nucleoplasmic pool is. Measured: against space-filling smooth structure the in-mask
    MAD inflates 10.3x where this moves 1.09x. Against ONE filament NEITHER moves, because a
    filament covers ~1.5% of in-mask voxels and a MAD tolerates 50% - so the inflation the MAD
    suffers is the diffuse pool, not the filaments.

    It does not cancel a filament entirely. At ~2 samples per FWHM a filament is about two
    pixels wide and some of it survives the difference; at 5 samples per FWHM almost none does.
    Validating this on well-sampled data alone would be optimistic.

    **The in-mask pair restriction is not optional, and not for the obvious reason.** The
    intuition is that the mask boundary is a cliff, so dropping the restriction would inflate
    the estimate. Measured, it does the opposite and by more: `flatten_in_mask` sets everything
    outside the support to exactly 0.0, 43% of x-pairs then lie entirely outside the mask with a
    difference of exactly zero, and those zeros flood the median - the estimate collapses to
    0.24x the in-mask value. Under-estimation is the dangerous direction, because gamma is
    5 * noise_sd and too small a gamma INFLATES the response and manufactures foreground.

    **Compute this on float32, before any integer round.** Recomputing it from a stored int16
    volume differs by up to 11.2%: a MAD of integer differences is itself an integer, quantised
    in 1-ADU steps, which is 12.5% of the statistic at the dim end.
    """
    a = np.asarray(a)
    if a.ndim == 0 or a.shape[axis] <= lag:
        return 0.0
    hi = [slice(None)] * a.ndim
    lo = [slice(None)] * a.ndim
    hi[axis] = slice(lag, None)
    lo[axis] = slice(None, -lag)
    d = a[tuple(hi)] - a[tuple(lo)]
    d = d[mask[tuple(hi)] & mask[tuple(lo)]] if mask is not None else d.ravel()
    if d.size == 0:
        return 0.0
    return float(MAD_TO_SD * np.median(np.abs(d - np.median(d))) / DIFF_TO_SINGLE)


def noise_sd(a: np.ndarray, mask: np.ndarray | None = None,
             estimator: str = "first_difference_x") -> float:
    """The one entry point. Which estimator is a PARAMETER, never a code choice."""
    if estimator == "first_difference_x":
        # axis -1 is x. Coordinates are (z, y, x) everywhere in this package.
        return noise_sd_first_difference(a, mask, axis=-1)
    if estimator == "in_mask_mad":
        return robust_sd(a, mask)
    raise ValueError(f"noise_estimator must be one of {NOISE_ESTIMATORS}, got {estimator!r}")


# ------------------------------------------------------------------------------------- the crops

@dataclass
class Crop:
    """One nucleus's actin, flattened, with the mask and the grid that place it in microns."""
    nucleus_uid: str
    label_id: int
    actin_flat: np.ndarray        # float32, native sampling, background removed
    mask: np.ndarray              # bool, same shape
    grid: Grid                    # native spacing, origin = the crop's corner in field microns
    noise_sd: float
    raw_median: float


def crop_bounds(bbox_vox, coarse_shape, bin_factor, margin_um, coarse_grid, native_shape):
    """A nucleus bbox on the coarse label grid -> native-grid slices with a physical margin."""
    (z0, z1), (y0, y1), (x0, x1) = bbox_vox
    mz = int(round(margin_um / coarse_grid.dz_um))
    my = int(round(margin_um / coarse_grid.dy_um))
    Z, Y, X = coarse_shape
    zz0, zz1 = max(0, z0 - mz), min(Z, z1 + mz)
    yy0, yy1 = max(0, y0 - my), min(Y, y1 + my)
    xx0, xx1 = max(0, x0 - my), min(X, x1 + my)
    nZ, nY, nX = native_shape
    return ((zz0, zz1), (yy0, yy1), (xx0, xx1),
            slice(max(0, zz0), min(nZ, zz1)),
            slice(yy0 * bin_factor, min(nY, yy1 * bin_factor)),
            slice(xx0 * bin_factor, min(nX, xx1 * bin_factor)))


class PlaneCache:
    """Destriped native planes of ONE field's actin channel, evicted least-recently-used.

    The read is expensive and it is shared. Reading one plane decodes EVERY channel of that
    frame to return one, and `destripe_plane` then runs on the full frame. Nuclei in one field
    overlap heavily in z, so reading inside the per-nucleus loop decodes and destripes every
    shared plane once per nucleus that touches it.

    The budget is real. One field here is 169 planes of 2280x2588 float32 = 4.0 GB, which must
    not be held whole. Ask for planes in increasing z and an LRU of a few tens of planes reads
    each one about once.
    """

    def __init__(self, field, role: str, destripe: bool, axis: int, budget_bytes: int):
        self._field, self._role = field, role
        self._destripe, self._axis = destripe, axis
        self._budget = int(budget_bytes)
        self._planes: OrderedDict[int, np.ndarray] = OrderedDict()
        self._bytes = 0
        self.reads = 0
        self.hits = 0

    def get(self, z: int) -> np.ndarray:
        """The destriped plane at native index `z`. The array is SHARED - never write to it."""
        hit = self._planes.get(z)
        if hit is not None:
            self._planes.move_to_end(z)
            self.hits += 1
            return hit
        pl = self._field.plane(z, self._role).astype(np.float32)
        if self._destripe:
            pl = destripe_plane(pl, axis=self._axis)
        self.reads += 1
        self._planes[z] = pl
        self._bytes += pl.nbytes
        while self._bytes > self._budget and len(self._planes) > 1:
            _, old = self._planes.popitem(last=False)
            self._bytes -= old.nbytes
        return pl

    @property
    def summary(self) -> str:
        total = self.reads + self.hits
        if not total:
            return "no plane reads"
        return f"{self.reads} plane reads for {total} requests ({self.hits / total:.0%} reused)"


def make_crop(cache: PlaneCache, labels: np.ndarray, label_id: int, bounds, bin_factor: int,
              native: Grid, params, nucleus_uid: str) -> Crop:
    """Cut one nucleus out of a field and flatten it inside its own mask."""
    (cz, _), (cy, _), (cx, _), zsl, ysl, xsl = bounds
    ny_c = (ysl.stop - ysl.start + bin_factor - 1) // bin_factor
    nx_c = (xsl.stop - xsl.start + bin_factor - 1) // bin_factor
    m_coarse = (labels[zsl, cy:cy + ny_c, cx:cx + nx_c] == int(label_id)).astype(np.float32)

    target = (m_coarse.shape[0], ysl.stop - ysl.start, xsl.stop - xsl.start)
    mask = zoom(m_coarse, (1.0, target[1] / m_coarse.shape[1], target[2] / m_coarse.shape[2]),
                order=1) > 0.5

    vol = np.empty(target, np.float32)
    for i, zi in enumerate(range(zsl.start, zsl.stop)):
        vol[i] = cache.get(zi)[ysl, xsl]

    flat = flatten_in_mask(vol, mask, native.sigma_vox(params.detect.background_sigma_um))
    # On float32 `flat`, BEFORE any integer round. See noise_sd_first_difference.
    sd = noise_sd(flat, mask, estimator=params.detect.noise_estimator)
    origin = (zsl.start * native.dz_um, ysl.start * native.dy_um, xsl.start * native.dx_um)
    g = Grid(native.dz_um, native.dy_um, native.dx_um, *origin)
    return Crop(nucleus_uid, int(label_id), flat, mask, g, sd,
                float(np.median(vol[mask])) if mask.any() else float("nan"))
