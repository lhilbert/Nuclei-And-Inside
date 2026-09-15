"""The 3D Hessian ridge filter, and the pure-noise response every threshold is stated against.

Two things here are not the library default, and both of them matter more than the choice of
filter:

1. **Anisotropic sigma.** A widefield PSF here is 0.23 um laterally and 0.59 um axially. An
   isotropic 3D Hessian on that data measures the PSF's own z-elongation and reports every
   filament as a plate. `sigma_z` matches the PSF anisotropy so the filter sees the structure
   and not the microscope.

2. **`gamma` is set in units of the measured noise sd.** `skimage`'s `frangi(gamma=None)` sets
   it per image to half the maximum Hessian norm, which makes the response incomparable between
   images. Measured on this data: with the default a PURE-NOISE IMAGE SCORES HIGHER THAN REAL
   DATA (data/noise p99.9 = 0.73); with gamma fixed in noise units the ratio is 19.7, and it is
   stable from 5 to 20 sd. Any frozen threshold on the library default is measuring the library,
   not the specimen.

**Why Frangi and not another ridge filter.** Measured on six phantom scenes at this optics:
frangi detects in 6/6, sato and meijering in 4/6 - both failing outright on the two faint ones
at 3.8 sigma peak SNR, and half of real nuclei sit below 4.8 sigma. The advantage is the
blobness term: on a synthetic tube-plus-point-source it suppresses the point relative to the
tube by 47x, against sato's 5x.

**And read `acceptance.py` before trusting that.** 47x is not enough when the false positives
are blobs four times wider than the analysis scale, because the flank of such a blob IS a ridge
at that scale - the blobness term rejects a blob measured at its own scale, not the shoulder of
a much wider one.

**Single scale, not multiscale.** Filaments here measure 1.39x the PSF, so the intrinsic width
is barely above the resolution and a scale sweep would recover the PSF rather than the object.
Apparent width is still measured per segment, and flagged as PSF-limited.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from skimage.feature import hessian_matrix_eigvals


def hessian_elems(vol: np.ndarray, sigma_vox) -> list[np.ndarray]:
    """Scale-normalised Hessian elements, upper triangle in skimage's canonical order.

    Normalisation is ``sigma_i * sigma_j`` for element (i, j) - the anisotropic generalisation
    of the usual ``sigma^2`` normalisation, which keeps the response of an ideal ridge
    independent of the scale it is measured at.
    """
    s = np.asarray(sigma_vox, dtype=float)
    v = vol.astype(np.float32, copy=False)
    out = []
    for i in range(3):
        for j in range(i, 3):
            order = [0, 0, 0]
            order[i] += 1
            order[j] += 1
            h = gaussian_filter(v, s, order=tuple(order), mode="nearest")
            out.append((h * (s[i] * s[j])).astype(np.float32))
    return out


def eigvals_by_magnitude(elems) -> np.ndarray:
    """Eigenvalues sorted by ASCENDING |lambda|: (l1, l2, l3) with |l1| <= |l2| <= |l3|."""
    ev = hessian_matrix_eigvals(elems)            # descending by value
    order = np.argsort(np.abs(ev), axis=0)
    return np.take_along_axis(ev, order, axis=0)


def frangi_3d(vol, sigma_vox, alpha=0.5, beta=0.5, gamma=None, noise_sd=None,
              gamma_in_noise_sd=5.0, black_ridges=False):
    """Frangi (1998) 3D vesselness. `gamma` is absolute, in intensity units.

    Pass `noise_sd` and let `gamma_in_noise_sd` build it. Passing neither raises, deliberately:
    the library default is not usable here and falling back to it silently is the failure this
    whole module exists to prevent.
    """
    if gamma is None:
        if noise_sd is None:
            raise ValueError("frangi_3d needs an explicit gamma or a noise_sd to build one; "
                             "the library's per-image default is not usable here - measured, "
                             "it scores pure noise ABOVE real data")
        gamma = gamma_in_noise_sd * float(noise_sd)
    l1, l2, l3 = eigvals_by_magnitude(hessian_elems(vol, sigma_vox))
    if black_ridges:
        l1, l2, l3 = -l1, -l2, -l3
    a2, b2, g2 = 2 * alpha ** 2, 2 * beta ** 2, 2 * gamma ** 2
    absl2, absl3 = np.abs(l2), np.abs(l3)
    Ra = np.divide(absl2, absl3, out=np.zeros_like(l2), where=absl3 > 0)
    denom = np.sqrt(absl2 * absl3)
    Rb = np.divide(np.abs(l1), denom, out=np.zeros_like(l1), where=denom > 0)
    S = np.sqrt(l1 ** 2 + l2 ** 2 + l3 ** 2)
    v = (1.0 - np.exp(-(Ra ** 2) / a2)) * np.exp(-(Rb ** 2) / b2) * (1.0 - np.exp(-(S ** 2) / g2))
    v[(l2 > 0) | (l3 > 0)] = 0.0          # a bright tube has l2, l3 < 0
    return np.nan_to_num(v, copy=False).astype(np.float32)


#: How the pure-noise response scales with `noise_sd`: ``null(s) = null(1) * s**d``.
#:
#: MEASURED, not assumed: 3.2e-5 maximum relative error over noise_sd 1-5000.
#:
#: The algebra says why, and the why is the point. `Ra` and `Rb` are ratios of eigenvalues and
#: are scale-free; the Hessian is linear, so `S` scales with the noise; and `gamma` is ITSELF in
#: units of noise_sd, so `S^2 / 2 gamma^2` cancels as well. Every term is invariant, so frangi's
#: null is CONSTANT in noise_sd. That is not a coincidence - it is exactly the property that
#: makes one frozen threshold mean the same thing on a bright nucleus and a faint one, written
#: in another form. A filter carrying a bare eigenvalue (sato, meijering) is degree 1 instead.
NULL_DEGREE_FRANGI = 0

_NULL_CACHE: dict = {}
_NULL_REF_SD = 1.0


def _null_uncached(noise_sd, sigma_vox, shape, percentile, seed, alpha, beta, gamma_in_noise_sd):
    rng = np.random.default_rng(seed)
    n = rng.normal(0.0, float(noise_sd), shape).astype(np.float32)
    r = frangi_3d(n, sigma_vox, alpha=alpha, beta=beta, noise_sd=noise_sd,
                  gamma_in_noise_sd=gamma_in_noise_sd)
    return float(np.percentile(r, percentile))


def noise_null(noise_sd: float, sigma_vox, shape=(24, 192, 192), percentile=99.9, seed=0,
               alpha=0.5, beta=0.5, gamma_in_noise_sd=5.0) -> float:
    """The response percentile of pure noise of the same sd at the same parameters.

    **Every threshold downstream is a multiple of THIS.** That is what makes a frozen threshold
    comparable across nuclei, fields and experiments rather than a statement about brightness.

    Computed once per parameter set and rescaled by `NULL_DEGREE_FRANGI`, because the noise
    volume is seeded and so ``n(s) = s * n(1)`` exactly. Recomputing per nucleus costs ~0.9 s
    each and returns the same number - on a few thousand nuclei that is most of an hour spent
    re-deriving a constant.
    """
    noise_sd = float(noise_sd)
    key = (tuple(round(float(s), 9) for s in np.atleast_1d(sigma_vox)), tuple(shape),
           float(percentile), int(seed), float(alpha), float(beta), float(gamma_in_noise_sd))
    if key not in _NULL_CACHE:
        _NULL_CACHE[key] = _null_uncached(_NULL_REF_SD, sigma_vox, shape, percentile, seed,
                                          alpha, beta, gamma_in_noise_sd)
    return _NULL_CACHE[key] * (noise_sd / _NULL_REF_SD) ** NULL_DEGREE_FRANGI


def enhance(crop, params) -> dict:
    """The response of one crop, and the noise null it is measured against.

    Deliberately fed the UN-DENOISED volume: the filter's own smoothing is the denoising, and a
    separate denoising step would be a second, undeclared scale.
    """
    d = params.detect
    s = params.scales()
    sigma_vox = crop.grid.sigma_vox(s["sigma_um"], s["sigma_z_um"])
    resp = frangi_3d(crop.actin_flat, sigma_vox, alpha=d.frangi_alpha, beta=d.frangi_beta,
                     noise_sd=crop.noise_sd, gamma_in_noise_sd=d.gamma_in_noise_sd)
    null = noise_null(crop.noise_sd, sigma_vox, percentile=d.null_percentile,
                      alpha=d.frangi_alpha, beta=d.frangi_beta,
                      gamma_in_noise_sd=d.gamma_in_noise_sd)
    in_mask = resp[crop.mask]
    p999 = float(np.percentile(in_mask, 99.9)) if in_mask.size else float("nan")
    return {"response": resp, "null_level": float(null), "sigma_vox": sigma_vox,
            "resp_p99_in_mask": float(np.percentile(in_mask, 99)) if in_mask.size else float("nan"),
            "resp_p999_in_mask": p999,
            # How far the brightest in-mask response sits above the noise floor. This is the one
            # number that says whether there is anything here to threshold - and on a control
            # packed with bright puncta it is LARGER than on real filaments, which is why it is
            # a diagnostic and never a detection criterion.
            "separation": p999 / null if null > 0 else float("nan")}


# ------------------------------------------------------------------------- across-ridge thinning

def principal_direction(elems, l1):
    """Unit eigenvector for `l1`, the smallest eigenvalue by magnitude.

    On a bright tube this points ALONG the filament, because the intensity barely changes along
    it. Computed from the adjugate of ``H - l1 I`` rather than by a full eigendecomposition,
    which would be minutes on an 8 M voxel volume.
    """
    a, b, c, d, e, f = elems
    m0 = np.stack([a - l1, b, c], axis=-1)
    m1 = np.stack([b, d - l1, e], axis=-1)
    m2 = np.stack([c, e, f - l1], axis=-1)
    cands = [np.cross(m0, m1), np.cross(m1, m2), np.cross(m0, m2)]
    norms = [np.linalg.norm(v, axis=-1) for v in cands]
    best = np.argmax(np.stack(norms, axis=0), axis=0)
    v = np.zeros_like(cands[0])
    for i in range(3):
        sel = best == i
        v[sel] = cands[i][sel]
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return np.divide(v, n, out=np.zeros_like(v), where=n > 1e-20)


def _orthonormal_pair(v):
    """Two unit vectors spanning the plane perpendicular to v, chosen stably."""
    a = np.zeros_like(v)
    a[..., 0] = 1.0
    flip = np.abs(v[..., 0]) > 0.9
    a[flip] = np.array([0.0, 1.0, 0.0])
    u = np.cross(v, a)
    u /= np.maximum(np.linalg.norm(u, axis=-1, keepdims=True), 1e-20)
    w = np.cross(v, u)
    w /= np.maximum(np.linalg.norm(w, axis=-1, keepdims=True), 1e-20)
    return u, w


def centerline_nms(response, vol, sigma_vox, spacing_um, step_um=0.09, min_response=0.0,
                   where=None):
    """Steger-style centreline: non-maximum suppression ACROSS the ridge, in 3D.

    It takes the maximum of the response in the plane perpendicular to the local ridge
    direction, so it is indifferent to how thick the detected object is - which matters here,
    because the axial PSF makes every filament a 4.35:1 ribbon and skeletonising a ribbon finds
    a medial axis that wanders inside it.

    **THE STEP IS PHYSICAL.** `principal_direction` returns a unit vector in the (z, y, x)
    MICRON frame, because `hessian_elems` is scale-normalised in physical units. Stepping one
    VOXEL INDEX along it would move 0.2 um in z and 0.046 um in x, comparing the response at
    incomparable distances - which suppressed every voxel and made this tracer detect nothing.

    `where` restricts the four interpolations to a candidate mask. The Hessian is still computed
    on the whole volume, because a derivative needs its neighbourhood, but the comparisons only
    ever matter inside the candidate set, which is usually 1-2% of the voxels.
    """
    sp = np.asarray(spacing_um, float)
    elems = hessian_elems(vol, sigma_vox)
    l1, l2, l3 = eigvals_by_magnitude(elems)
    cand = (response > min_response) & (l2 < 0) & (l3 < 0)
    if where is not None:
        cand &= where
    keep = np.zeros(response.shape, bool)
    idx = np.nonzero(cand)
    if idx[0].size == 0:
        return keep
    v = principal_direction(elems, l1)
    u, w = _orthonormal_pair(v)
    zz, yy, xx = (i.astype(np.float32) for i in idx)
    r0 = response[idx]
    ok = np.ones(r0.shape, bool)
    for d in (u, w):
        dd = d[idx]                                  # (N, 3) physical unit vectors
        off = (dd * step_um) / sp                    # -> index offsets, per axis
        for sgn in (1.0, -1.0):
            n = map_coordinates(response,
                                [zz + sgn * off[:, 0], yy + sgn * off[:, 1],
                                 xx + sgn * off[:, 2]],
                                order=1, mode="nearest")
            ok &= r0 >= n
    keep[idx] = ok
    return keep
