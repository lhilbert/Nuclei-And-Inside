"""Synthetic fields with known filament geometry, in real physical units.

**A phantom is parameterised by the optics it is imitating and says nothing about any other.**
The two regimes here are the widefield 100x this package's defaults were measured on, at 5.0
samples per lateral FWHM, and an iSIM regime at 1.94 samples per FWHM which exists only so that
`SamplingError` can be tested against a real sampling rather than an invented one.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter

from antenna3d.grid import Grid

REGIMES = {
    # The measured PSF of the data this package's defaults come from.
    "widefield-100x": dict(
        grid=Grid(dz_um=0.2, dy_um=0.04599853515625, dx_um=0.04599853515625),
        lateral_fwhm_um=0.2302, axial_fwhm_um=0.5923),
    # The THEORETICAL PSF at NA 1.49, which is what the frozen scale block was built from.
    "widefield-100x-theory": dict(
        grid=Grid(dz_um=0.2, dy_um=0.04599853515625, dx_um=0.04599853515625),
        lateral_fwhm_um=0.180, axial_fwhm_um=0.717),
    # Present so the sampling refusal is tested against a real acquisition, not an invented one.
    "isim-100x": dict(
        grid=Grid(dz_um=0.2, dy_um=0.065074530768122, dx_um=0.065074530768122),
        lateral_fwhm_um=0.1265, axial_fwhm_um=0.4093),
}
DEFAULT_REGIME = "widefield-100x"


def grid_for(regime=DEFAULT_REGIME) -> Grid:
    return REGIMES[regime]["grid"]


def psf_sigma_vox(regime=DEFAULT_REGIME):
    r = REGIMES[regime]
    g = r["grid"]
    return (r["axial_fwhm_um"] / g.dz_um / 2.3548,
            r["lateral_fwhm_um"] / g.dy_um / 2.3548,
            r["lateral_fwhm_um"] / g.dx_um / 2.3548)


def _grid_um(shape, grid):
    zz, yy, xx = np.mgrid[0:shape[0], 0:shape[1], 0:shape[2]].astype(np.float32)
    return zz * grid.dz_um, yy * grid.dy_um, xx * grid.dx_um


def segment_distance(shape, grid, p0_um, p1_um):
    """Exact distance in MICRONS from every voxel centre to the segment p0->p1."""
    Z, Y, X = _grid_um(shape, grid)
    p0 = np.asarray(p0_um, np.float32)
    p1 = np.asarray(p1_um, np.float32)
    d = p1 - p0
    L2 = float(d @ d)
    wz, wy, wx = Z - p0[0], Y - p0[1], X - p0[2]
    t = np.clip((wz * d[0] + wy * d[1] + wx * d[2]) / max(L2, 1e-12), 0.0, 1.0)
    return np.sqrt((wz - t * d[0]) ** 2 + (wy - t * d[1]) ** 2 + (wx - t * d[2]) ** 2)


def sphere_mask(shape, grid, centre_um, radius_um):
    Z, Y, X = _grid_um(shape, grid)
    return ((Z - centre_um[0]) ** 2 + (Y - centre_um[1]) ** 2
            + (X - centre_um[2]) ** 2) <= radius_um ** 2


def render(shape, grid, segments, radius_um=0.074, amplitude=3000.0, background=8000.0,
           noise_sd=300.0, seed=0, regime=DEFAULT_REGIME, blobs=()):
    """Draw the segments, blur by the PSF, add a pedestal and noise.

    The radial profile is a GAUSSIAN of sigma `radius_um`, not a hard-edged tube. A hard tube
    thinner than the 0.2 um z voxel aliases into a string of beads, which then traces as a dozen
    fragments - an artefact of the phantom, not of the pipeline. Real filaments here are
    ~0.175 um FWHM intrinsic, i.e. sigma ~0.074 um.

    `blobs` is a list of `(centre_um, sigma_um, amplitude_multiplier)`. They are what a
    probe-absent control is actually full of, and they are BRIGHTER than the filaments - 10.0
    sigma against 5.5 - which is why a test that only ever sees clean filaments proves nothing
    about specificity.
    """
    v = np.zeros(shape, np.float32)
    for p0, p1 in segments:
        d = segment_distance(shape, grid, p0, p1)
        np.maximum(v, np.exp(-0.5 * (d / radius_um) ** 2), out=v)
    if v.max() > 0:
        v *= 1.0 / v.max()
    for centre, sigma_um, amp in blobs:
        Z, Y, X = _grid_um(shape, grid)
        r2 = ((Z - centre[0]) ** 2 + (Y - centre[1]) ** 2 + (X - centre[2]) ** 2)
        np.maximum(v, amp * np.exp(-0.5 * r2 / sigma_um ** 2), out=v)
    v = gaussian_filter(v, psf_sigma_vox(regime))
    if v.max() > 0:
        v *= amplitude / v.max()
    v += background
    if noise_sd:
        v = v + np.random.default_rng(seed).normal(0, noise_sd, shape).astype(np.float32)
    return v.astype(np.float32)


def seg_length(p0, p1) -> float:
    return float(np.linalg.norm(np.asarray(p1, float) - np.asarray(p0, float)))


def one_filament_nucleus(regime=DEFAULT_REGIME, shape=(40, 160, 160), noise_sd=300.0,
                         amplitude=3000.0, seed=0, blobs=()):
    """A spherical nucleus with one straight oblique filament of known length inside it.

    Returns `(volume, mask, grid, segments, true_length_um)`. Oblique deliberately: an
    axis-aligned filament hides exactly the errors this package exists to avoid, because
    counting voxels happens to be right for it.
    """
    grid = grid_for(regime)
    c = (shape[0] * grid.dz_um / 2, shape[1] * grid.dy_um / 2, shape[2] * grid.dx_um / 2)
    r = min(c) * 0.8
    mask = sphere_mask(shape, grid, c, r)
    p0 = (c[0] - 0.35 * r, c[1] - 0.45 * r, c[2] - 0.5 * r)
    p1 = (c[0] + 0.35 * r, c[1] + 0.45 * r, c[2] + 0.5 * r)
    segs = [(p0, p1)]
    vol = render(shape, grid, segs, amplitude=amplitude, noise_sd=noise_sd, seed=seed,
                 regime=regime, blobs=blobs)
    return vol, mask, grid, segs, seg_length(p0, p1)
