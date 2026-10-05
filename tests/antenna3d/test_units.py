"""Invariant 1: every length is in microns, and `Grid` is the only converter.

The failure this guards against does not crash and does not look wrong - a number in voxels is
a plausible number.
"""
import numpy as np
import pytest

from antenna3d.grid import Grid, isotropic, sample_um
from antenna3d.trace import path_length_um, smooth_polyline

import phantom


def test_grid_converts_per_axis():
    g = Grid(0.2, 0.046, 0.046)
    sz, sy, sx = g.sigma_vox(0.092, 0.3585)
    assert sz == pytest.approx(0.3585 / 0.2)
    assert sy == sx == pytest.approx(0.092 / 0.046)
    # One micron is a different number of voxels on each axis. That is the whole point.
    assert g.um_to_vox(1.0, 0) != g.um_to_vox(1.0, 1)


def test_grid_round_trip_carries_the_origin():
    g = Grid(0.2, 0.046, 0.046, oz_um=3.0, oy_um=1.5, ox_um=0.5)
    assert Grid.from_dict(g.to_dict()) == g
    assert g.index_to_um((10, 20, 30)) == pytest.approx((3.0 + 2.0, 1.5 + 0.92, 0.5 + 1.38))


def test_voxel_volume_and_isotropy():
    assert isotropic(0.092).is_isotropic
    assert not Grid(0.3585, 0.092, 0.092).is_isotropic
    assert Grid(0.2, 0.046, 0.046).voxel_volume_um3 == pytest.approx(0.2 * 0.046 ** 2)


def test_oblique_length_is_not_a_voxel_count():
    """A digitised diagonal steps by up to sqrt(3) x spacing, so counting voxels is wrong.

    This is invariant 1's specific instance and the reason `path_length_um` sums real steps.
    """
    g = Grid(0.3585, 0.092, 0.092)
    p0, p1 = np.array([2.0, 5.0, 5.0]), np.array([17.0, 45.0, 55.0])
    t = np.linspace(0, 1, 4000)
    pts = np.round(p0[None, :] + t[:, None] * (p1 - p0)[None, :]).astype(int)
    sk = np.zeros((22, 60, 70), bool)
    sk[pts[:, 0], pts[:, 1], pts[:, 2]] = True

    true_um = float(np.linalg.norm((p1 - p0) * np.array(g.spacing)))
    stepped = path_length_um(sk, g.spacing)
    naive = sk.sum() * g.dx_um            # voxel count x ONE spacing: the bug

    assert stepped > 0
    assert abs(naive / true_um - 1.0) > 0.15, "the naive count should be visibly wrong here"
    # The quantisation zigzag is real and is what the polyline smoothing removes.
    poly = np.stack(np.nonzero(sk), 1).astype(float) * np.array(g.spacing)
    poly = poly[np.argsort(poly[:, 2])]
    raw = np.linalg.norm(np.diff(poly, axis=0), axis=1).sum()
    sm = np.linalg.norm(np.diff(smooth_polyline(poly, 3), axis=0), axis=1).sum()
    assert sm < raw, "smoothing must shorten a quantised line"
    assert sm / true_um == pytest.approx(1.0, abs=0.15)


def test_smoothing_preserves_a_straight_line_and_an_arc():
    straight = np.stack([np.zeros(40), np.zeros(40), np.linspace(0, 4, 40)], 1)
    assert np.allclose(smooth_polyline(straight, 5), straight, atol=1e-9)

    th = np.linspace(0, 1.2, 200)
    arc = np.stack([np.zeros_like(th), 3 * np.sin(th), 3 * np.cos(th)], 1)
    before = np.linalg.norm(np.diff(arc, axis=0), axis=1).sum()
    after = np.linalg.norm(np.diff(smooth_polyline(arc, 3), axis=0), axis=1).sum()
    assert after / before == pytest.approx(1.0, abs=0.01), "real curvature must survive"


def test_sample_um_respects_the_origin():
    g = Grid(1.0, 1.0, 1.0, oz_um=10.0, oy_um=20.0, ox_um=30.0)
    vol = np.arange(8 * 8 * 8, dtype=float).reshape(8, 8, 8)
    got = sample_um(vol, g, np.array([[10.0 + 2, 20.0 + 3, 30.0 + 4]]))
    assert got[0] == pytest.approx(vol[2, 3, 4])


def test_phantom_filament_length_is_what_it_says():
    _, _, grid, segs, L = phantom.one_filament_nucleus()
    (p0, p1), = segs
    assert L == pytest.approx(phantom.seg_length(p0, p1))
    assert L > 3.0, "the test filament must be long enough to survive pruning"
