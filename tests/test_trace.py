"""End-to-end on a phantom whose geometry is known exactly.

The one test in this suite that can fail because the SCIENCE is wrong rather than because the
plumbing is. If the recovered length of a straight filament drifts, something upstream of the
graph moved.
"""
import numpy as np
import pytest

from antenna3d.enhance import enhance
from antenna3d.grid import Grid
from antenna3d.params import Params
from antenna3d.preprocess import Crop, flatten_in_mask, noise_sd
from antenna3d.trace import (binarize, drop_short_components, path_length_um, reject_puncta,
                             resample_to, trace)

import phantom


def crop_from_phantom(vol, mask, grid, params, uid="phantom_p00_nuc001"):
    flat = flatten_in_mask(vol, mask, grid.sigma_vox(params.detect.background_sigma_um))
    return Crop(uid, 1, flat, mask, grid, noise_sd(flat, mask, params.detect.noise_estimator),
                float(np.median(vol[mask])))


def run_detect(vol, mask, grid, params):
    crop = crop_from_phantom(vol, mask, grid, params)
    enh = enhance(crop, params)
    fg, thr = binarize(enh["response"], enh["null_level"], crop.mask, params,
                       vol=crop.actin_flat, sigma_vox=enh["sigma_vox"],
                       spacing=crop.grid.spacing)
    return crop, enh, fg, trace(fg, crop.grid, params)


# --------------------------------------------------------------------------- the work grid

def test_the_work_grid_does_not_lose_a_one_voxel_centreline():
    """scipy's zoom interpolates without prefiltering, and a 2-3 voxel structure sampled at half
    rate can fall below any threshold EVERYWHERE and vanish outright. A whole phantom filament
    once did exactly that, which is why the resample is an exact block scatter."""
    native = Grid(0.2, 0.046, 0.046)
    fg = np.zeros((30, 200, 200), bool)
    t = np.linspace(0, 1, 3000)
    pts = np.round(np.array([5, 40, 40]) + t[:, None] * np.array([20, 120, 140])).astype(int)
    fg[pts[:, 0], pts[:, 1], pts[:, 2]] = True
    out, wg = resample_to(fg, native, (0.2962, 0.1151, 0.1151))
    assert out.sum() > 0, "the filament must survive the resample"
    assert wg.dz_um > wg.dy_um, "the work grid is anisotropic in microns, by design"


def test_the_resample_carries_the_origin():
    src = Grid(0.2, 0.046, 0.046, oz_um=4.0, oy_um=2.0, ox_um=1.0)
    _, wg = resample_to(np.zeros((8, 16, 16), bool), src, (0.3, 0.1, 0.1))
    assert wg.origin == src.origin


# ------------------------------------------------------------------------------ the filament

def test_a_known_filament_is_recovered_at_the_right_length():
    p = Params()
    vol, mask, grid, _, true_um = phantom.one_filament_nucleus(noise_sd=300.0, amplitude=3000.0)
    crop, enh, fg, tr = run_detect(vol, mask, grid, p)

    assert enh["separation"] > p.detect.high_k, (
        f"a {true_um:.1f} um filament at this contrast must clear the seed "
        f"(separation {enh['separation']:.1f} vs high_k {p.detect.high_k})")
    assert fg.sum() > 0
    got = tr["length_um"]
    assert got == pytest.approx(true_um, rel=0.10), (
        f"recovered {got:.3f} um for a {true_um:.3f} um filament")
    # The raw voxel-step sum is the QC number and is ~1.3x too long. Pinned so that the two
    # cannot silently swap: `length_um` is the measurement, `skeleton_length_um` is not.
    assert tr["skeleton_length_um"] > tr["length_um"] * 1.2


def test_the_filament_is_one_component_not_several():
    """A broken filament becomes two antennas with two false tips and half the length each -
    the failure mode that corrupts a graph and is invisible to a length statistic."""
    from scipy.ndimage import label as ndlabel
    p = Params()
    vol, mask, grid, _, _ = phantom.one_filament_nucleus(noise_sd=300.0, amplitude=3000.0)
    _, _, _, tr = run_detect(vol, mask, grid, p)
    _, n = ndlabel(tr["skeleton"], np.ones((3, 3, 3)))
    assert n <= 2, f"a single filament traced as {n} components"


def test_pure_noise_traces_almost_nothing():
    """The threshold is a multiple of the noise response, so noise must not clear it."""
    p = Params()
    grid = phantom.grid_for()
    shape = (40, 160, 160)
    mask = phantom.sphere_mask(shape, grid, (shape[0] * grid.dz_um / 2, shape[1] * grid.dy_um / 2,
                                             shape[2] * grid.dx_um / 2), 2.5)
    rng = np.random.default_rng(1)
    vol = (8000 + rng.normal(0, 300, shape)).astype(np.float32)
    _, _, _, tr = run_detect(vol, mask, grid, p)
    assert tr["length_um"] < 1.0, f"pure noise traced {tr['length_um']:.2f} um"


def test_a_lower_min_branch_is_not_the_same_graph_with_more_edges():
    """It re-prunes, so it changes topology. This is why the acceptance test RE-TRACES."""
    p = Params()
    vol, mask, grid, _, _ = phantom.one_filament_nucleus(noise_sd=400.0, amplitude=3000.0)
    crop, enh, fg, _ = run_detect(vol, mask, grid, p)
    short = trace(fg, crop.grid, p, min_branch_um=0.25)
    long_ = trace(fg, crop.grid, p, min_branch_um=2.0)
    assert short["length_um"] >= long_["length_um"]
    assert short["min_branch_um"] == 0.25 and long_["min_branch_um"] == 2.0


# -------------------------------------------------------------------------------- the puncta

def test_a_closed_ring_is_rejected_and_a_filament_is_not():
    g = Grid(0.3585, 0.092, 0.092)
    sk = np.zeros((10, 40, 40), bool)
    th = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    sk[5, (20 + 3 * np.sin(th)).round().astype(int),
       (20 + 3 * np.cos(th)).round().astype(int)] = True       # a small closed rim
    out, n = reject_puncta(sk, g.spacing, 2.0, 0.35, smooth_passes=3)
    assert n == 1 and out.sum() == 0

    line = np.zeros((10, 40, 40), bool)
    line[5, 20, 5:35] = True
    out2, n2 = reject_puncta(line, g.spacing, 2.0, 0.35, smooth_passes=3)
    assert n2 == 0 and out2.sum() == line.sum(), "a filament has free ends and is not a loop"


def test_the_loop_perimeter_threshold_is_a_real_perimeter():
    """`max_loop_perimeter_um` is named for a length in microns and must be given one.

    A digitised circle's raw voxel-step sum is ~1.27x its true perimeter, so comparing the raw
    sum against 2.0 um would silently make the threshold 1.5 um. This rim's true perimeter is
    1.734 um: inside the threshold, and it must be rejected.
    """
    g = Grid(0.3585, 0.092, 0.092)
    sk = np.zeros((10, 40, 40), bool)
    th = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    sk[5, (20 + 3 * np.sin(th)).round().astype(int),
       (20 + 3 * np.cos(th)).round().astype(int)] = True
    from antenna3d.trace import component_lengths
    true_perimeter = 2 * np.pi * 3 * 0.092
    raw = component_lengths(sk, g.spacing)[1][1]
    smoothed = component_lengths(sk, g.spacing, 3)[1][1]
    assert raw / true_perimeter > 1.2, "the raw step sum is expected to overestimate"
    assert smoothed == pytest.approx(true_perimeter, rel=0.10)
    assert reject_puncta(sk, g.spacing, 2.0, 0.35, 3)[1] == 1
    assert reject_puncta(sk, g.spacing, 2.0, 0.35, 0)[1] == 0, (
        "pinned so the difference the smoothing makes stays visible")


def test_a_broken_rim_survives_the_puncta_filter():
    """The measured limitation, pinned so it is not mistaken for a working filter.

    On a real probe-absent control this rejected 22 rims out of well over a hundred, because a
    rim that is broken or fused to its neighbour HAS free ends and so is not a closed loop.
    """
    g = Grid(0.3585, 0.092, 0.092)
    sk = np.zeros((10, 40, 40), bool)
    th = np.linspace(0, 1.7 * np.pi, 400)                      # an ARC, not a closed ring
    sk[5, (20 + 3 * np.sin(th)).round().astype(int),
       (20 + 3 * np.cos(th)).round().astype(int)] = True
    _, n = reject_puncta(sk, g.spacing, 2.0, 0.35)
    assert n == 0, "the filter is topological; a broken rim is not rejected and must not be"


def test_short_components_are_dropped_by_true_length():
    g = Grid(0.3585, 0.092, 0.092)
    sk = np.zeros((10, 40, 40), bool)
    sk[5, 10, 5:30] = True          # 24 * 0.092 = 2.2 um
    sk[5, 30, 5:8] = True           # 2 * 0.092 = 0.18 um
    out, n = drop_short_components(sk, g.spacing, 0.5)
    assert n == 1
    assert path_length_um(out, g.spacing) == pytest.approx(24 * 0.092, rel=0.02)
