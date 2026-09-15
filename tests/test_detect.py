"""Noise, gamma, and the noise null every threshold is a multiple of.

These four properties are what make one frozen threshold mean the same thing on a bright
nucleus and a faint one. If any of them breaks, nothing downstream announces it.
"""
import numpy as np
import pytest

from antenna3d.enhance import frangi_3d, noise_null
from antenna3d.preprocess import flatten_in_mask, noise_sd
from antenna3d.trace import binarize
from antenna3d.params import Params

SIGMA = (1.55, 1.7, 1.7)


def _tube_and_noise(sd=100.0, shape=(24, 96, 96), seed=0):
    rng = np.random.default_rng(seed)
    noise = rng.normal(0, sd, shape).astype(np.float32)
    zz, yy, xx = np.mgrid[0:shape[0], 0:shape[1], 0:shape[2]]
    tube = 5 * sd * np.exp(-(((yy - 48.) ** 2) + ((zz - 12.) * 4.) ** 2) / (2 * 1.7 ** 2))
    return (noise + tube).astype(np.float32), noise


# ------------------------------------------------------------------------------ the estimator

def test_both_estimators_agree_on_pure_noise():
    rng = np.random.default_rng(0)
    n = rng.normal(0, 100, (16, 64, 64)).astype(np.float32)
    m = np.ones_like(n, bool)
    assert noise_sd(n, m, "first_difference_x") == pytest.approx(100, rel=0.05)
    assert noise_sd(n, m, "in_mask_mad") == pytest.approx(100, rel=0.05)


def test_smooth_structure_inflates_the_mad_and_not_the_first_difference():
    """The measured defect: a MAD over an unordered voxel set cannot tell structure from noise.

    `gamma = 5 * noise_sd`, so an inflated denominator SUPPRESSES the response - and suppresses
    it most where there is most structure to detect.
    """
    rng = np.random.default_rng(0)
    shape = (16, 64, 64)
    n = rng.normal(0, 100, shape).astype(np.float32)
    zz, yy, xx = np.mgrid[0:16, 0:64, 0:64]
    smooth = 500 * np.sin(yy / 9.0) * np.cos(xx / 11.0)
    m = np.ones(shape, bool)
    mad = noise_sd(n + smooth, m, "in_mask_mad") / noise_sd(n, m, "in_mask_mad")
    fd = noise_sd(n + smooth, m, "first_difference_x") / noise_sd(n, m, "first_difference_x")
    assert mad > 2.0, "the MAD should absorb smooth structure"
    assert fd == pytest.approx(1.0, abs=0.1), "the first difference should not"


def test_an_unknown_estimator_raises_rather_than_defaulting():
    with pytest.raises(ValueError, match="noise_estimator"):
        noise_sd(np.zeros((4, 4, 4)), estimator="whatever")


# ------------------------------------------------------------------------------------- gamma

def test_gamma_in_noise_units_separates_data_from_noise_far_better():
    """The library's per-image default makes the response INCOMPARABLE between images.

    It sets gamma to half the maximum Hessian norm OF EACH IMAGE, so a pure-noise image gets a
    small gamma, saturates its own S term, and scores as if it were full of ridges. On the real
    data that puts noise ABOVE data (measured data/noise p99.9 = 0.73 against 19.7 with gamma in
    noise units). How far above depends on the scene, so what is asserted here is the part that
    is always true: one gamma for both images separates them and a per-image gamma does not.
    """
    from skimage.filters import frangi as sk_frangi
    data, noise = _tube_and_noise()
    ours = (np.percentile(frangi_3d(data, SIGMA, noise_sd=100.0), 99.9)
            / np.percentile(frangi_3d(noise, SIGMA, noise_sd=100.0), 99.9))
    theirs = (np.percentile(sk_frangi(data, sigmas=[1.7], black_ridges=False), 99.9)
              / np.percentile(sk_frangi(noise, sigmas=[1.7], black_ridges=False), 99.9))
    assert ours > 1.5, f"gamma in noise units must rank data well above noise (got {ours:.2f})"
    assert ours > 3 * theirs, (
        f"a per-image gamma should barely separate them: ours {ours:.2f}, library {theirs:.2f}")


def test_a_per_image_gamma_is_not_comparable_between_images():
    """The same structure at two brightnesses must not score differently. This is the property
    `high_k` depends on, and it is exactly what a per-image gamma destroys."""
    data, _ = _tube_and_noise(sd=100.0)
    dim, bright = data, data * 8.0
    ours_dim = np.percentile(frangi_3d(dim, SIGMA, noise_sd=100.0), 99.9)
    ours_bright = np.percentile(frangi_3d(bright, SIGMA, noise_sd=800.0), 99.9)
    assert ours_bright == pytest.approx(ours_dim, rel=0.02), (
        "scaling image and noise together must not move the response")


def test_frangi_refuses_to_guess_a_gamma():
    with pytest.raises(ValueError, match="not usable here"):
        frangi_3d(np.zeros((8, 16, 16), np.float32), SIGMA)


def test_the_noise_null_is_degree_zero_in_noise_sd():
    """Because Ra and Rb are ratios, the Hessian is linear, and gamma is itself in noise units.

    This is what makes `high_k = 3.0` mean the same thing on every nucleus.
    """
    assert noise_null(100.0, SIGMA) == pytest.approx(noise_null(1.0, SIGMA), rel=1e-6)
    assert noise_null(5000.0, SIGMA) == pytest.approx(noise_null(1.0, SIGMA), rel=1e-6)


def test_the_null_is_the_percentile_it_says_it_is():
    rng = np.random.default_rng(0)
    n = rng.normal(0, 100, (24, 192, 192)).astype(np.float32)
    r = frangi_3d(n, SIGMA, noise_sd=100.0)
    assert noise_null(100.0, SIGMA, percentile=99.9) == pytest.approx(
        float(np.percentile(r, 99.9)), rel=0.05)


# --------------------------------------------------------------------------------- thresholds

def test_the_threshold_is_in_null_units_not_an_image_percentile():
    """A percentile would force the same foreground fraction on every nucleus - which would make
    a probe-absent control pass the acceptance test by construction rather than by finding
    nothing."""
    p = Params()
    rng = np.random.default_rng(0)
    mask = np.ones((16, 64, 64), bool)
    null = 0.01
    empty = rng.uniform(0, null * 0.5, (16, 64, 64)).astype(np.float32)
    fg, thr = binarize(empty, null, mask, p.with_(
        detect=type(p.detect)(**{**p.detect.__dict__, "thin_by_nms": False})))
    assert fg.sum() == 0, "a response entirely below the null must yield NO foreground"
    assert thr["high"] == pytest.approx(p.detect.high_k * null)
    assert thr["low"] == pytest.approx(p.detect.low_k * null)


def test_flatten_in_mask_uses_only_in_mask_data():
    """A plain high-pass averages the surroundings into the estimate within ~sigma of the edge,
    and the surroundings here are BRIGHTER than the nucleus."""
    shape = (12, 48, 48)
    mask = np.zeros(shape, bool)
    mask[:, 12:36, 12:36] = True
    vol = np.full(shape, 100.0, np.float32)
    vol[~mask] = 10000.0                      # a very bright outside
    flat = flatten_in_mask(vol, mask, (2.0, 6.0, 6.0))
    inner = mask.copy()
    inner[:, :16, :] = inner[:, 32:, :] = inner[:, :, :16] = inner[:, :, 32:] = False
    assert np.abs(flat[inner]).max() < 1.0, "a flat interior must flatten to ~0 regardless of "\
                                            "how bright the outside is"
