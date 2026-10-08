"""
Tests for `total_persistence`.

Checked against landscapes whose answer can be worked out by hand, rather
than against a reference implementation -- the point of writing the
union-find out in full was to avoid a dependency, so pinning it to one would
defeat the exercise.
"""
import numpy as np
import pytest

from nucleus3d.core.quantify import total_persistence


def _flat(h, w, value=10.0):
    return np.full((h, w), value, dtype=np.float32)


def test_flat_image_has_no_finite_classes():
    """One island, never merges, excluded from the sum -> exactly zero."""
    img = _flat(20, 20)
    mask = np.ones_like(img, dtype=bool)
    assert total_persistence(img, mask, bg=0.0, smooth_px=0.0) == 0.0


def test_single_peak_has_no_finite_classes():
    """A lone hill is the surviving class, so it contributes nothing."""
    img = _flat(21, 21, 1.0)
    img[10, 10] = 5.0
    mask = np.ones_like(img, dtype=bool)
    assert total_persistence(img, mask, bg=0.0, smooth_px=0.0) == 0.0


def test_two_peaks_give_the_younger_peaks_prominence():
    """
    Two hills on a plateau, joined only through the background level.

    Peaks at 10 and 6 on a plateau of 1. Sweeping down, the second island is
    born at 6 and merges into the first at the plateau, 1, so its persistence
    is 6 - 1 = 5. Normalised by (brightest - bg) = 10, the total is 0.5.
    """
    img = _flat(9, 31, 1.0)
    img[4, 7] = 10.0
    img[4, 23] = 6.0
    mask = np.ones_like(img, dtype=bool)
    got = total_persistence(img, mask, bg=0.0, smooth_px=0.0)
    assert got == pytest.approx(0.5, abs=1e-6)


def test_fainter_second_peak_contributes_less():
    """Monotonic in prominence: a lower secondary peak scores lower."""
    base = _flat(9, 31, 1.0)
    base[4, 7] = 10.0
    mask = np.ones_like(base, dtype=bool)

    tall = base.copy(); tall[4, 23] = 8.0
    short = base.copy(); short[4, 23] = 3.0
    assert (total_persistence(tall, mask, bg=0.0, smooth_px=0.0)
            > total_persistence(short, mask, bg=0.0, smooth_px=0.0))


def test_invariant_to_intensity_scaling():
    """
    Normalisation is the whole reason the column is comparable between
    experiments: doubling exposure must not change the value.
    """
    img = _flat(9, 31, 1.0)
    img[4, 7] = 10.0
    img[4, 23] = 6.0
    mask = np.ones_like(img, dtype=bool)
    a = total_persistence(img, mask, bg=0.0, smooth_px=0.0)
    b = total_persistence(img * 3.0, mask, bg=0.0, smooth_px=0.0)
    assert a == pytest.approx(b, abs=1e-6)


def test_background_enters_only_through_the_normalisation():
    """A non-zero background rescales the span, hence the total."""
    img = _flat(9, 31, 1.0)
    img[4, 7] = 10.0
    img[4, 23] = 6.0
    mask = np.ones_like(img, dtype=bool)
    # span becomes 10 - 2 = 8, persistence still 5 -> 0.625
    assert total_persistence(img, mask, bg=2.0, smooth_px=0.0) == pytest.approx(
        0.625, abs=1e-6)


def test_mask_is_respected():
    """Pixels outside the mask take no part, so the second peak is invisible."""
    img = _flat(9, 31, 1.0)
    img[4, 7] = 10.0
    img[4, 23] = 6.0
    mask = np.zeros_like(img, dtype=bool)
    mask[:, :15] = True
    assert total_persistence(img, mask, bg=0.0, smooth_px=0.0) == 0.0


def test_too_small_or_degenerate_returns_nan():
    img = _flat(4, 2)
    assert np.isnan(total_persistence(img, np.ones_like(img, dtype=bool)))

    big = _flat(20, 20, 5.0)
    mask = np.ones_like(big, dtype=bool)
    # max equals bg -> span zero -> undefined rather than a divide-by-zero
    assert np.isnan(total_persistence(big, mask, bg=5.0, smooth_px=0.0))


def test_more_domains_score_higher_than_fewer():
    """
    The behaviour the metric exists for: many prominent domains score above
    few, which is the interphase-vs-condensed contrast in miniature.
    """
    mask = np.ones((41, 41), dtype=bool)

    many = _flat(41, 41, 1.0)
    for y in range(5, 40, 8):
        for x in range(5, 40, 8):
            many[y, x] = 9.0

    few = _flat(41, 41, 1.0)
    few[12, 12] = 9.0
    few[28, 28] = 9.0

    assert (total_persistence(many, mask, bg=0.0, smooth_px=0.0)
            > total_persistence(few, mask, bg=0.0, smooth_px=0.0))
