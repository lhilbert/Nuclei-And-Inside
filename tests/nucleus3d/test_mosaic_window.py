"""
Tests for tile-window sizing. Data-free: `auto_window_um` reads one column.
"""
import numpy as np
import pandas as pd
import pytest

from nucleus3d.analysis.mosaic import auto_window_um


def _areas(diameters_um):
    """A table of nuclei with the given equivalent diameters."""
    d = np.asarray(diameters_um, dtype=float)
    return pd.DataFrame({"max_area_um2": np.pi * (d / 2.0) ** 2})


def test_window_covers_the_largest_nucleus_by_default():
    df = _areas([5, 8, 10])
    w = auto_window_um(df, margin=1.15, floor=0.0)
    assert w >= 1.15 * 10 - 1e-9          # ceil, so never under the margin
    assert w == pytest.approx(np.ceil(1.15 * 10))


def test_quantile_ignores_a_single_large_outlier():
    """
    The reason the parameter exists: one merged pair should not set the
    scale for every tile in the figure.
    """
    df = _areas([9.0] * 99 + [30.0])
    full = auto_window_um(df, margin=1.15, floor=0.0, quantile=1.0)
    trimmed = auto_window_um(df, margin=1.15, floor=0.0, quantile=0.90)
    assert full > 30.0 and trimmed < 12.0


def test_window_is_monotonic_in_quantile():
    df = _areas(np.linspace(6, 15, 50))
    ws = [auto_window_um(df, margin=1.15, floor=0.0, quantile=q)
          for q in (0.5, 0.75, 0.9, 1.0)]
    assert all(b >= a for a, b in zip(ws, ws[1:]))


def test_floor_applies_and_missing_column_falls_back():
    assert auto_window_um(_areas([2.0]), floor=12.0) == 12.0
    assert auto_window_um(pd.DataFrame({"other": [1, 2]}), floor=12.0) == 12.0
    assert auto_window_um(pd.DataFrame({"max_area_um2": [np.nan]}),
                          floor=12.0) == 12.0
