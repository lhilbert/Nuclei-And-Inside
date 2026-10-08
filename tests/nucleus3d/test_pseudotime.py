"""
Tests for the angular pseudo-time coordinate.

All synthetic and data-free: points are placed on a known circle so the
expected angles can be written down, which is also the only way to test the
wrap-around behaviour that makes angles different from ordinary numbers.
"""
import numpy as np
import pandas as pd
import pytest

from nucleus3d.analysis.pseudotime import (angular_pseudotime, circular_mean,
                                           circular_R, angular_distance,
                                           phase_order_check)


def _ring(n=24, cx=0.0, cy=0.0, r=2.0):
    """n points evenly spaced on a circle, in raw (unstandardised) columns."""
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return pd.DataFrame({"mx": cx + r * np.cos(a), "my": cy + r * np.sin(a)})


# ---------------------------------------------------------------- circular


def test_circular_mean_handles_the_wrap():
    """
    The reason this module exists. Angles either side of zero have mean zero,
    but an ordinary mean of the same numbers lands near pi -- halfway round
    the circle from every one of the inputs.
    """
    a = np.array([0.1, 6.2, 0.05, 6.25])        # all near zero, straddling it
    got = circular_mean(a)
    assert min(got, 2 * np.pi - got) < 0.1
    assert abs(np.mean(a) - np.pi) < 0.2        # the trap, demonstrated


def test_circular_R_spans_zero_to_one():
    assert circular_R([1.0, 1.0, 1.0]) == pytest.approx(1.0)
    uniform = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    assert circular_R(uniform) < 1e-9


def test_circular_helpers_ignore_nan_and_empty():
    assert circular_mean([np.nan, 0.5, np.nan]) == pytest.approx(0.5)
    assert np.isnan(circular_mean([]))
    assert np.isnan(circular_R([np.nan]))


def test_angular_distance_is_the_short_way_round():
    assert angular_distance([0.1], 6.23)[0] == pytest.approx(0.1532, abs=1e-3)
    assert angular_distance([0.0], np.pi)[0] == pytest.approx(np.pi)
    assert (angular_distance(np.linspace(0, 2 * np.pi, 50), 1.0) <= np.pi + 1e-9).all()


# ------------------------------------------------------------ pseudotime


def test_pi_points_at_the_requested_direction():
    """The defining property: the reference point must land on theta = pi."""
    df = _ring()
    df.loc[len(df)] = {"mx": 3.0, "my": 0.0}          # due east of the centre
    out, _ = angular_pseudotime(df, centre=(0, 0), pi_toward=(3.0, 0.0),
                                x="mx", y="my", standardize=False)
    assert out["theta"].iloc[-1] == pytest.approx(np.pi, abs=1e-9)


def test_radius_is_distance_from_the_centre():
    out, _ = angular_pseudotime(_ring(r=2.0), centre=(0, 0), pi_toward=(1, 0),
                                x="mx", y="my", standardize=False)
    assert out["radius"].to_numpy() == pytest.approx(np.full(len(out), 2.0))


def test_theta_is_in_range_and_covers_the_circle():
    out, _ = angular_pseudotime(_ring(n=60), centre=(0, 0), pi_toward=(1, 0),
                                x="mx", y="my", standardize=False)
    t = out["theta"].to_numpy()
    assert t.min() >= 0.0 and t.max() < 2 * np.pi
    assert circular_R(t) < 0.05                   # evenly spread, as placed


def test_direction_reverses_the_order():
    """cw must traverse the ring the opposite way from ccw."""
    df = _ring(n=12)
    ccw, _ = angular_pseudotime(df, centre=(0, 0), pi_toward=(1, 0), x="mx",
                                y="my", direction="ccw", standardize=False)
    cw, _ = angular_pseudotime(df, centre=(0, 0), pi_toward=(1, 0), x="mx",
                               y="my", direction="cw", standardize=False)
    # theta_ccw + theta_cw == 0 (mod 2pi) everywhere. Compared with
    # `angular_distance`, not by subtracting floats: the sums straddle the
    # wrap, so some come back as ~1e-16 and others as ~6.2831 while all of
    # them mean the same angle. Writing this test naively reproduced the
    # exact bug the module docstring warns about.
    s = np.mod(ccw["theta"].to_numpy() + cw["theta"].to_numpy(), 2 * np.pi)
    assert (angular_distance(s, 0.0) < 1e-9).all()


def test_standardisation_removes_the_scale_imbalance():
    """
    Without standardising, an axis carrying larger numbers dominates the
    angle. Stretching y by 50 must not move theta once standardised.
    """
    df = _ring(n=16)
    a, _ = angular_pseudotime(df, centre=(0, 0), pi_toward=(1, 0), x="mx",
                              y="my", standardize=True)
    stretched = df.assign(my=df.my * 50.0)
    b, _ = angular_pseudotime(stretched, centre=(0, 0), pi_toward=(1, 0),
                              x="mx", y="my", standardize=True)
    assert a["theta"].to_numpy() == pytest.approx(b["theta"].to_numpy(), abs=1e-9)


def test_params_let_a_standardised_centre_be_reported_in_raw_units():
    df = _ring(n=20, cx=5.0, cy=100.0, r=3.0)
    _, p = angular_pseudotime(df, centre=(1.0, 0.0), pi_toward=(2.0, 0.0),
                              x="mx", y="my", standardize=True)
    assert p["centre_raw"][0] == pytest.approx(1.0 * p["sd"][0] + p["mean"][0])
    assert p["centre_raw"][1] == pytest.approx(0.0 * p["sd"][1] + p["mean"][1])


def test_bad_direction_and_missing_column_are_refused():
    df = _ring()
    with pytest.raises(ValueError):
        angular_pseudotime(df, (0, 0), (1, 0), x="mx", y="my", direction="up")
    with pytest.raises(KeyError):
        angular_pseudotime(df, (0, 0), (1, 0), x="mx", y="nope")


def test_phase_order_check_detects_order_and_counts():
    df = pd.DataFrame({"theta": [0.1, 0.2, 1.0, 1.1, 2.0],
                       "radius": [1.0] * 5,
                       "phase": ["a", "a", "b", "b", "c"]})
    good = phase_order_check(df, "phase", ["a", "b", "c"])
    assert good["monotonic"] and good["n_total"] == 5
    assert [r["n"] for r in good["per_phase"]] == [2, 2, 1]
    assert not phase_order_check(df, "phase", ["c", "b", "a"])["monotonic"]
