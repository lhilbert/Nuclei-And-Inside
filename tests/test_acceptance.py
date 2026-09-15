"""The probe-absent control test, and the ways it must refuse to answer."""
import numpy as np
import pandas as pd
import pytest

from antenna3d.acceptance import mannwhitney, run_acceptance, summarise


def curve(pos_density, neg_density, mb=0.5):
    rows = [{"nucleus_uid": f"p{i}", "condition": "probe", "min_branch_um": mb,
             "length_density_um_per_um3": float(v), "detected": bool(v > 0)}
            for i, v in enumerate(pos_density)]
    rows += [{"nucleus_uid": f"n{i}", "condition": "no-probe", "min_branch_um": mb,
              "length_density_um_per_um3": float(v), "detected": bool(v > 0)}
             for i, v in enumerate(neg_density)]
    return summarise(pd.DataFrame(rows), "probe", "no-probe").iloc[0]


def test_a_working_detector_scores_detect_neg_zero():
    rng = np.random.default_rng(0)
    r = curve(rng.gamma(3, 0.08, 40), np.zeros(40))
    assert r.detect_neg == 0.0
    assert r.rank_biserial == pytest.approx(1.0)


def test_a_broken_detector_is_caught_even_though_the_medians_differ():
    """THE FAILURE THIS TEST EXISTS FOR. A 1.8x ratio with an overlapping distribution is not a
    working detector, and a ratio alone would not say so."""
    rng = np.random.default_rng(0)
    r = curve(rng.gamma(3, 0.08, 40), rng.gamma(3, 0.05, 40))
    assert r.detect_neg == 1.0, "every control nucleus reports antennas"
    assert r.ratio_med > 1.5, "and the medians still differ"
    assert r.rank_biserial < 0.8, "while the distributions overlap"


def test_rank_biserial_is_zero_for_identical_distributions():
    rng = np.random.default_rng(0)
    v = rng.gamma(3, 0.08, 60)
    rb, p = mannwhitney(v, v.copy())
    assert rb == pytest.approx(0.0, abs=1e-9)
    assert p > 0.9


def test_the_statistic_refuses_below_n_equals_three():
    rb, p = mannwhitney([1.0, 2.0], [3.0, 4.0])
    assert np.isnan(rb) and np.isnan(p)


def test_the_curve_has_one_row_per_min_branch():
    rng = np.random.default_rng(0)
    rows = []
    for mb in (0.25, 0.5, 1.0, 2.0):
        for i in range(10):
            rows.append({"nucleus_uid": f"p{i}", "condition": "probe", "min_branch_um": mb,
                         "length_density_um_per_um3": float(rng.gamma(3, 0.08)), "detected": True})
            rows.append({"nucleus_uid": f"n{i}", "condition": "no-probe", "min_branch_um": mb,
                         "length_density_um_per_um3": 0.0, "detected": False})
    s = summarise(pd.DataFrame(rows), "probe", "no-probe")
    assert len(s) == 4
    assert list(s.min_branch_um) == sorted(s.min_branch_um), "reported as a CURVE, in order"


def test_it_refuses_an_output_tree_with_no_work_artefacts(tmp_path):
    pd.DataFrame([{"nucleus_uid": "a", "condition": "probe"}]).to_csv(
        tmp_path / "antenna_nuclei.csv", index=False)
    with pytest.raises(FileNotFoundError, match="keep_work"):
        run_acceptance(tmp_path, "probe", "no-probe")


def test_it_refuses_a_condition_that_is_not_there(tmp_path):
    (tmp_path / "work").mkdir()
    pd.DataFrame([{"nucleus_uid": "a", "condition": "probe"}]).to_csv(
        tmp_path / "antenna_nuclei.csv", index=False)
    with pytest.raises(ValueError, match="not in antenna_nuclei"):
        run_acceptance(tmp_path, "probe", "typo-condition")


def test_it_refuses_to_report_a_verdict_with_no_control(tmp_path):
    """Refusing beats reporting a test that was not performed."""
    (tmp_path / "work").mkdir()
    pd.DataFrame([{"nucleus_uid": "a", "condition": "probe"},
                  {"nucleus_uid": "b", "condition": "no-probe"}]).to_csv(
        tmp_path / "antenna_nuclei.csv", index=False)
    df = pd.read_csv(tmp_path / "antenna_nuclei.csv")
    df = df[df.condition == "probe"]
    df.to_csv(tmp_path / "antenna_nuclei.csv", index=False)
    with pytest.raises(ValueError, match="not in antenna_nuclei"):
        run_acceptance(tmp_path, "probe", "no-probe")
