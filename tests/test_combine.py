"""The join between the two tables.

The failure this guards against is specific: a mismatched `nucleus_uid`
produces an empty merge, and an empty merge is indistinguishable from a
dataset in which nothing was detected.
"""
import pandas as pd
import pytest

from combine import join_nuclei_and_antennas

UIDS = [f"SetC_Control_004_crop_p00_nuc{i:03d}" for i in range(1, 6)]


def write_pair(tmp_path, nucleus_uids=UIDS, antenna_uids=UIDS[:3]):
    nd, ad = tmp_path / "nuclei", tmp_path / "antennas"
    nd.mkdir(), ad.mkdir()
    pd.DataFrame({
        "nucleus_uid": nucleus_uids,
        "position": [0] * len(nucleus_uids),
        "label": list(range(1, len(nucleus_uids) + 1)),
        "file": ["SetC_Control_004_crop.nd2"] * len(nucleus_uids),
        "condition": ["Control"] * len(nucleus_uids),
        "volume_um3": [200.0] * len(nucleus_uids),
        "DAPI_median": [1000.0] * len(nucleus_uids),
    }).to_csv(nd / "nuclei_measurements.csv", index=False)
    pd.DataFrame({
        "nucleus_uid": antenna_uids,
        "position": [0] * len(antenna_uids),
        "label": list(range(1, len(antenna_uids) + 1)),
        "file": ["SetC_Control_004_crop"] * len(antenna_uids),
        "condition": ["Control"] * len(antenna_uids),
        "total_antenna_length_um": [12.0] * len(antenna_uids),
        "length_density_um_per_um3": [0.06] * len(antenna_uids),
    }).to_csv(ad / "antenna_nuclei.csv", index=False)
    return nd, ad


def test_it_joins_on_nucleus_uid(tmp_path):
    nd, ad = write_pair(tmp_path)
    both = join_nuclei_and_antennas(nd, ad, log=lambda *a: None)

    assert len(both) == 3
    assert {"DAPI_median", "length_density_um_per_um3"} <= set(both.columns)
    assert both.attrs["join"] == {"n_nuclei": 5, "n_antenna": 3, "n_joined": 3,
                                  "n_nucleus_only": 2, "n_antenna_only": 0}


def test_columns_that_must_agree_are_kept_once(tmp_path):
    nd, ad = write_pair(tmp_path)
    both = join_nuclei_and_antennas(nd, ad, log=lambda *a: None)

    assert "position" in both.columns and "position_nucleus" not in both.columns
    # `file` means different things on the two sides and must NOT be collapsed.
    assert {"file_nucleus", "file_antenna"} <= set(both.columns)


def test_an_empty_overlap_raises_rather_than_returning_nothing(tmp_path):
    """The whole reason this module exists."""
    nd, ad = write_pair(tmp_path, antenna_uids=["SetC_Control_004_crop#p00_nuc001"])
    with pytest.raises(ValueError, match="present in both tables"):
        join_nuclei_and_antennas(nd, ad, log=lambda *a: None)


def test_disagreeing_labels_under_one_uid_raise(tmp_path):
    nd, ad = write_pair(tmp_path)
    a = pd.read_csv(ad / "antenna_nuclei.csv")
    a.loc[0, "label"] = 99
    a.to_csv(ad / "antenna_nuclei.csv", index=False)
    with pytest.raises(ValueError, match="disagree on `label`"):
        join_nuclei_and_antennas(nd, ad, log=lambda *a: None)


def test_a_duplicated_uid_raises_rather_than_multiplying_rows(tmp_path):
    nd, ad = write_pair(tmp_path, nucleus_uids=UIDS[:4] + [UIDS[0]])
    with pytest.raises(ValueError, match="duplicated nucleus_uid"):
        join_nuclei_and_antennas(nd, ad, log=lambda *a: None)


def test_a_missing_table_says_which_pipeline_to_run(tmp_path):
    nd, ad = write_pair(tmp_path)
    (ad / "antenna_nuclei.csv").unlink()
    with pytest.raises(FileNotFoundError, match="run_antennas.py"):
        join_nuclei_and_antennas(nd, ad, log=lambda *a: None)


def test_failures_csv_is_reported(tmp_path):
    nd, ad = write_pair(tmp_path)
    pd.DataFrame({"field": ["x"], "error": ["boom"]}).to_csv(
        ad / "failures.csv", index=False)
    lines = []
    join_nuclei_and_antennas(nd, ad, log=lines.append)
    assert any("failures.csv" in ln and "1 field(s) failed" in ln for ln in lines)
