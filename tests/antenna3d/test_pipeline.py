"""`run_field` end to end on a synthetic field, with no .nd2 involved.

Covers what the per-module tests cannot: crop planning and the plane cache, the label TIFF round
trip, the work artefact the acceptance test reads, and the shape of the output tree.
"""
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile

from antenna3d import pipeline
from antenna3d.grid import Grid
from antenna3d.io import FieldRef
from antenna3d.params import Params

import phantom

SHAPE = (30, 200, 200)
BIN = 5


class FakeField:
    """Just enough Field for `run_field`: planes on demand, a grid, and an identity."""

    def __init__(self, vol, grid, stem, condition):
        self._vol, self.grid = vol, grid
        self.stem, self.position, self.condition = stem, 0, condition
        self.path = Path(f"/nowhere/{stem}.nd2")
        self.n_positions = 1

    @property
    def field_id(self):
        return f"{self.stem}_p{self.position:02d}"

    @property
    def shape(self):
        return self._vol.shape

    @property
    def slab_um(self):
        return self.shape[0] * self.grid.dz_um

    def plane(self, z, role="actin"):
        return self._vol[int(z)]

    def volume(self, role="actin", **kw):
        return self._vol


def build_field(tmp_path, stem="synthetic_field", condition="probe", seed=0):
    """Two spherical nuclei, each with one oblique filament, plus a label TIFF for them."""
    grid = phantom.grid_for()
    segs, mask = [], np.zeros(SHAPE, bool)
    labels = np.zeros((SHAPE[0], SHAPE[1] // BIN, SHAPE[2] // BIN), np.uint16)
    lgrid = Grid(grid.dz_um, grid.dy_um * BIN, grid.dx_um * BIN)
    for i, cx_um in enumerate((2.3, 6.9), start=1):
        c = (SHAPE[0] * grid.dz_um / 2, 4.6, cx_um)
        r = 2.1
        m = phantom.sphere_mask(SHAPE, grid, c, r)
        mask |= m
        labels[phantom.sphere_mask(labels.shape, lgrid, c, r)] = i
        segs.append(((c[0] - 0.5, c[1] - 1.0, c[2] - 1.0), (c[0] + 0.5, c[1] + 1.0, c[2] + 1.0)))
    vol = phantom.render(SHAPE, grid, segs, amplitude=3000.0, noise_sd=300.0, seed=seed)
    lp = tmp_path / "labels" / f"{stem}_p00.tif"
    lp.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(lp, labels, ome=True,
                     metadata={"axes": "ZYX", "PhysicalSizeZ": lgrid.dz_um,
                               "PhysicalSizeY": lgrid.dy_um, "PhysicalSizeX": lgrid.dx_um})
    return FakeField(vol, grid, stem, condition), segs


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("pipeline")
    field, segs = build_field(tmp)
    p = Params()
    # The synthetic nuclei are ~39 um3, well under the shipped 100 um3 gate.
    p = p.with_(segment=replace(p.segment, min_volume_um3=10.0, drop_xy_border_touching=False),
                graph=replace(p.graph, measure_width=False))

    from contextlib import contextmanager

    @contextmanager
    def fake_open(ref, *a, **kw):
        yield field

    real = pipeline.open_field
    pipeline.open_field = fake_open
    try:
        ref = FieldRef(field.path, 0, 1, field.condition)
        out = tmp / "out"
        r = pipeline.run_field(ref, p, out, labels_dir=tmp / "labels", keep_work="small",
                               qc=True, log=lambda *a: None)
    finally:
        pipeline.open_field = real
    return r, out, segs, p


def test_it_produced_a_graph_per_nucleus(ran):
    r, out, segs, _ = ran
    assert len(r["rows"]) == 2, "two nuclei in, two rows out"
    assert len(list((out / "graphs").glob("*.graphml"))) == 2


def test_the_output_tree_has_what_the_readme_says(ran):
    _, out, _, _ = ran
    assert (out / "graphs").is_dir() and (out / "work").is_dir()
    assert len(list((out / "qc" / "traces").glob("*.png"))) == 2
    assert (out / "qc" / "synthetic_field_p00.png").exists()


def test_the_uid_is_the_join_key_everywhere(ran):
    r, out, _, _ = ran
    for row in r["rows"]:
        uid = row["nucleus_uid"]
        assert uid.startswith("synthetic_field_p00_nuc")
        assert (out / "graphs" / f"{uid}.graphml").exists()
        assert (out / "work" / f"{uid}.npz").exists()
        assert (out / "qc" / "traces" / f"{uid}.png").exists()
    assert {e["nucleus_uid"] for e in r["edges"]} <= {q["nucleus_uid"] for q in r["rows"]}


def test_the_work_artefact_round_trips_and_is_small(ran):
    _, out, _, _ = ran
    for p in (out / "work").glob("*.npz"):
        w = pipeline.load_work(p)
        assert w["fg"].shape == w["mask"].shape and w["fg"].dtype == bool
        assert w["null_level"] > 0 and w["noise_sd"] > 0
        assert "actin_flat" not in w, "keep_work='small' must not store the image"
        assert p.stat().st_size < 400_000, "the small artefact must stay small"


def test_the_filaments_are_recovered_at_about_the_right_length(ran):
    r, _, segs, _ = ran
    true = sum(phantom.seg_length(a, b) for a, b in segs)
    got = sum(row["total_antenna_length_um"] for row in r["rows"])
    assert got == pytest.approx(true, rel=0.35), f"{got:.2f} um traced for {true:.2f} um drawn"


def test_the_field_summary_reports_what_it_promises(ran):
    r, _, _, _ = ran
    s = r["summary"]
    for k in ("field_id", "n_nuclei", "frac_detected", "separation_med", "noise_sd_med",
              "length_density_med", "axial_edge_fraction_med", "n_puncta_rejected"):
        assert k in s, k
    assert s["n_nuclei"] == 2 and s["condition"] == "probe"


def test_a_missing_label_volume_says_where_it_looked(tmp_path):
    p = Params()
    ref = FieldRef(Path("/nowhere/other_field.nd2"), 0, 1, None)
    with pytest.raises(FileNotFoundError) as e:
        pipeline._labels_for(ref, tmp_path, p, None, lambda *a: None)
    assert "other_field_p00.tif" in str(e.value)


def test_a_label_volume_with_the_wrong_depth_is_refused(tmp_path):
    """Masks from another stack would be measured against this one's pixels without error."""
    field, _ = build_field(tmp_path)
    lp = tmp_path / "labels" / "synthetic_field_p00.tif"
    lab = tifffile.imread(lp)
    tifffile.imwrite(lp, lab[:-3], ome=True,
                     metadata={"axes": "ZYX", "PhysicalSizeZ": 0.2,
                               "PhysicalSizeY": 0.23, "PhysicalSizeX": 0.23})
    from contextlib import contextmanager

    @contextmanager
    def fake_open(ref, *a, **kw):
        yield field

    real = pipeline.open_field
    pipeline.open_field = fake_open
    try:
        with pytest.raises(ValueError, match="same stack"):
            pipeline.run_field(FieldRef(field.path, 0, 1, None), Params(), tmp_path / "o",
                               labels_dir=tmp_path / "labels", log=lambda *a: None)
    finally:
        pipeline.open_field = real


def test_a_label_tiff_without_a_voxel_size_is_refused(tmp_path):
    """A label volume with no spacing is not a measurement, and guessing rescales every volume."""
    from antenna3d.segment import labels_from_tiff
    p = tmp_path / "bare.tif"
    tifffile.imwrite(p, np.zeros((4, 8, 8), np.uint16))
    with pytest.raises(ValueError, match="voxel_um"):
        labels_from_tiff(p)
    labels, g = labels_from_tiff(p, voxel_um=(0.2, 0.23, 0.23))
    assert labels.shape == (4, 8, 8) and g.dz_um == 0.2
