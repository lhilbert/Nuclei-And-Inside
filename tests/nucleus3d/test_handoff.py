"""The nucleus3d -> antenna3d handoff.

Two pipelines in one toolbox are only one toolbox if the masks travel and the
tables join. Everything here tests the seam between them, not either side.

No .nd2 and no real data: a Stack is a plain dataclass, so a synthetic one
exercises the whole handoff.
"""
import numpy as np
import pytest
from dataclasses import replace

from antenna3d.grid import Grid
from antenna3d.io import FieldRef, nucleus_uid as antenna_uid
from antenna3d.params import Params
from antenna3d.segment import labels_from_tiff, nuclei_from_labels
from nucleus3d.core.export import export_field_labels
from nucleus3d.core.io import Stack
from nucleus3d.core.quantify import nucleus_uid as nucleus_uid_of
from nucleus3d.core.segment import SegParams, segment_nuclei

# A stem with a DOUBLED separator. This shape of name exists in the project's own
# data and is exactly what broke an earlier id scheme that tried to parse itself
# back out of a path.
STEM = "130726_100xHAKactin-488_S5P-647_sphere_postfix__crop"
VOXEL = (0.2, 0.046, 0.046)


def make_stack(labels_shape=(14, 48, 48), position=3, voxel=VOXEL, stem=STEM):
    z, y, x = labels_shape
    data = np.zeros((z, 2, y, x), np.uint16)
    return Stack(data=data, channels=["DAPI", "GFP"], voxel_um=voxel,
                 name=f"{stem}.nd2#p{position:02d}", source_path=f"/data/{stem}.nd2",
                 position=position)


def two_nuclei(shape=(14, 48, 48)):
    lab = np.zeros(shape, np.int32)
    lab[3:11, 6:20, 6:20] = 1
    lab[4:10, 28:42, 28:42] = 2
    return lab


# --------------------------------------------------------------- the join key

def test_the_two_packages_build_the_same_nucleus_uid():
    """The join key is the whole point; if it disagrees the merge is empty."""
    stack = make_stack()
    assert nucleus_uid_of(stack, 7) == antenna_uid(STEM, 3, 7)


def test_the_uid_survives_a_stem_with_a_doubled_separator():
    assert nucleus_uid_of(make_stack(), 7).endswith("_crop_p03_nuc007")


def test_the_uid_is_safe_as_a_filename():
    """It is used verbatim as a .ome.tif and a .graphml name."""
    uid = nucleus_uid_of(make_stack(stem="odd name#with/chars"), 1)
    assert not (set(uid) & set(' #/\\?%*:|"<>'))


# ------------------------------------------------------- the file that travels

def test_the_label_file_is_named_where_antenna3d_looks_for_it(tmp_path):
    """antenna3d resolves <stem>_p<NN>.tif from its FieldRef, and nothing else."""
    stack = make_stack()
    written = export_field_labels(stack, two_nuclei(), tmp_path)
    ref = FieldRef(tmp_path / f"{STEM}.nd2", 3, 4, None)
    assert written.endswith(f"{ref.field_id}.tif")


def test_the_labels_round_trip_with_their_spacing(tmp_path):
    """A label volume without its spacing is not a measurement; antenna3d refuses
    one rather than guessing, so the spacing has to be in the file."""
    stack = make_stack()
    lab = two_nuclei()
    back, grid = labels_from_tiff(export_field_labels(stack, lab, tmp_path))

    assert np.array_equal(back, lab)
    assert (grid.dz_um, grid.dy_um, grid.dx_um) == pytest.approx(VOXEL)


def test_z_is_never_cropped(tmp_path):
    """antenna3d refuses labels whose plane count differs from the image's -- masks
    from a different stack would otherwise be measured against these pixels."""
    stack = make_stack()
    back, _ = labels_from_tiff(export_field_labels(stack, two_nuclei(), tmp_path))
    assert back.shape[0] == stack.data.shape[0]


def test_more_labels_than_uint16_widen_rather_than_wrap(tmp_path):
    """A wrapped label silently merges two nuclei into one."""
    lab = np.zeros((2, 4, 4), np.int32)
    lab[0, 0, 0] = 70000
    back, _ = labels_from_tiff(
        export_field_labels(make_stack(labels_shape=(2, 4, 4)), lab, tmp_path))
    assert back.max() == 70000


# ------------------------------------------- the masks are usable at the far end

@pytest.fixture(scope="module")
def segmented():
    """Three well-separated spheres, fully contained in z, segmented for real."""
    dz, dy, dx = 0.4, 0.2, 0.2
    shape = (40, 120, 120)
    zz, yy, xx = np.mgrid[0:shape[0], 0:shape[1], 0:shape[2]].astype(np.float32)
    vol = np.zeros(shape, np.float32)
    for cz, cy, cx in [(20, 30, 30), (20, 30, 85), (20, 85, 55)]:
        r2 = (((zz - cz) * dz) ** 2 + ((yy - cy) * dy) ** 2 + ((xx - cx) * dx) ** 2)
        vol += 3000.0 * np.exp(-r2 / (2 * 2.2 ** 2))
    rng = np.random.default_rng(0)
    vol += rng.normal(100.0, 20.0, shape).astype(np.float32)

    labels, props = segment_nuclei(vol, (dz, dy, dx),
                                   SegParams(sigma_large_um=6.0, min_volume_um3=15.0))
    return labels, props, (dz, dy, dx)


def test_segmentation_produces_three_nuclei(segmented):
    _, props, _ = segmented
    assert len(props) == 3


def test_the_masks_are_not_cut_in_z(segmented):
    """`axial_edge_fraction` is antenna3d's model-free covariate for a mask that
    stops mid-nucleus: 0 when it tapers away to nothing, 1 when it is cut at its
    widest. 2D-per-plane segmentation measured 0.8764 on data that is not axially
    truncated at all, and every density computed from it was divided by a volume
    ~23% too small. These nuclei are whole and this is a 3D segmentation, so the
    handed-over masks must taper. These score 0.087-0.088, against the 0.0548
    antenna3d measured for a 3D flow field on real data.
    """
    labels, _, voxel = segmented
    p = Params()
    p = p.with_(segment=replace(p.segment, min_volume_um3=10.0,
                                drop_xy_border_touching=False))
    nuclei = nuclei_from_labels(labels, Grid(*voxel), p)

    assert len(nuclei) == 3
    assert max(n.axial_edge_fraction for n in nuclei) < 0.3
    assert not any(n.z_truncated for n in nuclei)


def test_a_field_with_no_nuclei_still_gets_a_label_file(tmp_path):
    """An EMPTY label volume means "nothing here"; a MISSING one means the
    handoff broke. Omitting the file conflates the two, and makes the second
    pipeline record a field that simply had nothing in it as a failure.
    """
    stack = make_stack()
    empty = np.zeros(stack.data.shape[:1] + stack.data.shape[2:], np.int32)
    back, _ = labels_from_tiff(export_field_labels(stack, empty, tmp_path))
    assert back.shape == empty.shape and back.max() == 0


def test_save_labels_bypasses_the_measurement_cache():
    """The cache stores the TABLE, not the label image.

    A cached field returns early, before segmentation runs, so a cache hit
    with `save_labels` on would produce a measurement row and no label
    volume -- the handoff silently missing exactly the fields that were
    cheapest to redo. `save_labels` therefore joins `save_qc` and
    `save_boxes` in the bypass condition.
    """
    import inspect
    from nucleus3d.core import pipeline

    src = inspect.getsource(pipeline.process_field)
    assert "if cache_dir and not (save_qc or save_boxes or save_labels):" in src
    assert inspect.signature(pipeline.process_field).parameters[
        "save_labels"].default is True


def test_save_labels_reaches_the_parallel_worker():
    """run() dispatches to run_parallel for n_workers != 1, which rebuilds the
    call to process_field from a positional job tuple. A parameter added only
    to process_field would be silently dropped on that path."""
    import inspect
    from nucleus3d.core import parallel

    assert "save_labels" in inspect.signature(parallel.run_parallel).parameters
    src = inspect.getsource(parallel._process_one)
    assert "save_labels) = job" in src
    assert "save_labels=save_labels" in src
