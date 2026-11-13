"""spots3d: the MATLAB port, piece by piece, then end to end on a phantom.

No .nd2 and no real data. Every expected value is derived from what the MATLAB
code does (line numbers in docs/spots3d.md), not from what the port returns.
"""
import math

import numpy as np
import pytest
import tifffile
from scipy import ndimage as ndi
from skimage.morphology import disk

from nucleus3d.core.export import export_field_labels
from nucleus3d.core.io import Stack
from spots3d import ChannelSpec, SpotParams, detect_field
from spots3d.detect import (components_gated, dbscan_merge, erode_mask_union,
                            erosion_radius_px, matlab_round, threshold_in_nucleus)
from spots3d.features import elongation_matlab, spot_intensities
from spots3d.filters import dog_matlab, imgaussfilt2d
from spots3d.io import read_labels

VOXEL = (0.2, 0.0459985, 0.0459985)


# --------------------------------------------------------------------- filters

def test_gaussian_kernel_radius_is_matlab_ceil_2sigma():
    vol = np.zeros((3, 41, 41), np.float32)
    vol[1, 20, 20] = 1.0
    sigma = 1.52
    out = imgaussfilt2d(vol, sigma)
    r = math.ceil(2 * sigma)                       # MATLAB filter size 2*ceil(2s)+1
    nz = np.nonzero(out[1, 20])[0]
    assert nz.min() == 20 - r and nz.max() == 20 + r
    assert out[0].max() == 0 and out[2].max() == 0          # no smoothing along z
    assert out.sum() == pytest.approx(1.0, rel=1e-6)        # normalised kernel


def test_gaussian_pads_by_replication():
    # a step at the image edge: 'replicate' keeps the edge column at its own value
    vol = np.zeros((1, 9, 30), np.float32)
    vol[..., :3] = 10.0
    out = imgaussfilt2d(vol, 1.0)
    assert out[0, 4, 0] > 9.9                   # reflect would also give this, so:
    const = imgaussfilt2d(np.full((1, 9, 9), 7.0, np.float32), 2.0)
    assert np.allclose(const, 7.0)


def test_dog_small_sigma_branch():
    rng = np.random.default_rng(0)
    vol = rng.normal(100, 10, (2, 40, 40)).astype(np.float32)
    px = 0.046
    # 1e-4 um / px < 0.3 -> I - G(large), as for MATLAB's S5P channel
    d = dog_matlab(vol, 1e-4, 0.5, px)
    assert np.allclose(d, vol - imgaussfilt2d(vol, 0.5 / px), atol=1e-4)
    d2 = dog_matlab(vol, 0.07, 0.3, px)
    assert np.allclose(d2, imgaussfilt2d(vol, 0.07 / px) - imgaussfilt2d(vol, 0.3 / px),
                       atol=1e-4)


# --------------------------------------------------------------------- erosion

def test_matlab_round_and_radius():
    assert matlab_round(2.5) == 3 and matlab_round(3.5) == 4      # numpy: 2 and 4
    assert erosion_radius_px(0.5, VOXEL[2]) == 11                  # 10.87 -> 11


def test_erosion_equals_disk_erosion_with_foreground_border():
    rng = np.random.default_rng(1)
    m = ndi.binary_dilation(rng.random((2, 80, 90)) > 0.995, iterations=12)
    m[:, :, :5] = True                                   # touches the image edge
    r = 6
    got = erode_mask_union(m, r)
    for z in range(m.shape[0]):
        ref = ndi.binary_erosion(m[z], structure=disk(r), border_value=1)
        assert np.array_equal(got[z], ref)


# --------------------------------------------------------------- thresholding

def test_threshold_mean_plus_k_sd_ddof1_and_outside_as_zero():
    dog = np.array([[[1.0, 2.0, 3.0, 10.0, 5.0]]], np.float32)
    mask = np.array([[[True, True, True, True, False]]])
    binary, mu, sd, thr = threshold_in_nucleus(dog, mask, 1.0)
    vals = np.array([1, 2, 3, 10.0])
    assert mu == pytest.approx(vals.mean())
    assert sd == pytest.approx(vals.std(ddof=1))            # MATLAB std is N-1
    assert binary.tolist() == [[[False, False, False, True, False]]]  # 5.0 is outside
    # negative threshold: masked-out voxels are 0 > thr and pass, as in MATLAB
    b2, *_ = threshold_in_nucleus(dog - 100, mask, 0.0)
    assert b2[0, 0, 4]


# ---------------------------------------------------------------- components

def test_volume_gate_is_strict():
    b = np.zeros((1, 3, 40), bool)
    b[0, 1, 1:11] = True        # 10 voxels -- equal to the limit: dropped
    b[0, 1, 20:31] = True       # 11 voxels: kept
    lab, n = components_gated(b, (1.0, 1.0, 1.0), 10.0)
    assert n == 2 and lab.max() == 1 and (lab > 0).sum() == 11


def test_connectivity_is_18():
    b = np.zeros((2, 2, 2), bool)
    b[0, 0, 0] = b[1, 1, 1] = True          # corner contact only: 26- not 18-connected
    assert components_gated(b, (1, 1, 1), 0.0)[0].max() == 2
    b[1, 1, 1] = False
    b[1, 1, 0] = True                       # edge contact: 18-connected
    assert components_gated(b, (1, 1, 1), 0.0)[0].max() == 1


def test_component_order_is_matlab_column_major():
    b = np.zeros((1, 10, 10), bool)
    b[0, 5, 0] = True          # A: x=0, y=5  -> MATLAB linear index 5 (first)
    b[0, 0, 5] = True          # B: x=5, y=0  -> 50; numpy raster order puts B first
    lab, _ = components_gated(b, (1, 1, 1), 0.0)
    assert lab[0, 5, 0] == 1 and lab[0, 0, 5] == 2


def test_dbscan_chains_and_unions():
    lab = np.zeros((1, 1, 40), np.int32)
    lab[0, 0, [0, 10, 20, 39]] = [1, 2, 3, 4]   # spacing 10, 10, 19 (um with voxel 1)
    merged, n_merged = dbscan_merge(lab, (1, 1, 1), 10.0)
    # 1-2 and 2-3 within eps (<=, rangesearch), 1-3 is not: still one cluster
    assert merged[0, 0, [0, 10, 20]].tolist() == [1, 1, 1]
    assert merged[0, 0, 39] == 2
    assert n_merged.tolist() == [3, 1]
    off, _ = dbscan_merge(lab, (1, 1, 1), 0.0)
    assert off.max() == 4


def test_volume_gate_runs_before_merge():
    b = np.zeros((1, 3, 30), bool)
    b[0, 1, 0:12] = True      # 12 voxels, kept
    b[0, 1, 14:16] = True     # 2 voxels, 3 um away -- would merge, but is gated first
    lab, _ = components_gated(b, (1, 1, 1), 5.0)
    merged, n_merged = dbscan_merge(lab, (1, 1, 1), 20.0)
    assert merged.max() == 1 and (merged > 0).sum() == 12 and n_merged.tolist() == [1]


# ------------------------------------------------------------------- features

def test_elongation_of_rectangle_is_aspect_ratio():
    # with MATLAB's +1/12 pixel moment a w x h rectangle has axes 2w/sqrt3, 2h/sqrt3
    obj = np.zeros((3, 4, 12), bool)
    obj[:, 1:3, 1:11] = True
    assert elongation_matlab(obj) == pytest.approx(5.0, rel=1e-12)
    obj[1] = False                      # ceil(3/2) = plane 2 (1-based) is empty
    assert math.isnan(elongation_matlab(obj))


def test_intensity_ratio_is_median_over_nucleus_median():
    img = np.full((2, 10, 10), 100.0, np.float32)
    lab = np.zeros(img.shape, np.int32)
    lab[:, 4:6, 4:6] = 1
    img[lab == 1] = 300.0
    mask = np.ones(img.shape, bool)
    ratio, nuc_med = spot_intensities(lab, img, mask, 0.0, 0.046)
    assert nuc_med == 100.0 and ratio.tolist() == [3.0]


# ---------------------------------------------------------------- label input

def _stack(shape=(4, 32, 32), voxel=VOXEL):
    return Stack(data=np.zeros((shape[0], 1) + shape[1:], np.uint16), channels=["Dapi"],
                 voxel_um=voxel, name="f.nd2", source_path="/d/f.nd2", position=0)


def test_reads_nucleus3d_labels_and_refuses_other_grids(tmp_path):
    lab = np.zeros((4, 32, 32), np.uint16)
    lab[1:3, 5:20, 5:20] = 1
    p = export_field_labels(_stack(), lab, tmp_path)
    assert np.array_equal(read_labels(p, lab.shape, VOXEL), lab)
    with pytest.raises(ValueError, match="voxel size"):
        read_labels(p, lab.shape, (0.2, 0.092, 0.092))
    with pytest.raises(ValueError, match="shape"):
        read_labels(p, (4, 64, 64), VOXEL)
    raw = tmp_path / "uncalibrated.tif"
    tifffile.imwrite(raw, lab, photometric="minisblack")
    with pytest.raises(ValueError, match="PhysicalSize"):
        read_labels(raw, lab.shape, VOXEL)


# ------------------------------------------------------------------ end to end

def _gauss_spot(vol, zyx_px, sig_px, amp):
    z, y, x = np.indices(vol.shape, dtype=np.float32)
    vol += amp * np.exp(-0.5 * (((z - zyx_px[0]) / sig_px[0]) ** 2
                                + ((y - zyx_px[1]) / sig_px[1]) ** 2
                                + ((x - zyx_px[2]) / sig_px[2]) ** 2))


def _phantom():
    rng = np.random.default_rng(3)
    shape = (16, 120, 220)
    labels = np.zeros(shape, np.uint16)
    labels[:, 10:110, 10:100] = 1
    labels[:, 10:110, 120:210] = 2
    sig = (0.3 / 0.2, 0.1 / 0.046, 0.1 / 0.046)     # ~diffraction-limited spot
    truth = {1: [(8, 40, 40), (7, 80, 70)],
             2: [(8, 30, 140), (8, 30, 190), (9, 60, 160), (7, 90, 140), (8, 90, 190)]}
    ch = rng.normal(300, 10, shape).astype(np.float32)
    for pts in truth.values():
        for p in pts:
            _gauss_spot(ch, p, sig, 2000.0)
    edge_spot = (8, 60, 15)            # 5 px inside nucleus 1: inside the 11 px erosion rim
    _gauss_spot(ch, edge_spot, sig, 2000.0)
    other = rng.normal(300, 10, shape).astype(np.float32)
    return labels, {"A": ch, "B": other}, truth, edge_spot


def test_end_to_end_recovers_and_assigns_spots():
    labels, images, truth, edge_spot = _phantom()
    params = SpotParams(channels=(ChannelSpec("A", index=0), ChannelSpec("B", index=1)))
    res = detect_field(images, labels, VOXEL, params, "/d/field_x.nd2", 0)
    spots = res["spots"]
    a = spots[spots["channel"] == "A"]
    assert (spots["channel"] == "B").sum() == 0               # noise-only channel

    for lab_id, pts in truth.items():
        got = a[a["label"] == lab_id]
        assert len(got) == len(pts)
        assert set(got["nucleus_uid"]) == {f"field_x_p00_nuc{lab_id:03d}"}
        c = got[["centroid_z_px", "centroid_y_px", "centroid_x_px"]].to_numpy()
        for p in pts:
            assert np.abs(c - np.array(p)).max(axis=1).min() < 0.5
    # the edge spot lies in the eroded-away rim: excluded at 0.5 um, found at 0
    assert len(a) == 7
    res0 = detect_field(images, labels, VOXEL, params.with_(erosion_um=0.0),
                        "/d/field_x.nd2", 0)
    assert (res0["spots"]["channel"] == "A").sum() == 8

    nuc = res["nuclei"].set_index(["label", "channel"])
    assert not nuc.loc[(1, "A"), "qc_more_than_max_spots"]
    assert nuc.loc[(2, "A"), "qc_more_than_max_spots"]          # 5 > 4: flagged...
    assert nuc.loc[(2, "A"), "n_spots"] == 5                    # ...not removed
    assert (a["int_ratio_A"] > 2).all()
    assert res["spot_volumes"]["A"].max() == 7
