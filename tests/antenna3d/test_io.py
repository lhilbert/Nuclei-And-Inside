"""Channel identity, optical regime, and the position axis. Each guards a silent wrong answer."""
import os
from pathlib import Path

import pytest

from antenna3d.io import (IngestError, check_sampling, condition_of, describe_file,
                          n_positions, nucleus_uid, resolve_channels)
from antenna3d.params import Channels, Optics

ND2 = os.environ.get("ANTENNA3D_TEST_ND2")


def header(**kw):
    h = {"stem": "130726_100xHAKactin-488_sphere_postfix_002_crop",
         "axes": "ZCYX", "sizes": {"Z": 123, "C": 3, "Y": 2280, "X": 2588},
         "dxy_um": 0.04599853515625, "dy_um": 0.04599853515625, "dz_um": 0.2,
         "channels_stored": ["DAPI", "GFP", "RFP"], "na": 1.49}
    h.update(kw)
    return h


# ------------------------------------------------------------------ the channel is not an index

def test_channels_resolve_by_name():
    assert resolve_channels(header(), Channels(dna="DAPI", actin="GFP")) == {"dna": 0, "actin": 1}


def test_the_same_names_in_another_order_give_another_index():
    """The specific instance: actin at stored index 2 in one dataset and 0 in another."""
    a = resolve_channels(header(channels_stored=["DAPI", "Cy5", "GFP", "Black"]), Channels())
    b = resolve_channels(header(channels_stored=["GFP", "Black", "DAPI", "Cy5"]), Channels())
    assert a["actin"] == 2 and b["actin"] == 0
    assert a["dna"] == 0 and b["dna"] == 2


def test_a_missing_channel_raises_and_says_what_is_there():
    with pytest.raises(IngestError) as e:
        resolve_channels(header(channels_stored=["DAPI", "Cy5"]), Channels(actin="GFP"))
    assert "GFP" in str(e.value) and "describe_file" in str(e.value)


def test_an_ambiguous_channel_name_raises():
    with pytest.raises(IngestError, match="ambiguous"):
        resolve_channels(header(channels_stored=["DAPI", "GFP", "GFP"]), Channels())


# ------------------------------------------------------- the regime is the voxel size, not a name

def test_the_wrong_sampling_is_refused():
    with pytest.raises(IngestError) as e:
        check_sampling(header(dxy_um=0.114996337890625), Optics())
    assert "wrong physical scale" in str(e.value)


def test_the_wrong_axial_step_is_refused():
    with pytest.raises(IngestError, match="axial step"):
        check_sampling(header(dz_um=0.5), Optics())


def test_a_misnamed_file_fails_loudly():
    with pytest.raises(IngestError, match="misnamed at acquisition"):
        check_sampling(header(stem="130726_HAKactin_sphere_postfix_002_crop"), Optics())


def test_the_marker_can_be_switched_off():
    check_sampling(header(stem="anything at all"), Optics(stem_marker=None))


def test_a_disagreeing_na_is_refused():
    with pytest.raises(IngestError, match="NA"):
        check_sampling(header(na=1.15), Optics())


# ------------------------------------------------------------------- one .nd2 is many fields

def test_an_undeclared_position_axis_is_refused():
    h = header(axes="PZCYX", sizes={"P": 34, "Z": 27, "C": 4, "Y": 1036, "X": 1430})
    with pytest.raises(IngestError, match="multi_series"):
        n_positions(h, multi_series=False)
    assert n_positions(h, multi_series=True) == 34


def test_a_position_axis_inside_the_z_loop_is_refused():
    """`position * n_planes + z` addresses the wrong frame if P does not enclose Z."""
    h = header(axes="ZPCYX", sizes={"P": 34, "Z": 27, "C": 4, "Y": 1036, "X": 1430})
    with pytest.raises(IngestError, match="outside"):
        n_positions(h, multi_series=True)


def test_a_single_field_needs_no_declaration():
    assert n_positions(header(), multi_series=False) == 1


# ------------------------------------------------------------------------------- identifiers

def test_nucleus_uid_matches_the_nucleus3d_format():
    assert nucleus_uid("SetC_Control_004_crop", 0, 3) == "SetC_Control_004_crop_p00_nuc003"


def test_nucleus_uid_is_built_forward_and_survives_a_doubled_separator():
    """A stem can contain the separator. That is why nothing ever parses an id out of a path."""
    uid = nucleus_uid("130726_100xHAKactin-488_sphere_postfix__crop", 0, 30)
    assert uid.endswith("_p00_nuc030")
    assert "__crop" in uid, "the stem is carried verbatim, not normalised away"


def test_nucleus_uid_is_filename_safe():
    assert nucleus_uid("a/b:c d", 1, 2) == "a-b-c-d_p01_nuc002"


def test_condition_comes_from_the_subfolder():
    root = Path("/data/exp")
    assert condition_of(root / "Control" / "f.nd2", root) == "Control"
    assert condition_of(root / "f.nd2", root) is None
    # NOTHING is parsed out of it: a trailing token is not automatically a drug.
    assert condition_of(root / "LatB_009" / "f.nd2", root) == "LatB_009"


def test_condition_survives_a_symlinked_input_tree(tmp_path):
    """Staging a subset as symlinks is normal. Resolving the file jumps out of the tree to
    wherever the data really lives, which silently returned None for every condition."""
    real = tmp_path / "store"
    real.mkdir()
    (real / "a.nd2").write_bytes(b"")
    root = tmp_path / "staged"
    (root / "Control").mkdir(parents=True)
    (root / "Control" / "a.nd2").symlink_to(real / "a.nd2")
    assert condition_of(root / "Control" / "a.nd2", root) == "Control"


# ------------------------------------------------------------------------------ against real data

@pytest.mark.needs_data
@pytest.mark.skipif(not ND2, reason="set ANTENNA3D_TEST_ND2 to an .nd2 file")
def test_describe_file_on_a_real_stack():
    out = describe_file(ND2)
    assert "channels" in out and "voxel (z,y,x)" in out and "positions" in out
    assert "never the indices" in out
