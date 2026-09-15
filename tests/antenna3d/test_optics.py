"""The PSF, the ten scales derived from it, and the sampling that refuses to support them."""
import pytest

from antenna3d import optics
from antenna3d.optics import PSF, SamplingError
from antenna3d.params import Optics

import phantom

WIDEFIELD = phantom.REGIMES["widefield-100x"]
ISIM = phantom.REGIMES["isim-100x"]


def test_all_ten_scales_are_produced():
    s = optics.scales_for(PSF(**{k: WIDEFIELD[k] for k in ("lateral_fwhm_um", "axial_fwhm_um")}),
                          WIDEFIELD["grid"].dx_um, WIDEFIELD["grid"].dz_um)
    for k in optics.SCALE_KEYS:
        assert k in s and s[k] > 0
    assert len(optics.SCALE_KEYS) == 10


def test_work_grid_is_nyquist_per_axis_not_isotropic():
    psf = PSF(0.2302, 0.5923)
    s = optics.scales_for(psf, 0.04599853515625, 0.2)
    assert s["work_xy_um"] == pytest.approx(0.2302 / 2)
    assert s["work_z_um"] == pytest.approx(0.5923 / 2)
    assert s["work_z_um"] > s["work_xy_um"], "the work grid is NOT isotropic in microns"


def test_three_scales_follow_the_work_grid():
    psf = PSF(0.2302, 0.5923)
    s = optics.scales_for(psf, 0.04599853515625, 0.2)
    assert s["junction_merge_um"] == pytest.approx(3 * s["work_xy_um"])
    assert s["tangent_skip_um"] == pytest.approx(1 * s["work_z_um"])
    assert s["tangent_window_um"] == pytest.approx(3 * s["work_z_um"])


def test_isim_sampling_is_refused():
    """0.826 native pixels is not a Gaussian derivative, and running it anyway is silent."""
    psf = PSF(ISIM["lateral_fwhm_um"], ISIM["axial_fwhm_um"])
    with pytest.raises(SamplingError) as e:
        optics.scales_for(psf, ISIM["grid"].dx_um, ISIM["grid"].dz_um)
    msg = str(e.value)
    assert "0.8" in msg and "samples per FWHM" in msg
    # Non-strict reports it instead of raising, so a diagnostic can still see the numbers.
    s = optics.scales_for(psf, ISIM["grid"].dx_um, ISIM["grid"].dz_um, strict=False)
    assert s["sigma_lateral_px"] < optics.MIN_SIGMA_PX and s["notes"]


def test_both_axes_are_checked_not_only_the_lateral_one():
    """A filter that aliases in z aliases as badly as one that aliases in x."""
    psf = PSF(lateral_fwhm_um=1.0, axial_fwhm_um=0.4)     # lateral fine, axial not
    with pytest.raises(SamplingError) as e:
        optics.scales_for(psf, 0.05, 0.5)
    assert "axial" in str(e.value)


def test_widefield_sampling_is_accepted():
    psf = PSF(0.2302, 0.5923)
    s = optics.scales_for(psf, 0.04599853515625, 0.2)
    assert not s["notes"]
    assert s["samples_per_lateral_fwhm"] == pytest.approx(5.0, abs=0.01)


def test_frozen_block_is_all_or_nothing():
    with pytest.raises(ValueError, match="all-or-nothing"):
        Optics(frozen_scales={"sigma_um": 0.078}).scales()
    with pytest.raises(ValueError, match="not derived scales"):
        Optics(frozen_scales={**optics.E1_100X_FROZEN, "nope": 1.0}).scales()


def test_a_freeze_stays_visible():
    o = Optics(frozen_scales=optics.E1_100X_FROZEN)
    s = o.scales()
    assert s["sigma_um"] == optics.E1_100X_FROZEN["sigma_um"]
    assert any("FROZEN" in n for n in s["notes"]), "a freeze must not become invisible"


def test_frozen_block_is_within_the_documented_range_of_derivation():
    """0.78x-1.35x. If this moves, the README's number is wrong."""
    o = Optics()
    derived = optics.scales_for(o.psf(), o.dxy_um, o.dz_um)
    r = [optics.E1_100X_FROZEN[k] / derived[k] for k in optics.SCALE_KEYS]
    assert min(r) == pytest.approx(0.78, abs=0.01)
    assert max(r) == pytest.approx(1.35, abs=0.01)


def test_digest_moves_with_the_psf_and_not_with_anything_else():
    a = Optics()
    assert a.digest() == Optics().digest()
    assert a.digest() != Optics(lateral_fwhm_um=0.25).digest()
    # A different PSF SOURCE string is still a different PSF record, deliberately.
    assert a.digest() != Optics(psf_source="theory").digest()


def test_theoretical_psf_is_the_widefield_limit():
    p = optics.theoretical(525.0, 1.49, 1.515)
    assert p.lateral_fwhm_um == pytest.approx(0.51 * 0.525 / 1.49)
    assert p.axial_fwhm_um == pytest.approx(2 * 1.515 * 0.525 / 1.49 ** 2)
    assert p.source == "theory"
