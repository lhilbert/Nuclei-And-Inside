"""The PSF, and everything derived from it.

TEN numbers in this pipeline are not free parameters - they are the PSF wearing different
units, and every one of them changes with the optics:

    sigma_um              lateral Hessian scale           = lateral FWHM / 2.355 * sigma_scale
    sigma_z_um            axial   Hessian scale           = axial   FWHM / 2.355 * sigma_scale
    work_xy_um            lateral work-grid spacing       = lateral FWHM / 2, clamped to native
    work_z_um             axial   work-grid spacing       = axial   FWHM / 2, clamped to native
    nms_step_um           across-ridge NMS step           = lateral FWHM / 2
    junction_merge_um     one junction, not two           = 3 x work_xy
    tangent_skip_um       arc discarded before a fit      = 1 x work_z
    tangent_window_um     arc the tangent is fitted over  = 3 x work_z
    validate_nms_step_um  the same NMS step, in validation
    validate_tolerance_um a hit against a reference       = lateral FWHM

The last three of the first group are multiples of the WORK GRID rather than of the PSF
directly, because what they are about is the lattice the skeleton lives on. They still move
with the optics, through it.

**`min_branch_um` is NOT in this list and must never join it.** It is a statement about what
counts as an antenna - biology, not optics - so deriving it would make an antenna on one
microscope a different object from an antenna on another. It is declared, swept, and reported
as a curve instead.

The defaults shipped with this package are widefield 100x, NA 1.49, measured at 0.2302 um
lateral and 0.5923 um axial FWHM - 5.0 samples per lateral FWHM at the 0.046 um pixel. They do
not transfer. An image-scanning microscope on the same objective has a different PSF at a
different sampling, and hard-coding these values there would run a lateral Hessian at 1.5x the
right scale and decide topology on a grid that is not Nyquist. `scales_for` raises
`SamplingError` rather than doing that quietly.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

FWHM_PER_SIGMA = 2.3548200450309493


@dataclass(frozen=True)
class PSF:
    lateral_fwhm_um: float
    axial_fwhm_um: float
    source: str = "theory"

    @property
    def lateral_sigma_um(self) -> float:
        return self.lateral_fwhm_um / FWHM_PER_SIGMA

    @property
    def axial_sigma_um(self) -> float:
        return self.axial_fwhm_um / FWHM_PER_SIGMA

    @property
    def anisotropy(self) -> float:
        return self.axial_fwhm_um / self.lateral_fwhm_um

    def nyquist(self) -> tuple[float, float]:
        """(lateral, axial) spacing at which the PSF is critically sampled."""
        return (self.lateral_fwhm_um / 2.0, self.axial_fwhm_um / 2.0)

    def to_dict(self) -> dict:
        return {"lateral_fwhm_um": self.lateral_fwhm_um, "axial_fwhm_um": self.axial_fwhm_um,
                "source": self.source}


def theoretical(emission_nm: float, na: float, immersion_ri: float = 1.515,
                resolution_gain: float = 1.0) -> PSF:
    """Widefield diffraction limit, optionally scaled by a modality's resolution gain.

    `resolution_gain` is 1.0 for widefield. Image-scanning modalities such as iSIM improve on
    it, but by how much is an instrument property and a claim about the data - so it is
    DECLARED and CHECKED against a measurement, never assumed here.
    """
    lam = emission_nm / 1000.0
    lateral = 0.51 * lam / na / resolution_gain
    axial = 2.0 * immersion_ri * lam / (na ** 2) / resolution_gain
    return PSF(lateral, axial, "theory")


#: The ten derived scales, in the order they are reported. Held as a tuple because
#: `run_parameters.json` and `scales_digest` both need a stable order, and because a list of
#: names in one place is how you notice that only four of them were updated for new optics.
SCALE_KEYS = (
    "sigma_um", "sigma_z_um",
    "work_xy_um", "work_z_um",
    "nms_step_um",
    "junction_merge_um",
    "tangent_skip_um", "tangent_window_um",
    "validate_nms_step_um", "validate_tolerance_um",
)

#: Three scales are multiples of the WORK GRID rather than of the PSF, because what they are
#: about is the lattice the skeleton lives on. The skeleton climbs out of a junction over about
#: one work-grid z voxel, and a junction splits across ~3 lateral work voxels.
JUNCTION_MERGE_IN_WORK_XY = 3.0
TANGENT_SKIP_IN_WORK_Z = 1.0
TANGENT_WINDOW_IN_WORK_Z = 3.0

#: Below this, a Gaussian derivative is not a derivative: the kernel is narrower than the
#: sample spacing, so what comes out is a difference of neighbouring pixels dressed up as a
#: measurement at a scale that was never resolved.
MIN_SIGMA_PX = 0.9


class SamplingError(RuntimeError):
    """The sampling cannot support a PSF-matched Hessian. Not a bug - a fact about the data."""


def samples_per_fwhm(psf: PSF, native_dxy_um: float, native_dz_um: float) -> tuple[float, float]:
    return (psf.lateral_fwhm_um / native_dxy_um, psf.axial_fwhm_um / native_dz_um)


def _scales(psf: PSF, work_xy: float, work_z: float, sigma_scale: float) -> dict:
    """The ten derived scales, given an already-resolved work grid. In microns."""
    half_fwhm = psf.lateral_fwhm_um / 2.0
    return {
        "sigma_um": psf.lateral_sigma_um * sigma_scale,
        "sigma_z_um": psf.axial_sigma_um * sigma_scale,
        "work_xy_um": work_xy,
        "work_z_um": work_z,
        "nms_step_um": half_fwhm,
        "junction_merge_um": JUNCTION_MERGE_IN_WORK_XY * work_xy,
        "tangent_skip_um": TANGENT_SKIP_IN_WORK_Z * work_z,
        "tangent_window_um": TANGENT_WINDOW_IN_WORK_Z * work_z,
        "validate_nms_step_um": half_fwhm,
        "validate_tolerance_um": psf.lateral_fwhm_um,
    }


def derived_scales(psf: PSF, sigma_scale: float = 1.0) -> dict:
    """The ten scales at the PSF's own Nyquist grid, with no sampling clamp. The theory.

    `scales_for` is what a run must use, because it also knows what the detector sampled.
    """
    ny_xy, ny_z = psf.nyquist()
    return _scales(psf, ny_xy, ny_z, sigma_scale)


def scales_for(psf: PSF, native_dxy_um: float, native_dz_um: float,
               sigma_scale: float = 1.0, strict: bool = True) -> dict:
    """The ten PSF-derived scales, clamped to what the sampling can actually support.

    Two things are checked:

    * **The Hessian scale must be at least ~1 pixel, ON BOTH AXES.** Widefield 100x at a
      0.046 um pixel measures 2.13 lateral px and is fine. An iSIM stack at 0.065 um measures
      **0.83 lateral and 0.87 axial** - and a filter that aliases in z aliases just as badly as
      one that aliases in x, which is why both are checked and not only the lateral one.
    * **The work grid cannot be finer than the acquisition.** Where Nyquist for the measured
      PSF is finer than the pixel, the grid is clamped to native and topology is decided at
      slightly coarser than Nyquist. That is a fact about the acquisition, not a choice, and it
      is recorded in `notes` rather than hidden.

    The clamp happens BEFORE the work-grid-derived scales are computed, so a clamped grid moves
    `junction_merge_um` and the two tangent windows with it.
    """
    ny_xy, ny_z = psf.nyquist()
    notes = []
    if ny_xy < native_dxy_um:
        notes.append(f"work grid clamped to native lateral sampling: Nyquist {ny_xy:.4f} um is "
                     f"finer than the {native_dxy_um:.4f} um pixel")
        ny_xy = native_dxy_um
    if ny_z < native_dz_um:
        notes.append(f"work grid clamped to native axial sampling: Nyquist {ny_z:.4f} um is "
                     f"finer than the {native_dz_um:.4f} um step")
        ny_z = native_dz_um

    out = _scales(psf, ny_xy, ny_z, sigma_scale)
    sig_px = out["sigma_um"] / native_dxy_um
    sig_z_px = out["sigma_z_um"] / native_dz_um

    under = [(name, sig, px, d, f) for name, sig, px, d, f in (
        ("lateral", out["sigma_um"], sig_px, native_dxy_um, psf.lateral_fwhm_um),
        ("axial", out["sigma_z_um"], sig_z_px, native_dz_um, psf.axial_fwhm_um),
    ) if px < MIN_SIGMA_PX]
    if under:
        parts = [f"{name} Hessian scale is {sig:.4f} um = {px:.3f} native "
                 f"{'pixels' if name == 'lateral' else 'steps'} ({d:.4f} um sampling against a "
                 f"{f:.4f} um FWHM = {f / d:.2f} samples per FWHM)"
                 for name, sig, px, d, f in under]
        msg = ("at sigma_scale=%.4f the %s, below the %s px at which a Gaussian derivative stops "
               "being a derivative. This sampling cannot support single-scale PSF-matched "
               "enhancement." % (sigma_scale, "; and the ".join(parts), MIN_SIGMA_PX))
        if strict:
            raise SamplingError(msg)
        notes.append(msg)

    out.update({"sigma_lateral_px": sig_px, "sigma_axial_px": sig_z_px,
                "samples_per_lateral_fwhm": psf.lateral_fwhm_um / native_dxy_um,
                "samples_per_axial_fwhm": psf.axial_fwhm_um / native_dz_um,
                "notes": notes})
    return out


#: The scales the published widefield-100x results were produced with.
#:
#: They are NOT `scales_for(PSF(0.2302, 0.5923), 0.046, 0.2)`. They were built from the
#: THEORETICAL PSF at NA 1.49 (0.1797 / 0.7165 um) and then hand-rounded - 3 x work_xy = 0.276
#: shipped as 0.30, 1 x work_z = 0.3585 as 0.40, 3 x work_z = 1.0755 as 1.0, and one lateral
#: FWHM 0.1797 as 0.184 (= 4 native pixels). Against the measured PSF they are 0.78x-1.35x what
#: derivation gives.
#:
#: Pass this to `DetectParams.from_scales(...)` to reproduce those numbers exactly. For new
#: data, derive from your own measured PSF instead - a literal block is all-or-nothing, and
#: half of one is how a pipeline ends up running two different optics at once.
E1_100X_FROZEN = {
    "sigma_um": 0.078,
    "sigma_z_um": 0.31,
    "work_xy_um": 0.092,
    "work_z_um": 0.3585,
    "nms_step_um": 0.09,
    "junction_merge_um": 0.30,
    "tangent_skip_um": 0.40,
    "tangent_window_um": 1.0,
    "validate_nms_step_um": 0.09,
    "validate_tolerance_um": 0.184,
}


def scales_digest(psf: PSF, scales: dict) -> str:
    """A digest of the PSF and the ten scales, and of nothing else.

    It answers one question: were these two runs detected at the same physical scales? A hash
    over the whole parameter set would move when a comment moves; this moves only when a number
    that reaches the detector does.
    """
    blob = json.dumps({"psf": psf.to_dict(),
                       "scales": {k: round(float(scales[k]), 12) for k in SCALE_KEYS}},
                      sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()[:12]
