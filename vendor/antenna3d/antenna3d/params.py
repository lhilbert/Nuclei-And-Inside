"""Every parameter in the pipeline, in one place.

`run_analysis.py` holds parameters and no logic; the modules hold logic and no parameters. This
module is the hinge between them: five frozen dataclasses, each field carrying the measurement
or the decision that set its default.

Two rules govern what is allowed to be a field here.

**The ten PSF-derived scales are not fields.** They live on `Optics` and are computed by
`optics.scales_for`, because they are the PSF wearing different units and moving one without the
other nine is how a pipeline ends up running two different optics at once. `Optics.scales()` is
the only way to get them.

**`TraceParams.min_branch_um` IS a field, and deliberately.** It is a statement about what counts
as an antenna - biology, not optics. Deriving it would make an antenna on one microscope a
different object from an antenna on another. It is declared, and the acceptance test reports a
curve over it rather than a number at one value.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, replace
from typing import Any

from . import optics as _optics
from .optics import PSF, SamplingError  # noqa: F401  (re-exported for callers)


@dataclass(frozen=True)
class Optics:
    """The microscope, declared. Checked against every file before its pixels are read.

    `dxy_um` and `dz_um` are the sampling this configuration is FOR, not a description of
    whatever file turns up. A file whose voxel size disagrees is refused rather than analysed at
    the wrong physical scale - one folder in the validation data holds 100x and 40x stacks
    together, and the only thing that distinguishes them reliably is the voxel size.
    """

    #: Measured PSF. The defaults are widefield 100x NA 1.49, 324 lateral and 134 axial profiles
    #: over two fields. `optics.theoretical()` gives the diffraction limit if you have no
    #: measurement, but a measurement wins whenever you have one.
    lateral_fwhm_um: float = 0.2302
    axial_fwhm_um: float = 0.5923
    psf_source: str = "measured"

    #: The sampling this configuration is for.
    dxy_um: float = 0.04599853515625
    dz_um: float = 0.2

    #: Corroborating evidence only. If set, a file whose stem does not contain it is refused -
    #: but the voxel size above is what actually decides, and the marker only catches a file
    #: misnamed at acquisition.
    stem_marker: str | None = "100x"

    na: float | None = 1.49
    immersion_ri: float = 1.515
    emission_actin_nm: float = 525.0

    #: Multiplies both Hessian scales. 1.0 is PSF-matched. Raising it is the only legitimate
    #: response to a sampling that cannot support a PSF-matched derivative, and it must then be
    #: described as what it is - coarser than the PSF - not as "filament-matched".
    sigma_scale: float = 1.0

    #: Set to `optics.E1_100X_FROZEN` to reproduce the published widefield-100x numbers exactly.
    #: `None` derives all ten from the PSF above, which is what new data should do.
    frozen_scales: dict[str, float] | None = None

    def psf(self) -> PSF:
        return PSF(self.lateral_fwhm_um, self.axial_fwhm_um, self.psf_source)

    def scales(self, strict: bool = True) -> dict:
        """The ten derived scales, plus diagnostics and `notes`.

        With `frozen_scales` set, the frozen values are returned and what derivation WOULD have
        given is reported in `notes` as a ratio - so a freeze stays visible instead of becoming
        invisible. A partial freeze raises: a literal block is all-or-nothing.
        """
        derived = _optics.scales_for(self.psf(), self.dxy_um, self.dz_um,
                                     sigma_scale=self.sigma_scale, strict=strict)
        if self.frozen_scales is None:
            return derived

        missing = [k for k in _optics.SCALE_KEYS if k not in self.frozen_scales]
        if missing:
            raise ValueError(
                f"frozen_scales is partial: {len(missing)} of {len(_optics.SCALE_KEYS)} keys "
                f"missing ({', '.join(missing)}). A literal block is all-or-nothing - half of "
                f"one runs some scales on this PSF and the rest on whatever produced the block.")
        extra = [k for k in self.frozen_scales if k not in _optics.SCALE_KEYS]
        if extra:
            raise ValueError(f"frozen_scales has keys that are not derived scales: {extra}")

        out = dict(derived)
        ratios = []
        for k in _optics.SCALE_KEYS:
            out[k] = float(self.frozen_scales[k])
            ratios.append(out[k] / derived[k])
        out["sigma_lateral_px"] = out["sigma_um"] / self.dxy_um
        out["sigma_axial_px"] = out["sigma_z_um"] / self.dz_um
        out["notes"] = list(derived["notes"]) + [
            f"all {len(_optics.SCALE_KEYS)} scales are FROZEN literals "
            f"({min(ratios):.2f}x-{max(ratios):.2f}x what this PSF would give)"]
        return out

    def digest(self) -> str:
        return _optics.scales_digest(self.psf(), self.scales(strict=False))


#: The optics every shipped default in this package was measured on, with the published
#: numbers' frozen scale block. This is the `OPTICS = ...` a new user starts from and edits.
WIDEFIELD_100X = Optics(frozen_scales=_optics.E1_100X_FROZEN)


@dataclass(frozen=True)
class Channels:
    """Channel roles, by OME NAME. Never by index.

    One dataset in this project stores actin at index 2 and another at index 0 - same
    instrument, five days apart. Resolving by index analyses the DNA channel as actin and never
    errors. Run `describe_file` on every new dataset before you set these.
    """
    dna: str = "DAPI"
    actin: str = "GFP"


@dataclass(frozen=True)
class SegParams:
    """Nucleus segmentation, and the gates on what counts as a nucleus.

    Only used when you let this package segment. The intended path is to pass labels you
    already have - from `nucleus3d`, or any other segmenter - in which case only the grid and
    the gates below apply.
    """
    #: Nuclei are found on a coarse grid. A 10 um object does not need 0.046 um sampling.
    nuclei_grid_um: float = 0.23

    #: Equivalent diameter of the nuclei being looked for.
    diameter_um: float = 9.0
    model: str = "cpsam"
    flow_threshold: float = 0.4
    cellprob_threshold: float = 0.0
    gpu: bool = True

    #: 3D flow field over the volume, with `anisotropy` taken from the grid. The 2D-per-plane
    #: alternative was measured cutting every nucleus in z: median axial edge fraction 0.577
    #: against 0.053 here, and NOT ONE nucleus in any 2D setting tapered away in z. It costs
    #: ~3.4x the runtime and it is worth it.
    do_3d: bool = True

    min_volume_um3: float = 100.0
    max_volume_um3: float = 4000.0
    drop_xy_border_touching: bool = True
    #: Recorded, never dropped: whether a nucleus is cut in z is a covariate, not a reason to
    #: discard it. `axial_edge_fraction` carries it into every table.
    drop_z_truncated: bool = False

    #: marching-cubes step_size for the boundary mesh.
    boundary_surface_step: int = 2


@dataclass(frozen=True)
class DetectParams:
    """Flattening, noise, enhancement and thresholding.

    Every threshold here is a multiple of the response of PURE NOISE at the same parameters,
    which is what makes one frozen value mean the same thing on a bright nucleus and a faint
    one. None of them is a percentile of the image: a percentile forces the same foreground
    fraction on every nucleus and would make a probe-absent control pass by construction.
    """
    #: Context kept around the nuclear bounding box in each crop.
    margin_um: float = 1.5
    #: Much larger than a filament (0.25 um) and much smaller than a nucleus (10 um).
    background_sigma_um: float = 1.5

    #: WHICH noise estimator. Every threshold is referenced to the noise, so this key decides
    #: the meaning of all of them - which is why it is a declared parameter and not a line in a
    #: function body.
    #:   first_difference_x  MAD of the lag-1 difference along x. Cancels SMOOTH structure.
    #:   in_mask_mad         MAD of the values themselves. Measured inflating 10.3x against
    #:                       space-filling smooth structure where the first difference moves
    #:                       1.09x, so the response falls with the structure it should detect.
    noise_estimator: str = "first_difference_x"

    #: Vertical fixed-pattern striping: column-median high-pass sd 27-73 ADU against a
    #: per-pixel noise sd of 412-2122, i.e. 3-6% of noise. Small, coherent, nearly free.
    destripe: bool = True
    destripe_axis: int = 0
    #: Destriped native planes held while cropping one field. One read serves every nucleus
    #: overlapping that plane. NOT the whole field - a 169-plane 2280x2588 field is 4.0 GB.
    plane_cache_mb: int = 2048

    #: Frangi. `alpha` is plate-vs-line, `beta` is blob-vs-line. On a synthetic tube-plus-point
    #: the blobness term suppresses the point relative to the tube by 47x, against Sato's 5x -
    #: which is why this filter and not another one.
    frangi_alpha: float = 0.5
    frangi_beta: float = 0.5
    #: gamma, IN UNITS OF THE MEASURED NOISE SD. With skimage's per-image default a pure-noise
    #: image scores HIGHER than real data (ratio 0.73); fixed in noise units the ratio is 19.7
    #: and is stable from 5 to 20 sd. A frozen threshold on the library default measures the
    #: library, not the specimen.
    gamma_in_noise_sd: float = 5.0
    #: The response percentile of a pure-noise image of the same sd at the same parameters.
    #: Every threshold below is a multiple of THIS.
    null_percentile: float = 99.9

    #: Hysteresis, in units of null_level.
    low_k: float = 1.0
    high_k: float = 3.0
    min_object_voxels: int = 20

    #: Steger-style across-ridge non-maximum suppression before skeletonisation. It costs 2.4x
    #: the runtime and 0.017 um of distance error - a tenth of the resolution - and buys a 14x
    #: reduction in spurious breaks, which is the topology error that corrupts a graph.
    thin_by_nms: bool = True


@dataclass(frozen=True)
class TraceParams:
    """Centreline extraction, junction resolution and gap joining."""

    #: A traced segment shorter than this is not an antenna. 2.8x the lateral PSF FWHM: long
    #: enough that one punctum cannot qualify, short enough to keep the real short segments.
    #: DECLARED, not derived - see the module docstring - and swept by `run_acceptance.py`.
    min_branch_um: float = 0.5

    #: Laplacian passes over each digitised centreline before any length is measured. Summing
    #: raw voxel steps overestimates a straight oblique filament's arc length by 10.6%;
    #: smoothing removes the quantisation zigzag and leaves real curvature alone.
    polyline_smooth_passes: int = 3

    #: The false positives of a probe-absent control are RINGS: the rim of a bright flat-topped
    #: blob is a genuine ridge, so a punctum yields a small closed loop. Rejected by that
    #: signature, and COUNTED per nucleus rather than silently dropped - the count is the
    #: evidence about whether the filter is working on the thing it was designed for.
    reject_puncta: bool = True
    max_loop_perimeter_um: float = 2.0
    max_loop_area_um2: float = 0.35

    #: Explicit crossing resolution, not a side effect of morphological cleanup.
    reconnect: bool = True
    #: Furthest a pair of fragment ends may be joined across. Widening it was measured joining
    #: OTHER objects rather than reuniting one: 83% of joins off-object at 0.6 um, 97% at the
    #: widest tested, with the on-object count flat throughout.
    max_gap_um: float = 0.6
    #: Orientation continuity: the turn between the two end tangents.
    max_angle_deg: float = 35.0
    #: Intensity continuity: the ratio of the two fragments' mean intensity must lie in [1/r, r].
    max_intensity_ratio: float = 2.0


@dataclass(frozen=True)
class GraphParams:
    """What the graph counts as a root, a tip and a cut end."""

    #: 0.0 = the sub-envelope shell is INCLUDED. That shell is where the probe-present /
    #: probe-absent contrast lives: 2.24x, 2.72x, 0.91x at 0.2-1.2 um depth against 0.80x,
    #: 0.87x, 0.20x deeper than 1.6 um.
    interior_erosion_um: float = 0.0
    #: A node within this of the nuclear surface is a root.
    root_tolerance_um: float = 0.5

    #: An end against a MANUFACTURED cut face is neither a root nor a tip. A face counts as cut
    #: where the mask's cross-section there is still above this fraction of its peak - a mask
    #: that tapers away was not cut, one that stops at 80% of its widest was.
    cut_area_fraction: float = 0.30

    #: Measure the apparent cross-section FWHM per segment. It is PSF-limited on this data
    #: (0.251 um measured against a 0.180 um PSF), so both the raw and the quadrature-
    #: deconvolved value are emitted and `width_is_psf_limited` is set. Read it as an upper
    #: bound. Turning it off is the single biggest saving in graph building.
    measure_width: bool = True


@dataclass(frozen=True)
class Params:
    """Everything, bundled. This is what the pipeline takes and what is written to disk."""
    optics: Optics = field(default_factory=lambda: WIDEFIELD_100X)
    channels: Channels = field(default_factory=Channels)
    segment: SegParams = field(default_factory=SegParams)
    detect: DetectParams = field(default_factory=DetectParams)
    trace: TraceParams = field(default_factory=TraceParams)
    graph: GraphParams = field(default_factory=GraphParams)

    def scales(self, strict: bool = True) -> dict:
        return self.optics.scales(strict=strict)

    def to_dict(self) -> dict[str, Any]:
        """Every parameter used, plus the resolved scales. Written to `run_parameters.json`."""
        d = asdict(self)
        s = self.scales(strict=False)
        d["resolved_scales"] = {k: s[k] for k in _optics.SCALE_KEYS}
        d["resolved_scales_notes"] = s["notes"]
        d["sampling"] = {"sigma_lateral_px": s["sigma_lateral_px"],
                         "sigma_axial_px": s["sigma_axial_px"],
                         "samples_per_lateral_fwhm": s["samples_per_lateral_fwhm"],
                         "samples_per_axial_fwhm": s["samples_per_axial_fwhm"]}
        d["scales_digest"] = self.optics.digest()
        return d

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)

    def with_(self, **kw) -> "Params":
        """A copy with some sub-blocks replaced. `p.with_(trace=replace(p.trace, ...))`."""
        return replace(self, **kw)
