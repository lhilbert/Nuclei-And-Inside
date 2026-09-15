"""Response -> binary -> centreline. Thresholding, skeletonisation, pruning, punctum rejection.

Everything here is in microns. A skeleton is an index lattice and lengths taken off it are the
easiest place in the whole pipeline to produce a number in voxels wearing microns, so every
length below sums TRUE EUCLIDEAN STEPS - never voxel count times spacing.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import convolve, label as ndlabel, zoom
from skan import Skeleton, summarize
from skimage.filters import apply_hysteresis_threshold
from skimage.morphology import skeletonize

from .enhance import centerline_nms
from .grid import Grid

# skan branch types
END_END, JUNCTION_END, JUNCTION_JUNCTION, CYCLE = 0, 1, 2, 3


def _sp(spacing_um):
    """A scalar or a 3-tuple, always returned as a 3-tuple."""
    return (spacing_um,) * 3 if np.isscalar(spacing_um) else tuple(spacing_um)


def _degree_map(sk: np.ndarray) -> np.ndarray:
    k = np.ones((3, 3, 3), np.uint8)
    k[1, 1, 1] = 0
    return convolve(sk.astype(np.uint8), k, mode="constant") * sk


def smooth_polyline(pts: np.ndarray, passes: int) -> np.ndarray:
    """Laplacian smoothing of a digitised centreline, endpoints fixed.

    A centreline quantised onto a voxel grid zigzags between the 26 available step directions,
    and summing those steps OVERestimates arc length - measured 10.6% on a straight oblique
    phantom filament. A few passes remove the quantisation without touching real curvature: a
    straight line converges to straight, and an arc keeps its length.
    """
    if passes <= 0 or len(pts) < 3:
        return pts
    p = pts.astype(float).copy()
    for _ in range(passes):
        q = p.copy()
        q[1:-1] = 0.25 * p[:-2] + 0.5 * p[1:-1] + 0.25 * p[2:]
        p = q
    return p


def drop_isolated(sk: np.ndarray) -> np.ndarray:
    """Remove skeleton voxels with no 26-neighbour.

    skan builds a path graph and a lone voxel contributes no path, which makes the CSR index
    pointer empty and raises from deep inside scipy. A lone voxel is not a filament under any
    definition here, so it is dropped before skan sees it.
    """
    return sk & (_degree_map(sk) > 0)


def safe_skeleton(sk: np.ndarray, spacing):
    """`skan.Skeleton` or None. Returns None when there is no path to build."""
    if not sk.any():
        return None
    try:
        return Skeleton(sk, spacing=_sp(spacing))
    except (ValueError, IndexError):
        return None


def remove_small(mask: np.ndarray, min_voxels: int) -> np.ndarray:
    """Drop connected components with fewer than `min_voxels` voxels, **26-connected**.

    Written out rather than calling `skimage.morphology.remove_small_objects` for two reasons.
    Its size parameter was renamed AND its boundary changed in 0.26, so the same call means
    different things in different versions. And its default connectivity is 1, which in 3D is
    6-neighbour: a filament running diagonally is 26-connected but not 6-connected, so the
    default would chop every oblique filament into pieces and then delete them.
    """
    lab, n = ndlabel(mask, np.ones((3, 3, 3)))
    if n == 0:
        return mask
    sizes = np.bincount(lab.ravel())
    keep = sizes >= int(min_voxels)
    keep[0] = False
    return keep[lab]


def binarize(response, null_level, mask, params, vol=None, sigma_vox=None, spacing=None):
    """Hysteresis in units of the NOISE NULL, then across-ridge thinning.

    `lo` and `hi` are multiples of the response of pure noise at the same parameters. They are
    never a percentile of this image: a percentile forces the same foreground fraction onto
    every nucleus, which would make a probe-absent control pass the acceptance test by
    construction rather than by detecting nothing.
    """
    d = params.detect
    lo = float(d.low_k) * float(null_level)
    hi = float(d.high_k) * float(null_level)
    fg = apply_hysteresis_threshold(response, lo, hi) & mask
    fg = remove_small(fg, int(d.min_object_voxels))
    used_nms = False
    if d.thin_by_nms and vol is not None:
        # Thins the PSF's 4.35:1 axial ribbon by taking the RESPONSE maximum instead of the
        # shape's medial axis. Measured: 14x fewer spurious breaks than plain skeletonisation,
        # for 2.4x the runtime and 0.017 um of distance error - a tenth of the resolution.
        fg = fg & centerline_nms(response, vol, sigma_vox, spacing,
                                 step_um=float(params.scales()["nms_step_um"]),
                                 min_response=lo, where=fg)
        used_nms = True
    return fg, {"low": lo, "high": hi, "null": float(null_level), "nms": used_nms}


def resample_to(fg: np.ndarray, src: Grid, dst_spacing_um) -> tuple[np.ndarray, Grid]:
    """Binary from the native grid onto the RESOLUTION-NORMALISED work grid.

    `skeletonize` (Lee) has no spacing parameter: it treats every voxel as a unit cube and
    decides topology on the index grid. The condition it actually needs is not that the grid be
    isotropic in MICRONS, but that the STRUCTURE be isotropic in VOXELS.

    Measured on a phantom: after enhancement a filament's detected cross-section is 1.000 um in
    z by 0.230 um in x - 4.35:1, against a PSF anisotropy of 4:1. On a 0.092 um isotropic grid that
    is an 8 x 2 voxel RIBBON, and skeletonising a ribbon gives a centreline that wanders inside
    it: 7.86 um of skeleton with a spurious branch point for a 6.55 um filament, 20% too long,
    which then loses 3 um to spur pruning.

    So the work grid is Nyquist-matched to the PSF ON EACH AXIS. In those units the blurred
    filament is round, `skeletonize` sees a tube, and no axial information is lost because the
    axial sampling is still at Nyquist. Lengths are then measured with the true anisotropic
    spacing, so geometry stays in microns while topology is decided where the data is isotropic.
    """
    sp = _sp(dst_spacing_um)
    factors = (src.dz_um / sp[0], src.dy_um / sp[1], src.dx_um / sp[2])
    # EXACT BLOCK SCATTER, NOT AN INTERPOLATION AND NOT A FILTER.
    #
    # scipy's zoom interpolates at the new sample positions and does not prefilter, so a
    # structure 2-3 voxels wide sampled at half rate can fall below any threshold everywhere and
    # vanish outright - a whole phantom filament did exactly that. A Gaussian prefilter then
    # erases a one-voxel centreline instead (its peak falls to ~0.15), and a maximum filter
    # over-covers: the smallest integer max filter spans 0.4 um in z, which fused the four arms
    # of an X phantom into a degree-3 junction.
    #
    # So each source voxel is scattered into exactly the target voxel that contains it, OR-ed
    # with a nearest-neighbour resample which handles any axis being upsampled. Nothing is lost
    # and nothing is thickened.
    out_shape = tuple(max(1, int(round(n * f))) for n, f in zip(fg.shape, factors))
    out = zoom(fg.astype(np.uint8), factors, order=0, mode="nearest").astype(bool)
    if out.shape != out_shape:
        out = np.zeros(out_shape, bool)
    zz, yy, xx = np.nonzero(fg)
    if len(zz):
        tz = np.clip((zz * factors[0]).astype(np.int64), 0, out_shape[0] - 1)
        ty = np.clip((yy * factors[1]).astype(np.int64), 0, out_shape[1] - 1)
        tx = np.clip((xx * factors[2]).astype(np.int64), 0, out_shape[2] - 1)
        out[tz, ty, tx] = True
    return out, Grid(sp[0], sp[1], sp[2], *src.origin)


def prune_spurs(sk, spacing_um, min_branch_um: float, max_iter: int = 12):
    """Remove junction-to-tip branches shorter than `min_branch_um`, repeatedly.

    Iterative because removing one spur exposes another: a junction that loses two of its three
    branches becomes an ordinary path point, and the branch beyond it is then a spur.
    """
    sk = sk.copy()
    removed = 0
    for _ in range(max_iter):
        if not sk.any():
            break
        S = safe_skeleton(sk, spacing_um)
        if S is None:
            break
        df = summarize(S, separator="_")
        short = df[(df.branch_type == JUNCTION_END) & (df.branch_distance < min_branch_um)]
        if short.empty:
            break
        deg = _degree_map(sk)
        drop = np.zeros_like(sk)
        for pid in short.index:
            for z, y, x in S.path_coordinates(pid).astype(int):
                if deg[z, y, x] <= 2:            # keep the junction itself
                    drop[z, y, x] = True
        if not drop.any():
            break
        sk &= ~drop
        removed += int(drop.sum())
    return sk, removed


def path_length_um(sk: np.ndarray, spacing_um) -> float:
    """Total centreline length, summing real Euclidean steps between skeleton voxels.

    NOT voxel count times spacing. A digitised diagonal steps by sqrt(3) * spacing per voxel, so
    counting voxels underestimates an oblique filament by up to 1.73x - a length in voxels
    wearing microns, which is the error class this package exists to avoid.

    **This is a QC number, not the measurement.** Summing raw voxel steps goes the OTHER way and
    OVERestimates, because a digitised line zigzags between the 26 available step directions:
    measured 1.33x on a straight oblique phantom filament. Use `polyline_length_um`, which is
    what the graph's edge lengths are, and what `trace()` reports as `length_um`.
    """
    S = safe_skeleton(sk, spacing_um)
    if S is None:
        return 0.0
    return float(summarize(S, separator="_").branch_distance.sum())


def polyline_length_um(sk: np.ndarray, spacing_um, passes: int) -> float:
    """Total centreline length over SMOOTHED polylines. **This is the measurement.**

    The quantisation zigzag that `path_length_um` sums is an artefact of the lattice, not
    curvature of the object. Measured on a straight oblique phantom filament of known length:
    raw steps give 1.33x truth, and 1, 3, 5, 8 smoothing passes give 1.11x, 1.03x, 1.01x, 1.00x.
    Three passes is the shipped value and leaves a +3% bias, which is recorded rather than
    tuned away - more passes would also start shortening real curvature.

    Needs no image, only the skeleton, so the acceptance test reports the same length whether
    or not the flattened actin was kept.
    """
    S = safe_skeleton(sk, spacing_um)
    if S is None:
        return 0.0
    total = 0.0
    sp = np.array(_sp(spacing_um))
    for pid in range(S.n_paths):
        pts = smooth_polyline(S.path_coordinates(pid).astype(float) * sp, passes)
        total += float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
    return total


def component_lengths(sk: np.ndarray, spacing_um, smooth_passes: int = 0):
    """(component label image, {label: centreline length in um}).

    `smooth_passes > 0` measures over SMOOTHED polylines, which is the physical length. Pass it
    whenever the result is compared against a threshold stated in microns: the raw voxel-step
    sum is ~1.3x too long on an oblique path, so `max_loop_perimeter_um = 2.0` compared against
    it is really a 1.5 um perimeter, and `min_branch_um` is really a shorter minimum than it
    says. Both thresholds are named for a physical length and must be given one.
    """
    lab, n = ndlabel(sk, np.ones((3, 3, 3)))
    out = {i: 0.0 for i in range(1, n + 1)}
    if n == 0 or not sk.any():
        return lab, out
    S = safe_skeleton(sk, spacing_um)
    if S is None:
        return lab, out
    df = summarize(S, separator="_")
    sp = np.array(_sp(spacing_um))
    for pid in range(S.n_paths):
        c = S.path_coordinates(pid).astype(float)
        if smooth_passes > 0:
            pts = smooth_polyline(c * sp, smooth_passes)
            L = float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
        else:
            L = float(df.branch_distance.iloc[pid])
        i0 = c.astype(int)[0]
        out[int(lab[i0[0], i0[1], i0[2]])] += L
    return lab, out


def reject_puncta(sk, spacing_um, max_loop_perimeter_um: float, max_loop_area_um2: float,
                  smooth_passes: int = 0):
    """Drop small closed loops.

    Measured: a probe-absent control's false positives are RINGS. A bright, flat-topped punctum
    has a genuine ridge around its rim, so a ridge filter traces a small closed curve. A
    filament does not close on itself at that scale. The test is topological - no endpoints -
    plus a size bound, never an intensity threshold.

    **This filter is too narrow, and that is measured, not suspected.** On a real unstained
    control it rejected 22 rims out of well over a hundred, because a rim that is broken, or
    fused to a neighbouring rim, HAS free ends and so is not a closed loop at all. The count is
    returned rather than swallowed precisely so that this stays visible.
    """
    lab, lens = component_lengths(sk, spacing_um, smooth_passes)
    deg = _degree_map(sk)
    drop = np.zeros_like(sk)
    n_rej = 0
    equiv_d = 2.0 * np.sqrt(max_loop_area_um2 / np.pi)
    sp = _sp(spacing_um)
    for i in sorted(lens):
        m = lab == i
        if ((deg == 1) & m).any():
            continue                              # has a free end: not a closed loop
        zz, yy, xx = np.nonzero(m)
        span = float(np.sqrt((np.ptp(zz) * sp[0]) ** 2 + (np.ptp(yy) * sp[1]) ** 2
                             + (np.ptp(xx) * sp[2]) ** 2))
        if lens[i] <= max_loop_perimeter_um and span <= equiv_d * 1.5:
            drop |= m
            n_rej += 1
    return sk & ~drop, n_rej


def drop_short_components(sk, spacing_um, min_um: float, smooth_passes: int = 0):
    """Drop whole components shorter than `min_um` of TRUE centreline length."""
    lab, lens = component_lengths(sk, spacing_um, smooth_passes)
    if not lens:
        return sk, 0
    keep = np.zeros(max(lens) + 1, bool)
    for i, L in lens.items():
        keep[i] = L >= min_um
    return keep[lab] & sk, int(sum(1 for L in lens.values() if L < min_um))


def trace(fg: np.ndarray, src: Grid, params, min_branch_um: float | None = None) -> dict:
    """Binary on the native grid -> pruned skeleton on the work grid.

    `min_branch_um` overrides the parameter block. The acceptance test sweeps it, and it has to
    RE-TRACE at each value rather than filter afterwards: pruning a spur changes the topology,
    so a shorter minimum branch is not the same graph with more edges kept.
    """
    t = params.trace
    s = params.scales()
    work_spacing = (s["work_z_um"], s["work_xy_um"], s["work_xy_um"])
    min_branch = float(t.min_branch_um if min_branch_um is None else min_branch_um)

    work_fg, wgrid = resample_to(fg, src, work_spacing)
    sp = wgrid.spacing
    sk = drop_isolated(skeletonize(work_fg))
    n0 = int(sk.sum())
    sk, n_pruned = prune_spurs(sk, sp, min_branch)
    n_rej = 0
    if t.reject_puncta:
        sk, n_rej = reject_puncta(sk, sp, t.max_loop_perimeter_um, t.max_loop_area_um2,
                                  t.polyline_smooth_passes)
    sk, n_short = drop_short_components(sk, sp, min_branch, t.polyline_smooth_passes)
    sk = drop_isolated(sk)
    return {"skeleton": sk, "grid": wgrid, "fg_work": work_fg,
            "min_branch_um": min_branch,
            "n_skel_voxels_raw": n0, "n_pruned_voxels": n_pruned,
            "n_puncta_rejected": n_rej, "n_short_components_dropped": n_short,
            # THE measurement: smoothed polylines, the same lengths the graph's edges carry.
            "length_um": polyline_length_um(sk, sp, t.polyline_smooth_passes),
            # The raw voxel-step sum. A QC number - it is ~1.33x too long on an oblique
            # filament because a digitised line zigzags. Reported so the gap stays visible.
            "skeleton_length_um": path_length_um(sk, sp)}
