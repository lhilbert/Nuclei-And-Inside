"""
Per-nucleus intensity quantification.

Every channel is measured inside every nucleus mask, background-corrected,
and reported alongside enough identifying columns to trace any row back to
the field, the file and the exported crop it came from.
"""

import os

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.measure import regionprops
from skimage.morphology import convex_hull_image


def background_level(vol, labels, dilate_iter=6):
    """
    Median intensity outside all nuclei.

    The nuclear mask is DILATED before the background is taken. Without
    that, the halo of out-of-focus light hugging each nucleus counts as
    background, biases the estimate upward, and depresses every
    background-corrected value in the field.
    """
    fg = ndi.binary_dilation(labels > 0, iterations=dilate_iter)
    bg_pixels = vol[~fg]
    if bg_pixels.size < 100:              # densely packed field -- fall back
        bg_pixels = vol[labels == 0]
    return float(np.median(bg_pixels)) if bg_pixels.size else 0.0


def nucleus_uid(stack, lab):
    """
    Stable identifier for one nucleus: '<file stem>#p<pos>_nuc<label>'.

    Unique across an entire run, and reused as the stem of the exported
    OME-TIFF, so a row in the measurement table and a file on disk carry
    the same name.
    """
    stem = os.path.splitext(os.path.basename(stack.source_path))[0]
    return f"{stem}#p{stack.position:02d}_nuc{int(lab):03d}"


def weighted_radius(coords_um, weights):
    """
    Intensity-weighted mean radius about an object's geometric centre.

    `coords_um` is (N, D) voxel coordinates already scaled to microns;
    `weights` the matching intensities. Returns (r_weighted, r_normalised):

        r_weighted    sum(I * r) / sum(I), in microns
        r_normalised  r_weighted / mean(r), dimensionless

    The normalised form is the one to compare across nuclei: r_weighted
    alone scales with nucleus size, so a big nucleus and a small one with
    identical intensity architecture score differently. Against the
    unweighted mean radius of the same mask, 1.0 means intensity is spread
    exactly like the mask itself, >1 means it is pushed toward the rim,
    <1 means it is concentrated centrally.

    Weights are clipped at zero. After background subtraction a few dim
    voxels go negative, and a negative weight would pull the centroid of
    the intensity distribution in the wrong direction.
    """
    centre = coords_um.mean(axis=0)
    r = np.linalg.norm(coords_um - centre, axis=1)
    r_unweighted = float(r.mean())

    w = np.clip(np.asarray(weights, dtype=np.float64), 0.0, None)
    total = float(w.sum())
    if total <= 0 or r_unweighted <= 0:
        return float("nan"), float("nan")

    r_w = float((w * r).sum() / total)
    return r_w, r_w / r_unweighted


def midplane_metrics(mask3d, dna_crop, voxel_um, bg=0.0):
    """
    Shape and texture on the single z-plane where the nucleus is widest.

    A 2D read-out taken at the nucleus's own mid-section, rather than an
    average over the stack. In a thin slab the top and bottom planes are
    grazing cuts through the nuclear cap: they are small, ragged, and
    dominated by partial-volume effects, so any solidity or CV averaged
    over z is a mixture of real chromatin structure and how much of the
    nucleus the slab happened to catch. The widest plane is the closest
    thing to a reproducible equatorial section.

    Parameters
    ----------
    mask3d : (z, y, x) bool
        Nucleus mask, cropped to its bounding box.
    dna_crop : (z, y, x) float
        DNA channel over the same crop, raw (not background-subtracted).
    voxel_um : (dz, dy, dx)
    bg : float
        Field background for the DNA channel.

    Returns
    -------
    dict of columns, all prefixed `mid_`.
    """
    dz, dy, dx = voxel_um
    zi = int(np.argmax(mask3d.sum(axis=(1, 2))))     # widest plane
    m2 = mask3d[zi]
    vals = np.asarray(dna_crop[zi][m2], dtype=np.float64)

    hull_area = float(convex_hull_image(m2).sum())
    area_px = float(m2.sum())

    # CV is reported both raw and background-corrected because the two
    # answer different questions. The camera offset inflates the mean but
    # not the standard deviation, so the raw CV is systematically damped
    # and depends on detector settings; the corrected one estimates the
    # contrast of the chromatin itself and is what to compare between
    # experiments. Neither is more "true" -- they differ by a known factor.
    mean_raw = float(vals.mean())
    sd = float(vals.std())
    mean_corr = mean_raw - float(bg)

    coords_um = np.argwhere(m2) * np.array([dy, dx], dtype=np.float64)
    r_w, r_norm = weighted_radius(coords_um, vals - bg)

    return dict(
        mid_z_offset=zi,
        # no area column here: it would equal `max_area_um2` from
        # segment.py exactly, and two names for one measurement is how a
        # table grows a feature that then gets five votes in a PCA

        mid_solidity=(area_px / hull_area) if hull_area > 0 else np.nan,
        mid_dna_cv=(sd / mean_raw) if mean_raw > 0 else np.nan,
        mid_dna_cv_corr=(sd / mean_corr) if mean_corr > 0 else np.nan,
        mid_dna_radial_um=r_w,
        mid_dna_radial_norm=r_norm,
    )


def quantify_nuclei(stack, labels, dna_channel, bg_dilate_iter=6,
                    extra_columns=None):
    """
    Measure every channel inside every nucleus.

    Parameters
    ----------
    stack : io.Stack
    labels : (Z, Y, X) int label image from segment.segment_nuclei
    dna_channel : str or int
        Which channel was segmented; texture metrics use this one.
    extra_columns : dict or None
        Constant columns to attach to every row (e.g. condition, replicate).

    Returns
    -------
    DataFrame, one row per nucleus.

    Identifying columns
        nucleus_uid, file, position, field, source_path, label

    Per channel <C>
        <C>_median, <C>_mean, <C>_p90, <C>_std     raw
        <C>_bg                                     field background
        <C>_median_corr, <C>_mean_corr             background-subtracted
        <C>_integrated_corr                        sum of (I - bg) over the
                                                   nucleus, in intensity*um^3.
                                                   USE THIS when comparing
                                                   nuclei of different size --
                                                   it is the quantity
                                                   proportional to total
                                                   molecule number.
        <C>_enrichment                             median / background

    Chromatin texture on the DNA channel
        dna_cv, dna_p90_p50, dna_p90_p10, dna_top10_fraction

    Radial distribution of DNA intensity, whole 3D mask
        dna_radial_um        intensity-weighted mean radius from the mask
                             centre, in microns
        dna_radial_norm      the same divided by the unweighted mean radius.
                             USE THIS to compare nuclei: 1.0 = intensity
                             spread like the mask, >1 = rim-biased,
                             <1 = centrally concentrated.

    Mid-plane read-out (widest z-plane of each nucleus; see midplane_metrics)
        mid_z                absolute z index of that plane in the stack
                             (its cross-sectional area is `max_area_um2`,
                             from segment.py -- same measurement)
        mid_solidity         2D solidity, area / convex-hull area
        mid_dna_cv           CV of DNA intensity, raw
        mid_dna_cv_corr      CV after background subtraction
        mid_dna_radial_um    intensity-weighted radius within that plane
        mid_dna_radial_norm  normalised, as dna_radial_norm above
    """
    voxel_volume = stack.voxel_volume_um3
    dz, dy, dx = stack.voxel_um

    bg = {c: background_level(stack.channel(c).astype(np.float32), labels,
                              bg_dilate_iter)
          for c in stack.channels}

    dna_idx = (stack.channels.index(dna_channel)
               if isinstance(dna_channel, str) else dna_channel)
    dna_name = stack.channels[dna_idx]

    rows = []
    for r in regionprops(labels):
        mask, sl = r.image, r.slice

        row = dict(
            nucleus_uid=nucleus_uid(stack, r.label),
            file=os.path.basename(stack.source_path),
            position=stack.position,
            field=stack.name,
            source_path=stack.source_path,
            label=int(r.label),
            dna_channel=dna_name,
            volume_um3=r.area * voxel_volume,
        )
        if extra_columns:
            row.update(extra_columns)

        zc, yc, xc = r.centroid
        row.update(centroid_z_um=zc * dz, centroid_y_um=yc * dy,
                   centroid_x_um=xc * dx)

        for ci, cname in enumerate(stack.channels):
            vals = stack.data[:, ci][sl][mask].astype(np.float32)
            b = bg[cname]
            corrected = vals - b
            row.update({
                f"{cname}_median": float(np.median(vals)),
                f"{cname}_mean": float(vals.mean()),
                f"{cname}_p90": float(np.percentile(vals, 90)),
                f"{cname}_std": float(vals.std()),
                f"{cname}_bg": b,
                f"{cname}_median_corr": float(np.median(vals) - b),
                f"{cname}_mean_corr": float(corrected.mean()),
                f"{cname}_integrated_corr": float(corrected.sum() * voxel_volume),
                f"{cname}_enrichment": float(np.median(vals) / b) if b > 0 else np.nan,
            })

        # chromatin texture, DNA channel only
        dna_crop = stack.data[:, dna_idx][sl].astype(np.float32)
        dna = dna_crop[mask]
        p10, p50, p90 = np.percentile(dna, [10, 50, 90])
        row.update(
            dna_cv=float(dna.std() / dna.mean()),
            dna_p90_p50=float(p90 / p50),
            dna_p90_p10=float(p90 / p10) if p10 > 0 else np.nan,
            dna_top10_fraction=float(dna[dna >= p90].sum() / dna.sum()),
        )

        # intensity-weighted radial distribution, whole 3D mask
        coords3d_um = np.argwhere(mask) * np.array([dz, dy, dx])
        r_w3, r_n3 = weighted_radius(coords3d_um, dna - bg[dna_name])
        row.update(dna_radial_um=r_w3, dna_radial_norm=r_n3)

        # 2D read-out on the widest plane of this nucleus
        mid = midplane_metrics(mask, dna_crop, stack.voxel_um, bg[dna_name])
        mid["mid_z"] = int(sl[0].start + mid.pop("mid_z_offset"))
        # The widest plane is only an equatorial section if the nucleus is
        # fully inside the slab. When it lands on the first or last plane of
        # the stack the nucleus is cut off in z and its true widest section
        # was never imaged, so mid_solidity and mid_dna_cv describe a
        # grazing cut. Filter on this before comparing conditions.
        mid["mid_at_z_border"] = bool(mid["mid_z"] in (0, stack.data.shape[0] - 1))
        row.update(mid)

        rows.append(row)

    return pd.DataFrame(rows)
