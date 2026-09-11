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
            n_voxels=int(r.area),
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
        dna = stack.data[:, dna_idx][sl][mask].astype(np.float32)
        p10, p50, p90 = np.percentile(dna, [10, 50, 90])
        row.update(
            dna_cv=float(dna.std() / dna.mean()),
            dna_p90_p50=float(p90 / p50),
            dna_p90_p10=float(p90 / p10) if p10 > 0 else np.nan,
            dna_top10_fraction=float(dna[dna >= p90].sum() / dna.sum()),
        )
        rows.append(row)

    return pd.DataFrame(rows)
