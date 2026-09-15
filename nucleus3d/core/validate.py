"""
Per-field validation figures.

The failure modes of this pipeline are visual, not exceptional: merged
neighbours, over-split nuclei, boundaries cutting inside the chromatin edge,
and missed dim nuclei all produce clean-looking tables and no error. Look at
these figures for every run.

Each figure has five panels:

    1. DNA max projection                 what the data looks like
    2. segmentation outlines + labels     what was found; label IDs match
                                          the `label` column in the table
    3. XZ side view                       shows z-clipping and axial extent
    4. unsegmented foreground             DNA signal not assigned to any
                                          nucleus -- the missed-nucleus check
    5. size vs solidity scatter           outliers here are merges or debris
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label as sklabel, regionprops  # noqa: E402
from skimage.segmentation import find_boundaries       # noqa: E402


def _stretch(img, lo_pct=2, hi_pct=99.7):
    lo, hi = np.percentile(img, [lo_pct, hi_pct])
    return np.clip((img - lo) / max(hi - lo, 1e-9), 0, 1)


def unsegmented_fraction(vol, labels, min_blob_um3=5.0, voxel_um=(1, 1, 1),
                         rim_tolerance_um=0.5):
    """
    Fraction of DNA foreground not covered by any nucleus, plus the list of
    sizeable uncovered blobs.

    Foreground is an Otsu threshold on the smoothed raw image -- deliberately
    permissive. A percentile threshold only finds BRIGHT signal and will
    report "nothing missed" while dim nuclei sit unsegmented. That mistake
    is easy to make and gives a falsely clean result.

    Nuclear masks are dilated by `rim_tolerance_um` before the comparison.
    Without it, the thin halo where the permissive threshold extends a little
    past each boundary counts as "missed" and dominates the statistic,
    burying the signal that actually matters: a whole nucleus with no label.
    """
    sm = gaussian(vol.astype(np.float32), sigma=(2, 4, 4))
    fg = sm > threshold_otsu(sm)

    rim = tuple(max(int(round(rim_tolerance_um / s)), 0) for s in voxel_um)
    covered = (ndi.binary_dilation(labels > 0,
                                   structure=np.ones((3, 3, 3)),
                                   iterations=max(rim[1], 1))
               if rim_tolerance_um > 0 else labels > 0)

    uncovered = fg & ~covered
    frac = float(uncovered.sum()) / max(float(fg.sum()), 1.0)

    vvox = float(np.prod(voxel_um))
    blobs = [r.area * vvox for r in regionprops(sklabel(uncovered))
             if r.area * vvox > min_blob_um3]
    return frac, sorted(blobs, reverse=True), fg, uncovered


def validation_figure(stack, labels, props, dna_channel, outpath,
                      min_blob_um3=5.0):
    """
    Write the five-panel validation figure for one field.

    Returns a dict of QC numbers, which the runner collects into a
    per-field summary table so problem fields can be found without
    opening every image.
    """
    vol = stack.channel(dna_channel)
    base = _stretch(vol.max(axis=0))
    dz, dy, _ = stack.voxel_um

    frac, blobs, fg, uncovered = unsegmented_fraction(
        vol, labels, min_blob_um3, stack.voxel_um)

    # The XZ view is a wide, short strip; giving it a full equal-width column
    # leaves most of that column empty. Three square image panels, then a
    # narrow column holding the XZ strip above the scatter.
    fig = plt.figure(figsize=(19, 5.4))
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.9],
                          height_ratios=[0.30, 1.0], hspace=0.22, wspace=0.10,
                          left=0.010, right=0.975, top=0.90, bottom=0.11)
    ax0 = fig.add_subplot(gs[:, 0])
    ax1 = fig.add_subplot(gs[:, 1])
    ax2 = fig.add_subplot(gs[:, 2])
    ax_xz = fig.add_subplot(gs[0, 3])
    ax_sc = fig.add_subplot(gs[1, 3])

    ax0.imshow(base, cmap="gray")
    ax0.set_title(f"{dna_channel} max projection", fontsize=10)

    ov = np.dstack([base] * 3)
    if labels.max():
        ov[find_boundaries(labels.max(axis=0), mode="outer")] = [1, 0.55, 0]
    ax1.imshow(ov)
    for r in regionprops(labels.max(axis=0)):
        ax1.text(r.centroid[1], r.centroid[0], str(r.label), color="yellow",
                 fontsize=11, ha="center", va="center", fontweight="bold")
    ax1.set_title(f"segmentation — {len(props)} nuclei", fontsize=10)

    unc2d = uncovered.max(axis=0)
    ax2.imshow(base, cmap="gray")
    ax2.imshow(np.ma.masked_where(~unc2d, unc2d), cmap="cool", alpha=0.9)
    ax2.set_title(f"unsegmented DNA — {frac:.1%} of foreground\n"
                  f"{len(blobs)} blob(s) > {min_blob_um3:g} µm³", fontsize=10)

    ax_xz.imshow(_stretch(vol.max(axis=2)), cmap="gray", aspect=dz / dy)
    ax_xz.set_title(f"XZ side view — slab {stack.z_extent_um:.1f} µm",
                    fontsize=10)

    for a in (ax0, ax1, ax2, ax_xz):
        a.axis("off")

    if len(props):
        ax_sc.scatter(props.max_area_um2, props.solidity, s=42,
                      c="#3a7ca5", edgecolor="white", zorder=3)
        for r in props.itertuples():
            ax_sc.annotate(str(r.label), (r.max_area_um2, r.solidity),
                           fontsize=8, xytext=(3, 3), textcoords="offset points")
    ax_sc.set_xlabel("max cross-section (µm²)", fontsize=9)
    ax_sc.set_ylabel("solidity", fontsize=9)
    ax_sc.set_title("size vs shape", fontsize=10)
    ax_sc.tick_params(labelsize=8)
    ax_sc.grid(alpha=0.3)

    fig.suptitle(stack.name, fontsize=12)
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fig.savefig(outpath, dpi=95)
    plt.close(fig)

    return dict(field=stack.name, n_nuclei=len(props),
                unsegmented_fraction=round(frac, 4),
                n_uncovered_blobs=len(blobs),
                largest_uncovered_um3=round(blobs[0], 1) if blobs else 0.0,
                median_area_um2=round(float(props.max_area_um2.median()), 2)
                if len(props) else np.nan,
                min_solidity=round(float(props.solidity.min()), 3)
                if len(props) else np.nan,
                qc_figure=os.path.basename(outpath))
