"""QC figures. Look at them before using the tables.

    qc/<stem>_pNN_overview.png   per channel: z max projection, nucleus
                                 outlines and label ids, every detected spot
    qc/nuclei/<nucleus_uid>.png  per nucleus, per channel: raw MIP with the
                                 (eroded) mask outline, DoG MIP, spot outlines
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from scipy import ndimage as ndi         # noqa: E402
from skimage.segmentation import find_boundaries   # noqa: E402

COLOURS = ("#ff4fd8", "#3cff3c", "#3fd2ff", "#ffd23f")


def _stretch(img, lo=1.0, hi=99.9):
    a, b = np.percentile(img, [lo, hi])
    return np.clip((img.astype(np.float32) - a) / max(b - a, 1e-6), 0, 1)


def _outline_rgba(mask2d, rgb=(1.0, 1.0, 0.0)):
    edge = find_boundaries(mask2d, mode="inner")
    out = np.zeros(mask2d.shape + (4,), np.float32)
    out[edge] = (*rgb, 1.0)
    return out


def field_overview(images: dict, labels, spots, params, voxel_um, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lab_mip = labels.max(axis=0)
    objs = ndi.find_objects(lab_mip)
    n = len(params.channels)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 7.4), squeeze=False)
    for ax, s, col in zip(axes[0], params.channels, COLOURS):
        ax.imshow(_stretch(images[s.name].max(axis=0)), cmap="gray", interpolation="nearest")
        ax.imshow(_outline_rgba(lab_mip, (1.0, 0.85, 0.0)), interpolation="nearest")
        for lid, sl in enumerate(objs, start=1):
            if sl is not None:
                ax.text(sl[1].start, sl[0].start, str(lid), color="yellow", fontsize=6)
        d = spots[spots["channel"] == s.name] if len(spots) else spots
        if len(d):
            ax.scatter(d["centroid_x_px"], d["centroid_y_px"], s=40, facecolors="none",
                       edgecolors=col, linewidths=0.8)
        ax.set_title(f"{s.name} (nd2 idx {s.index}) - {len(d)} spots, z-MIP", fontsize=10)
        ax.set_axis_off()
    fig.suptitle(f"{path.stem}   erosion {params.erosion_um} um", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def nucleus_panels(crops, params, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    by_nuc = defaultdict(dict)
    for c in crops:
        by_nuc[c["nucleus_uid"]][c["channel"]] = c
    for uid, chans in by_nuc.items():
        names = [s.name for s in params.channels if s.name in chans]
        fig, axes = plt.subplots(len(names), 3, figsize=(10.5, 3.5 * len(names)),
                                 squeeze=False)
        for row, name in zip(axes, names):
            c = chans[name]
            raw = _stretch(c["raw"], 1, 99.95)
            row[0].imshow(raw, cmap="gray", interpolation="nearest")
            row[0].imshow(_outline_rgba(c["mask"]), interpolation="nearest")
            row[0].set_title(f"{name} raw MIP + mask", fontsize=9)
            row[1].imshow(c["dog"], cmap="magma", interpolation="nearest",
                          vmin=0, vmax=max(c["threshold"] * 2, 1e-6))
            row[1].set_title(f"DoG MIP (0..2x thr={c['threshold']:.1f})", fontsize=9)
            row[2].imshow(raw, cmap="gray", interpolation="nearest")
            row[2].imshow(_outline_rgba(c["spots"] > 0, (0.2, 1.0, 1.0)),
                          interpolation="nearest")
            oz, oy, ox = c["origin"]
            for f in c["feats"]:
                row[2].text(f["centroid_x_px"] - ox + 6, f["centroid_y_px"] - oy - 6,
                            f"{f['spot_id']}: {f['volume_um3']:.2f}", color="cyan",
                            fontsize=7)
            row[2].set_title(f"{len(c['feats'])} spots (id: volume um3)", fontsize=9)
            for ax in row:
                ax.set_axis_off()
        fig.suptitle(uid, fontsize=10)
        fig.tight_layout()
        fig.savefig(outdir / f"{uid}.png", dpi=110)
        plt.close(fig)
