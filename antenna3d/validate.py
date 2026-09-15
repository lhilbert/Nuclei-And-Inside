"""QC figures and the per-field QC numbers. The review step is a human looking at these.

There are two levels, and they exist for different failures.

The **per-field figure** is for "did this field work at all" - segmentation, threshold placement,
slab geometry. The **per-nucleus trace overlay** is for "is the thing it traced a filament", and
it is not optional decoration: on a probe-absent control it is what showed the tracer following
rings on the rims of autofluorescent puncta, which no summary statistic revealed. Panel 2 and
the trace overlays are the two that catch a detector that is confidently wrong.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.ndimage import binary_erosion, gaussian_filter, zoom


def norm(a, lo=1.0, hi=99.7):
    a = np.asarray(a, np.float32)
    v0, v1 = np.percentile(a, lo), np.percentile(a, hi)
    return np.clip((a - v0) / max(v1 - v0, 1e-9), 0, 1)


def norm_sd(a, noise_sd, lo_sd=-1.0, hi_sd=8.0):
    """Scale a flattened image in units of ITS OWN NOISE, not by percentiles.

    A percentile stretch is wrong for this data and it once hid a result: the flattened actin
    has a long bright tail from out-of-nucleus puncta, so the 99.7th percentile sits far above
    the filament level and renders a perfectly good network almost black. Reviewing overlays on
    that stretch would have led to rejecting a detector that was working.
    """
    a = np.asarray(a, np.float32)
    v0, v1 = lo_sd * noise_sd, hi_sd * noise_sd
    return np.clip((a - v0) / max(v1 - v0, 1e-9), 0, 1)


def to_native(sk_work, native_shape):
    """A work-grid boolean back onto the native crop grid, for overlaying on the image."""
    f = tuple(n / s for n, s in zip(native_shape, sk_work.shape))
    return zoom(sk_work.astype(np.float32), f, order=0, mode="nearest") > 0.5


def trace_overlay(crop, tr, out_path, slab_um: float = 1.0, dpi: int = 110):
    """Thin-MIP of the flattened actin, the trace in yellow, the envelope in red.

    **Look at these.** A length density cannot tell you that what was traced is the rim of a
    blob; this can.
    """
    from matplotlib import pyplot as plt
    actin, mask = crop.actin_flat, crop.mask
    skel = to_native(tr["skeleton"], actin.shape)
    Z = actin.shape[0]
    cnt = mask.reshape(Z, -1).sum(1)
    mid = int(np.argmax(cnt)) if cnt.any() else Z // 2
    half = max(1, int(round(slab_um / crop.grid.dz_um / 2)))
    lo, hi = max(0, mid - half), min(Z, mid + half + 1)
    img = np.stack([gaussian_filter(actin[i], 1.0) for i in range(lo, hi)]).max(0)
    g = norm_sd(img, crop.noise_sd) if crop.noise_sd else norm(img)
    rgb = np.dstack([g, g, g])
    rgb[skel[lo:hi].any(0)] = [1.0, 0.95, 0.1]
    mm = mask[mid]
    if mm.any():
        rgb[mm & ~binary_erosion(mm, np.ones((3, 3)))] = [1.0, 0.25, 0.25]
    fig, ax = plt.subplots(figsize=(6, 6 * rgb.shape[0] / max(rgb.shape[1], 1)))
    ax.imshow(rgb, interpolation="nearest")
    ax.axis("off")
    ax.set_title(f"{crop.nucleus_uid}\n{tr['length_um']:.1f} um traced, "
                 f"{tr['n_puncta_rejected']} puncta rejected", fontsize=7)
    fig.tight_layout(pad=0.2)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def field_figure(field, labels, lgrid, graphs, params, out_path, dpi: int = 110):
    """Five panels per field. Read them in order; each answers a different question."""
    from matplotlib import pyplot as plt

    # Plane by plane rather than reading the whole volume: a field here is 4.0 GB.
    mip = None
    for z in range(field.shape[0]):
        pl = field.plane(z, "actin").astype(np.float32)
        mip = pl if mip is None else np.maximum(mip, pl)
    lab2 = labels.max(axis=0)

    fig, axs = plt.subplots(1, 5, figsize=(26, 5.4))

    # 1 - the raw data
    axs[0].imshow(norm(mip), cmap="gray", interpolation="nearest")
    axs[0].set_title(f"1. actin max projection\n{field.field_id}", fontsize=9)

    # 2 - nuclei and what was traced in them
    g = norm(zoom(mip, (lab2.shape[0] / mip.shape[0], lab2.shape[1] / mip.shape[1]), order=1))
    rgb = np.dstack([g, g, g])
    edge = (lab2 > 0) & ~binary_erosion(lab2 > 0, np.ones((3, 3)))
    rgb[edge] = [0.2, 0.9, 1.0]
    axs[1].imshow(rgb, interpolation="nearest")
    for G in graphs:
        a = G.graph
        # Polylines are in FIELD microns; divide by the label spacing to draw on this panel.
        for _, _, e in G.edges(data=True):
            p = np.array([[float(v) for v in q.split(",")] for q in e["polyline"].split(";")])
            axs[1].plot(p[:, 2] / lgrid.dx_um, p[:, 1] / lgrid.dy_um, "-",
                        color="#ffe400", lw=0.6)
        # The label id, so a nucleus in the figure can be found in the `label` column.
        axs[1].text(a["nucleus_centroid_x_um"] / lgrid.dx_um,
                    a["nucleus_centroid_y_um"] / lgrid.dy_um, str(a["label"]),
                    color="#ff7043", fontsize=6, ha="center", va="center")
    axs[1].set_title("2. nuclei (cyan) + traced centrelines (yellow)\n"
                     "IS WHAT IT TRACED A FILAMENT?", fontsize=9)

    # 3 - where the threshold sits relative to the response
    seps = np.array([G.graph.get("separation", np.nan) for G in graphs], float)
    seps = seps[np.isfinite(seps)]
    if seps.size:
        axs[2].hist(np.log10(np.maximum(seps, 1e-6)), bins=18, color="#607d8b")
    axs[2].axvline(np.log10(params.detect.high_k), color="#d32f2f", lw=2,
                   label=f"high_k = {params.detect.high_k:g}")
    axs[2].axvline(np.log10(params.detect.low_k), color="#f9a825", lw=2,
                   label=f"low_k = {params.detect.low_k:g}")
    axs[2].set_xlabel("log10( in-mask response p99.9 / noise null )")
    axs[2].set_ylabel("nuclei")
    axs[2].legend(fontsize=7)
    axs[2].set_title("3. response vs threshold\nnuclei left of the red line cannot seed",
                     fontsize=9)

    # 4 - the slab, and how cut the nuclei are
    xz = labels.max(axis=1) > 0
    axs[3].imshow(xz, cmap="gray", aspect="auto", interpolation="nearest")
    axs[3].set_xlabel("x (label voxels)")
    axs[3].set_ylabel("z (planes)")
    ef = np.array([G.graph.get("axial_edge_fraction", np.nan) for G in graphs], float)
    med_ef = float(np.nanmedian(ef)) if np.isfinite(ef).any() else float("nan")
    axs[3].set_title(f"4. XZ side view - slab {field.slab_um:.1f} um\n"
                     f"median axial_edge_fraction {med_ef:.3f} "
                     f"({'CUT' if med_ef > 0.3 else 'closes in z'})", fontsize=9)

    # 5 - the outlier check
    vol = np.array([G.graph["nucleus_volume_um3"] for G in graphs], float)
    dens = np.array([G.graph["length_density_um_per_um3"] for G in graphs], float)
    axs[4].scatter(vol, dens, s=18, c="#37474f")
    axs[4].set_xlabel("nucleus volume (um3)")
    axs[4].set_ylabel("length density (um / um3)")
    axs[4].set_title("5. density vs volume\noutliers are puncta-driven nuclei", fontsize=9)

    for a in axs[:2]:
        a.axis("off")
    fig.tight_layout(pad=0.6)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def field_summary(ref, rows, seconds: float) -> dict:
    """The per-field QC figure, numerically. Sort by these to find the fields worth opening."""
    import pandas as pd
    df = pd.DataFrame(rows)

    def med(c):
        return float(df[c].median()) if len(df) and c in df else float("nan")

    return {
        "field_id": ref.field_id, "file": ref.stem, "position": ref.position,
        "condition": ref.condition, "seconds": round(seconds, 1),
        "n_nuclei": int(len(df)),
        "n_detected": int(df.detected.sum()) if len(df) else 0,
        "frac_detected": float(df.detected.mean()) if len(df) else float("nan"),
        "noise_sd_med": med("noise_sd"), "null_level_med": med("null_level"),
        # The single most diagnostic number: how far the brightest in-mask response sits above
        # the noise floor. Low means nothing can seed. HIGH DOES NOT MEAN IT IS FILAMENTS - on a
        # control packed with bright puncta this is larger than on real signal.
        "separation_med": med("separation"),
        "fg_fraction_med": med("fg_fraction"),
        "length_density_med": med("length_density_um_per_um3"),
        "total_length_med_um": med("total_antenna_length_um"),
        "n_puncta_rejected": int(df.n_puncta_rejected.sum()) if len(df) else 0,
        "branch_points_per_um_med": med("branch_points_per_um"),
        "n_gap_joins": int(df.n_gap_joins.sum()) if len(df) else 0,
        "n_crossings_resolved": int(df.n_crossings_resolved.sum()) if len(df) else 0,
        "axial_edge_fraction_med": med("axial_edge_fraction"),
        "cut_ends_frac": (float(df.n_cut_ends.sum() / max(df.n_ends_total.sum(), 1))
                          if len(df) else float("nan")),
    }
