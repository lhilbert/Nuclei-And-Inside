"""
Point-based figures from the measurement table.

    midplane_scatter    one panel per condition, every nucleus a point
    zclip_diagnostics   how badly the slab truncates the population --
                        read this BEFORE comparing conditions
    pca_summary         scree, loadings and scores for a feature_pca result

Image tiles live in `mosaic.py`; the feature decomposition in `features.py`.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from .style import CONTEXT_GREY, PALETTE, sizes

def pca_summary(result, out_png, conditions=None, condition_labels=None,
                n_loadings=None, base_fontsize=9, dpi=300, point_size=13):
    """
    Scree, loadings and scores for a `pca.feature_pca` result.

    Three panels, because a PCA is not interpretable from any one of them:
    how much variance each component carries, which features build the
    first two, and where the conditions sit in that plane. `n_loadings`
    limits the loading panel to the strongest features on PC1 (default:
    all of them).
    """
    ex, load = result["explained"], result["loadings"]
    scores = result["scores"]

    fs = sizes(base_fontsize)
    fig = plt.figure(figsize=(10.5, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.15], hspace=0.45,
                          wspace=0.32)
    ax_scree, ax_sc = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    ax_load = fig.add_subplot(gs[:, 1])

    n = len(ex)
    ax_scree.bar(range(1, n + 1), 100 * ex["variance_ratio"], color="0.55",
                 width=0.68)
    ax_scree.plot(range(1, n + 1), 100 * ex["cumulative"], color=PALETTE[0],
                  marker="o", ms=3.5, lw=1.2)
    for i, (v, c) in enumerate(zip(ex["variance_ratio"], ex["cumulative"]), 1):
        ax_scree.annotate(f"{100 * v:.0f}", xy=(i, 100 * v), xytext=(0, 2),
                          textcoords="offset points", ha="center",
                          fontsize=fs["annot"], color="0.25")
    ax_scree.set_xticks(range(1, n + 1))
    ax_scree.set_xticklabels(ex.index, fontsize=fs["tick"])
    ax_scree.set_ylabel("variance explained (%)", fontsize=fs["base"])
    ax_scree.set_title(f"{result['n_samples']} nuclei, "
                       f"{len(result['features'])} features",
                       fontsize=fs["base"], loc="left")
    ax_scree.annotate("cumulative", xy=(n, 100 * ex["cumulative"].iloc[-1]),
                      xytext=(-4, -12), textcoords="offset points", ha="right",
                      fontsize=fs["annot"], color=PALETTE[0])
    ax_scree.tick_params(labelsize=fs["tick"])
    for s in ("top", "right"):
        ax_scree.spines[s].set_visible(False)

    order = load["PC1"].abs().sort_values(ascending=False).index
    if n_loadings:
        order = order[:n_loadings]
    ypos = np.arange(len(order))
    ax_load.barh(ypos + 0.19, load.loc[order, "PC1"], height=0.36,
                 color=PALETTE[0], label=f"PC1 ({100 * ex['variance_ratio'].iloc[0]:.0f}%)")
    ax_load.barh(ypos - 0.19, load.loc[order, "PC2"], height=0.36,
                 color=PALETTE[1], label=f"PC2 ({100 * ex['variance_ratio'].iloc[1]:.0f}%)")
    ax_load.axvline(0, color="0.3", lw=0.7)
    ax_load.set_yticks(ypos)
    ax_load.set_yticklabels(order, fontsize=fs["tick"])
    ax_load.invert_yaxis()
    ax_load.set_xlabel("loading", fontsize=fs["base"])
    ax_load.tick_params(labelsize=fs["tick"])
    ax_load.legend(frameon=False, fontsize=fs["annot"], loc="lower right")
    for s in ("top", "right"):
        ax_load.spines[s].set_visible(False)

    if "condition" not in scores.columns:
        scores = scores.assign(condition="all")
    conds = list(conditions) if conditions else sorted(scores["condition"].unique())
    for i, c in enumerate(conds):
        sub = scores[scores["condition"] == c]
        lab = (condition_labels or {}).get(c, str(c).split("/")[-1])
        col = PALETTE[i % len(PALETTE)]
        ax_sc.scatter(sub["PC1"], sub["PC2"], s=point_size, c=col, alpha=0.7,
                      linewidths=0, label=f"{lab} (n={len(sub)})")
        ax_sc.plot(sub["PC1"].median(), sub["PC2"].median(), marker="+",
                   ms=10, mew=2.0, color=col)
    ax_sc.axhline(0, color="0.85", lw=0.7, zorder=0)
    ax_sc.axvline(0, color="0.85", lw=0.7, zorder=0)
    ax_sc.set_xlabel(f"PC1 ({100 * ex['variance_ratio'].iloc[0]:.0f}%)",
                     fontsize=fs["base"])
    ax_sc.set_ylabel(f"PC2 ({100 * ex['variance_ratio'].iloc[1]:.0f}%)",
                     fontsize=fs["base"])
    ax_sc.legend(frameon=False, fontsize=fs["annot"], loc="best")
    ax_sc.tick_params(labelsize=fs["tick"])
    for s in ("top", "right"):
        ax_sc.spines[s].set_visible(False)

    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    return dict(path=out_png, figure=fig,
                loading_order=list(order),
                medians={c: (float(scores.loc[scores["condition"] == c, "PC1"].median()),
                             float(scores.loc[scores["condition"] == c, "PC2"].median()))
                         for c in conds})


def zclip_diagnostics(table, out_png, n_z=None, condition_labels=None,
                      exclude_xy_border=True, base_fontsize=9, dpi=300):
    """
    How badly the slab thickness truncates the mid-plane read-out.

    Three panels: where each nucleus's widest plane sits in the stack, what
    fraction of each condition is clipped at a slab face, and whether being
    clipped biases the measurements.

    Run this before trusting a mid-plane comparison. Two things decide
    whether the read-out is usable: the clipping rate (nuclei clipped in z
    have no imaged equatorial section) and whether that rate differs between
    conditions -- if it does, filtering on `mid_at_z_border` leaves the
    conditions sampled differently, which is a selection effect rather than
    a biological one.

    Returns a dict with the per-condition rates and the clipped/inside
    median ratio for each metric, so the numbers can go in a caption.
    """
    df = table.copy()
    if "condition" not in df.columns:
        df["condition"] = "all"
    if exclude_xy_border and "touches_xy_border" in df.columns:
        # isolate the z effect: a laterally cut nucleus is a fragment, and
        # letting fragments into the comparison would mix the two
        # truncations rather than measure either
        df = df[~df["touches_xy_border"]]
    df["cond"] = df["condition"].map(
        lambda c: (condition_labels or {}).get(c, str(c).split("/")[-1]))
    if n_z is None:
        n_z = int(df["mid_z"].max()) + 1

    mets = [("max_area_um2", "mid-plane area"),
            ("mid_dna_cv_corr", "CV"),
            ("mid_solidity", "solidity"),
            ("mid_dna_radial_norm", "radial index")]
    fs = sizes(base_fontsize)
    fig, axs = plt.subplots(1, 3, figsize=(11.4, 3.1),
                            gridspec_kw=dict(width_ratios=[1.15, 1.0, 1.35]))

    # --- (a) where the widest plane sits -----------------------------------
    ax = axs[0]
    counts = df["mid_z"].value_counts().reindex(range(n_z)).fillna(0)
    faces = [0, n_z - 1]
    ax.bar([z for z in range(n_z) if z not in faces],
           [counts[z] for z in range(n_z) if z not in faces],
           color=CONTEXT_GREY, width=0.9)
    ax.bar(faces, [counts[z] for z in faces], color="#D55E00", width=0.9)
    ax.set_xlabel("z-plane of widest cross-section", fontsize=fs["base"])
    ax.set_ylabel("nuclei", fontsize=fs["base"])
    ax.annotate("clipped at\nslab face", xy=(0, counts[0]),
                xytext=(n_z * 0.22, counts.max() * 0.86),
                fontsize=fs["annot"], color="#D55E00",
                arrowprops=dict(arrowstyle="-", color="#D55E00", lw=0.8))

    # --- (b) clipping rate per condition -----------------------------------
    ax = axs[1]
    rate = (df.groupby("cond")["mid_at_z_border"]
              .agg(["size", "sum"]).sort_values("size", ascending=False))
    rate["pct"] = 100 * rate["sum"] / rate["size"]
    ypos = np.arange(len(rate))
    ax.hlines(ypos, 0, rate["pct"], color=CONTEXT_GREY, lw=1.2)
    ax.plot(rate["pct"], ypos, "o", color="#D55E00", ms=6)
    for y, (pct, n) in enumerate(zip(rate["pct"], rate["size"])):
        ax.annotate(f"{pct:.0f}%  (n={n})", (pct, y), xytext=(6, 0),
                    textcoords="offset points", va="center",
                    fontsize=fs["annot"])
    ax.set_yticks(ypos); ax.set_yticklabels(rate.index, fontsize=fs["tick"])
    ax.set_xlabel("nuclei clipped in z (%)", fontsize=fs["base"])
    ax.set_xlim(0, max(100, rate["pct"].max() * 1.55))
    # keep the lowest marker's value label clear of the x-axis title
    ax.set_ylim(-0.65, len(rate) - 0.35)

    # --- (c) does clipping bias the metrics? -------------------------------
    ax = axs[2]
    ratios = {}
    rng = np.random.default_rng(0)
    for i, (col, name) in enumerate(mets):
        inside = df.loc[~df["mid_at_z_border"], col].dropna()
        clipped = df.loc[df["mid_at_z_border"], col].dropna()
        ref = float(np.median(inside))
        ratios[col] = float(np.median(clipped) / ref) if ref else np.nan
        for j, (vals, colour) in enumerate([(inside, "#0072B2"),
                                            (clipped, "#D55E00")]):
            x = i + (j - 0.5) * 0.34 + rng.normal(0, 0.035, len(vals))
            ax.scatter(x, vals / ref, s=7, c=colour, alpha=0.55, linewidths=0,
                       rasterized=True)
            ax.plot([i + (j - 0.5) * 0.34 - 0.1, i + (j - 0.5) * 0.34 + 0.1],
                    [np.median(vals) / ref] * 2, color="black", lw=1.6)
    ax.axhline(1.0, color="black", lw=0.6, ls=":", zorder=0)
    ax.set_xticks(range(len(mets)))
    ax.set_xticklabels([m[1] for m in mets], fontsize=fs["tick"])
    ax.set_ylabel("value / median of unclipped", fontsize=fs["base"])
    ax.scatter([], [], s=18, c="#0072B2", label="inside slab")
    ax.scatter([], [], s=18, c="#D55E00", label="z-clipped")
    ax.legend(frameon=False, fontsize=fs["annot"], loc="upper right",
              handletextpad=0.2, borderpad=0.1)

    for a in axs:
        a.tick_params(labelsize=fs["tick"])
        for side in ("top", "right"):
            a.spines[side].set_visible(False)

    fig.tight_layout(w_pad=1.4)
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")

    return dict(path=out_png, figure=fig,
                clipped_pct=rate["pct"].round(1).to_dict(),
                n_by_condition=rate["size"].to_dict(),
                clipped_over_inside=ratios)


def midplane_scatter(table, out_png, metrics=None, conditions=None,
                     condition_labels=None, exclude_z_border=True,
                     exclude_xy_border=True, base_fontsize=9, dpi=300,
                     point_size=13):
    """
    One panel per condition, every nucleus a point.

    Parameters
    ----------
    table : DataFrame
        `nuclei_measurements.csv`, or the in-memory equivalent.
    out_png : output path; parent directories are created.
    metrics : list of (x_column, y_column, x_label, y_label), one row of
        panels each. Defaults to CV vs solidity, then CV vs radial index.
    conditions : order of the condition panels; default sorted unique.
    exclude_z_border : drop nuclei whose widest plane is the first or last
        plane of the stack. Those nuclei are cut off in z, so their widest
        imaged section is a grazing cut rather than an equatorial one and
        their solidity and CV are not comparable with the rest. Counts of
        what was dropped are returned.
    exclude_xy_border : drop nuclei whose mask reaches an image edge in x or
        y (`touches_xy_border`). These are cut off laterally, which is the
        larger of the two truncations: on the example data their median
        mid-plane area is 0.61x and their volume 0.58x that of interior
        nuclei, and their solidity, CV and radial index all differ
        significantly too. A laterally cut nucleus is a fragment, so every
        shape and texture read-out on it describes the fragment.

    Returns
    -------
    dict with the output path and per-condition counts kept/dropped.
    """
    if metrics is None:
        metrics = [
            ("mid_dna_cv_corr", "mid_solidity",
             "CV of DNA intensity", "Solidity of DNA mask"),
            ("mid_dna_cv_corr", "mid_dna_radial_norm",
             "CV of DNA intensity", "Intensity-weighted radius / mean radius"),
            # CV stays on x in every row so a column reads as one nucleus
            # population seen against three different shape/size read-outs.
            ("mid_dna_cv_corr", "max_area_um2",
             "CV of DNA intensity", "Mid-plane area (µm²)"),
        ]

    df = table.copy()
    if "condition" not in df.columns:
        df["condition"] = "all"

    n_before = df.groupby("condition").size()
    dropped, dropped_xy = {}, {}
    if exclude_z_border and "mid_at_z_border" in df.columns:
        dropped = df[df["mid_at_z_border"]].groupby("condition").size().to_dict()
        df = df[~df["mid_at_z_border"]]
    if exclude_xy_border and "touches_xy_border" in df.columns:
        dropped_xy = (df[df["touches_xy_border"]].groupby("condition")
                      .size().to_dict())
        df = df[~df["touches_xy_border"]]
    if df.empty:
        raise ValueError("every nucleus was excluded; relax the border filters")

    conds = list(conditions) if conditions else sorted(df["condition"].unique())
    colours = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(conds)}
    fs = sizes(base_fontsize)

    nrows, ncols = len(metrics), len(conds)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.35 * ncols + 0.6,
                                                    2.5 * nrows + 0.5),
                             sharex="row", sharey="row", squeeze=False)

    for ri, (xcol, ycol, xlab, ylab) in enumerate(metrics):
        allx, ally = df[xcol].to_numpy(), df[ycol].to_numpy()

        for ci, cond in enumerate(conds):
            ax = axes[ri][ci]
            sub = df[df["condition"] == cond]

            # every other nucleus in grey, so each panel is read against the
            # whole dataset rather than against its own axis limits
            ax.scatter(allx, ally, s=point_size * 0.7, c=CONTEXT_GREY,
                       linewidths=0, zorder=1, rasterized=True)
            ax.scatter(sub[xcol], sub[ycol], s=point_size, c=colours[cond],
                       linewidths=0, alpha=0.85, zorder=3, rasterized=True)

            # median of this condition, drawn as a glyph that cannot be
            # mistaken for a nucleus
            mx, my = np.nanmedian(sub[xcol]), np.nanmedian(sub[ycol])
            ax.plot(mx, my, marker="+", ms=11, mew=2.0, color="black", zorder=5)

            if ri == 0:
                # condition folders are often long paths; show the leaf name
                short = (condition_labels or {}).get(cond, cond.split("/")[-1])
                ax.set_title(f"{short}\nn = {len(sub)}", fontsize=fs["base"],
                             loc="left")
            if ci == 0:
                ax.set_ylabel(ylab, fontsize=fs["base"])
            ax.set_xlabel(xlab, fontsize=fs["base"])
            ax.tick_params(labelsize=fs["tick"])
            ax.margins(0.06)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

    # one key for the median glyph, in the whitespace of the first panel
    axes[0][0].plot([], [], marker="+", ms=9, mew=2.0, color="black", ls="none",
                    label="condition median")
    axes[0][0].legend(frameon=False, fontsize=fs["annot"], loc="lower left",
                      handletextpad=0.3, borderpad=0.1)

    fig.tight_layout(w_pad=0.6, h_pad=1.0)
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")

    return dict(path=out_png, figure=fig,
                kept={c: int((df["condition"] == c).sum()) for c in conds},
                dropped_z_border=dropped,
                dropped_xy_border=dropped_xy,
                n_before={k: int(v) for k, v in n_before.items()})
