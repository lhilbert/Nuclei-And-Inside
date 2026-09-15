"""
Image galleries and mosaics laid out on the measurement plane.

Both figures here answer "what do the nuclei at this coordinate actually
look like", and they differ only in density: `nucleus_gallery` shows a
dozen captioned examples beside a scatter, `nucleus_mosaic` composites
many tiles into one image whose axes are the metrics. The three steps they
share -- select drawable nuclei, pick one per cell without replacement,
crop mid-plane tiles from the raw stacks -- are the private helpers below.

Crops are read from the .nd2 files on demand via `source_path`, `position`
and `mid_z`, so the image data must be reachable; each field is opened once
however many nuclei it contributes. No mask outline is drawn: the label
image is not stored in the table and re-segmenting a field to draw one
costs ~40 s.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from ..core.io import load_field
from .style import CONTEXT_GREY, HIGHLIGHT, sizes

TILE_MARGIN = 1.2
"""
Tile size as a multiple of the largest drawn nucleus's diameter.

All tiles must share one physical scale to stay comparable, so the window
is set by the largest nucleus in the selection. A fixed window looks tidy
until a nucleus exceeds it, and then the tile clips it in exactly the way
the xy-border filter exists to avoid -- which reads as a filtering failure
when it is only a framing one.
"""


def auto_window_um(df, margin=TILE_MARGIN, floor=12.0):
    """Tile size that cannot cut the nuclei being drawn."""
    if "max_area_um2" not in df.columns or df["max_area_um2"].dropna().empty:
        return floor
    diam = 2.0 * np.sqrt(df["max_area_um2"].max() / np.pi)
    return float(max(floor, np.ceil(margin * diam)))


def _drawable(table, needed, condition=None, exclude_z_border=True,
              exclude_xy_border=True):
    """
    Nuclei that can be drawn: filtered, and with the columns a crop needs.

    `exclude_xy_border` drops nuclei whose mask reaches an image edge. That
    is the larger of the two truncations -- on the example data their median
    mid-plane area is 0.61x and their volume 0.58x that of interior nuclei
    -- and a laterally cut nucleus is a fragment, so every shape and texture
    read-out on it describes the fragment. `exclude_z_border` drops nuclei
    whose widest plane is a slab face, whose widest section was never
    imaged.
    """
    df = table.copy()
    if condition is not None:
        df = df[df["condition"].astype(str).str.contains(condition)]
    if exclude_z_border and "mid_at_z_border" in df.columns:
        df = df[~df["mid_at_z_border"]]
    if exclude_xy_border and "touches_xy_border" in df.columns:
        df = df[~df["touches_xy_border"]]
    df = df.dropna(subset=list(needed) + ["mid_z", "centroid_y_um",
                                          "centroid_x_um"])
    if df.empty:
        raise ValueError("no nuclei left to draw after filtering")
    return df


def _pick_one_per_cell(df, x, y, xedges, yedges, tol):
    """
    Nearest nucleus to each cell centre, WITHOUT replacement.

    Distances are in units of each cell's own size, so the two metrics
    weigh equally whatever their ranges, and uneven (quantile) cells are
    handled correctly. Without replacement because `tol` lets one nucleus
    be nearest to two adjacent centres, and the same nucleus drawn twice
    wastes a tile and reads as two examples. A cell further than `tol` from
    every nucleus stays empty, which is itself informative: that
    combination does not occur in the data.
    """
    nrow, ncol = len(yedges) - 1, len(xedges) - 1
    picks, used = {}, set()
    for r in range(nrow):
        for c in range(ncol):
            avail = df.drop(index=list(used), errors="ignore")
            if avail.empty:
                return picks
            xc, yc = xedges[c:c + 2].mean(), yedges[r:r + 2].mean()
            wx = max(xedges[c + 1] - xedges[c], 1e-12)
            wy = max(yedges[r + 1] - yedges[r], 1e-12)
            d = np.hypot((avail[x] - xc) / wx, (avail[y] - yc) / wy)
            if d.min() <= tol:
                idx = d.idxmin()
                picks[(r, c)] = avail.loc[idx]
                used.add(idx)
    return picks


def _crop_tiles(picks, window_um, dna_channel, percentiles=None):
    """
    Mid-plane crop per picked nucleus, grouped so each field is read once.

    Every tile gets the same pixel dimensions and a crop running off the
    image edge is zero-padded rather than trimmed, so the tiles form a
    regular grid at one physical scale. With `percentiles`, each tile is
    contrast-stretched to 0..1 individually -- right for judging chromatin
    texture, wrong for comparing brightness between tiles.

    Returns (tiles, dy, (ty, tx)).
    """
    by_field = {}
    for key, row in picks.items():
        by_field.setdefault((row["source_path"], int(row["position"])), []).append(key)

    tiles, dy_out, ty, tx = {}, None, None, None
    for (path, pos), keys in by_field.items():
        stack = load_field(path, pos)
        _, dy, dx = stack.voxel_um
        vol = stack.channel(dna_channel)
        if ty is None:
            ty, tx = int(round(window_um / dy)), int(round(window_um / dx))
            dy_out = dy
        for key in keys:
            row = picks[key]
            y0 = int(round(row["centroid_y_um"] / dy - ty / 2))
            x0 = int(round(row["centroid_x_um"] / dx - tx / 2))
            tile = np.zeros((ty, tx), dtype=np.float32)
            sy0, sx0 = max(0, y0), max(0, x0)
            sy1, sx1 = min(vol.shape[1], y0 + ty), min(vol.shape[2], x0 + tx)
            tile[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = \
                vol[int(row["mid_z"]), sy0:sy1, sx0:sx1]
            if percentiles is not None:
                lo, hi = np.percentile(tile, percentiles)
                tile = np.clip((tile - lo) / max(hi - lo, 1e-9), 0, 1)
            tiles[key] = tile
        del stack, vol
    return tiles, dy_out, (ty, tx)


def nucleus_gallery(table, out_png, x="mid_dna_cv_corr", y="mid_solidity",
                    grid=(3, 4), window_um=None, dna_channel="DAPI",
                    condition=None, exclude_z_border=True,
                    exclude_xy_border=True, percentiles=(1, 99.5),
                    base_fontsize=9, dpi=300):
    """
    What does a nucleus in each region of the scatter actually look like?

    Tiles the `x`-`y` plane into `grid` = (rows, cols) cells over the 2nd to
    98th percentile range, picks the nucleus closest to each cell centre,
    and shows its mid-plane xy section next to a scatter marking which
    nuclei were picked. The tiles are laid out in the same geometry as the
    scatter -- x increasing rightwards, y increasing upwards -- so the
    gallery reads as the plane itself.
    """
    df = _drawable(table, (x, y), condition, exclude_z_border,
                   exclude_xy_border)
    window_um = auto_window_um(df) if window_um is None else window_um

    nrow, ncol = grid
    xedges = np.linspace(*np.percentile(df[x], [2, 98]), ncol + 1)
    yedges = np.linspace(*np.percentile(df[y], [2, 98]), nrow + 1)
    picks = _pick_one_per_cell(df, x, y, xedges, yedges, tol=0.75)
    crops, dy, _ = _crop_tiles(picks, window_um, dna_channel)

    fs = sizes(base_fontsize)
    fig = plt.figure(figsize=(3.6 + 1.55 * ncol, 1.0 + 1.55 * nrow))
    gs = fig.add_gridspec(nrow, ncol + 3, wspace=0.08, hspace=0.08)
    ax_s = fig.add_subplot(gs[:, :3])

    ax_s.scatter(df[x], df[y], s=10, c=CONTEXT_GREY, linewidths=0,
                 rasterized=True)
    for (r, c), row in picks.items():
        ax_s.scatter([row[x]], [row[y]], s=42, facecolors="none",
                     edgecolors=HIGHLIGHT, linewidths=1.3)
        ax_s.annotate(f"{chr(ord('a') + r * ncol + c)}", (row[x], row[y]),
                      xytext=(4, 3), textcoords="offset points",
                      fontsize=fs["annot"], color=HIGHLIGHT)
    for e in xedges:
        ax_s.axvline(e, color="black", lw=0.3, alpha=0.25, zorder=0)
    for e in yedges:
        ax_s.axhline(e, color="black", lw=0.3, alpha=0.25, zorder=0)
    ax_s.set_xlabel(x, fontsize=fs["base"])
    ax_s.set_ylabel(y, fontsize=fs["base"])
    ax_s.tick_params(labelsize=fs["tick"])
    for side in ("top", "right"):
        ax_s.spines[side].set_visible(False)

    for r in range(nrow):
        for c in range(ncol):
            ax = fig.add_subplot(gs[nrow - 1 - r, c + 3])
            ax.set_xticks([]); ax.set_yticks([])
            letter = chr(ord("a") + r * ncol + c)
            if (r, c) not in crops:
                ax.text(0.5, 0.5, f"{letter}\nno nuclei", ha="center",
                        va="center", color=CONTEXT_GREY, fontsize=fs["tick"],
                        linespacing=1.5, transform=ax.transAxes)
                for side in ax.spines.values():
                    side.set_color(CONTEXT_GREY)
                    side.set_linestyle((0, (2, 2)))
                continue
            img = crops[(r, c)]
            lo, hi = np.percentile(img, percentiles)
            ax.imshow(img, cmap="gray", vmin=lo, vmax=hi,
                      interpolation="nearest")
            row = picks[(r, c)]
            ax.set_title(f"{letter}   {row[x]:.2f} / {row[y]:.3f}",
                         fontsize=fs["tick"], pad=2, loc="left")
            if (r, c) == min(crops):              # one scale bar for the set
                px = 5.0 / dy
                ax.plot([img.shape[1] - px - 4, img.shape[1] - 4],
                        [img.shape[0] - 5] * 2, color="white", lw=2.0)
                ax.annotate("5 µm", (img.shape[1] - px / 2 - 4,
                                     img.shape[0] - 8), color="white",
                            fontsize=fs["tick"], ha="center", va="bottom")

    fig.suptitle(f"tile label: {x} / {y}", fontsize=fs["annot"], x=0.99,
                 ha="right", y=0.995)
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")

    return dict(path=out_png, figure=fig, window_um=window_um,
                picked={f"{chr(ord('a') + r * ncol + c)}":
                        dict(nucleus_uid=row["nucleus_uid"],
                             **{x: float(row[x]), y: float(row[y])})
                        for (r, c), row in picks.items()},
                n_cells=nrow * ncol, n_filled=len(crops))


def nucleus_mosaic(table, out_png, x="mid_dna_cv_corr", y="mid_solidity",
                   grid=(6, 9), window_um=None, dna_channel="DAPI",
                   condition=None, exclude_z_border=True,
                   exclude_xy_border=True, percentiles=(1, 99.5),
                   binning="quantile", sort_by=None,
                   xlabel=None, ylabel=None, base_fontsize=9, dpi=200):
    """
    A dense image mosaic laid out ON the measurement plane.

    Where the gallery shows a dozen captioned examples, this shows the
    continuum -- the layout to use when the plane appears to order nuclei
    along something monotonic, such as chromatin condensation.

    The tiles are pasted into one array and drawn with a single `imshow`
    under `extent`, not as one Axes per cell: sixty-odd Axes would each
    carry their own transform and margins, whereas one array makes the
    mosaic a real image, puts the axis ticks on the metric values by
    construction, and leaves empty cells (drawn white) exactly where the
    plane is unpopulated.

    `binning`
        "quantile" (default) puts cell edges at equal-count quantiles, so
        the mosaic is dense where the nuclei are; the tick labels carry the
        real values at those edges, but axis SPACING is then rank, not
        value. "linear" keeps evenly spaced edges and an honest plane at
        the cost of empty cells -- on the example data a linear 6x9 grid
        fills 30 of 54 cells against 52 for quantile, because the nuclei
        occupy a narrow diagonal band.

    `sort_by`
        A column name gives a one-dimensional mosaic instead: nuclei ranked
        by that column, the grid filled in reading order with an even
        sample across the whole range. `x`/`y` are then ignored.
    """
    df = _drawable(table, (x, y) if sort_by is None else (sort_by,),
                   condition, exclude_z_border, exclude_xy_border)
    window_um = auto_window_um(df) if window_um is None else window_um
    nrow, ncol = grid

    if sort_by is not None:
        ranked = df.sort_values(sort_by)
        take = np.unique(np.linspace(0, len(ranked) - 1, nrow * ncol).astype(int))
        picks = {(nrow - 1 - n // ncol, n % ncol): ranked.iloc[p]
                 for n, p in enumerate(take)}
        sorted_vals = ranked[sort_by].to_numpy()[take]
        xedges, yedges = np.arange(ncol + 1), np.arange(nrow + 1)
    else:
        sorted_vals = None
        if binning == "quantile":
            xedges = np.quantile(df[x], np.linspace(0, 1, ncol + 1))
            yedges = np.quantile(df[y], np.linspace(0, 1, nrow + 1))
        else:
            xedges = np.linspace(*np.percentile(df[x], [1, 99]), ncol + 1)
            yedges = np.linspace(*np.percentile(df[y], [1, 99]), nrow + 1)
        picks = _pick_one_per_cell(df, x, y, xedges, yedges, tol=0.7)

    tiles, _, (ty, tx) = _crop_tiles(picks, window_um, dna_channel,
                                     percentiles=percentiles)

    canvas = np.full((nrow * ty, ncol * tx), np.nan, dtype=np.float32)
    for (r, c), tile in tiles.items():
        r_img = nrow - 1 - r                       # row 0 at the bottom
        canvas[r_img * ty:(r_img + 1) * ty, c * tx:(c + 1) * tx] = tile

    fs = sizes(base_fontsize)
    fig, ax = plt.subplots(figsize=(ncol * 1.05, nrow * 1.05))
    cmap = plt.get_cmap("gray").copy()
    cmap.set_bad("white")

    # Tiles occupy equal shares of the canvas, so the axes must be in CELL
    # units whenever the cell edges are not equally spaced -- otherwise the
    # extent stretches the canvas linearly over the value range and the grid
    # lines no longer coincide with the tile boundaries they mark. Quantile
    # and sorted modes therefore plot in rank space and carry the real
    # values as tick labels; linear mode can use the values directly.
    rank_space = (sort_by is not None) or (binning == "quantile")
    if rank_space:
        extent = (0, ncol, 0, nrow)
        xgrid, ygrid = np.arange(1, ncol), np.arange(1, nrow)
    else:
        extent = (xedges[0], xedges[-1], yedges[0], yedges[-1])
        xgrid, ygrid = xedges[1:-1], yedges[1:-1]

    ax.imshow(canvas, cmap=cmap, vmin=0, vmax=1, aspect="auto",
              interpolation="nearest", extent=extent, origin="upper")
    for e in xgrid:
        ax.axvline(e, color="white", lw=0.5, alpha=0.55)
    for e in ygrid:
        ax.axhline(e, color="white", lw=0.5, alpha=0.55)

    if sort_by is not None:
        # label the first tile of each row, so the ramp is readable without
        # one label per tile
        ax.set_xticks([])
        rows = range(0, len(sorted_vals), ncol)
        ax.set_yticks([nrow - 1 - n // ncol + 0.5 for n in rows])
        ax.set_yticklabels([f"{sorted_vals[n]:.2f}" for n in rows])
        ax.set_ylabel(ylabel or f"{sort_by} (row start)", fontsize=fs["base"])
        ax.set_xlabel(xlabel or f"ranked by {sort_by}, left to right",
                      fontsize=fs["base"])
    elif rank_space:
        # ticks on cell boundaries, labels the quantile edges. No reversal
        # on y: canvas row 0 holds grid row nrow-1 and origin="upper" puts
        # it at y = nrow, so ascending positions carry ascending edges.
        stepx, stepy = (2 if ncol > 6 else 1), (2 if nrow > 6 else 1)
        ax.set_xticks(np.arange(0, ncol + 1, stepx))
        ax.set_xticklabels([f"{v:.2f}" for v in xedges[::stepx]])
        ax.set_yticks(np.arange(0, nrow + 1, stepy))
        ax.set_yticklabels([f"{v:.3f}" for v in yedges[::stepy]])
        ax.set_xlabel(xlabel or f"{x} (equal-count cells)", fontsize=fs["base"])
        ax.set_ylabel(ylabel or f"{y} (equal-count cells)", fontsize=fs["base"])
    else:
        ax.set_xticks(np.round(xedges, 2))
        ax.set_yticks(np.round(yedges, 3))
        ax.set_xlabel(xlabel or x, fontsize=fs["base"])
        ax.set_ylabel(ylabel or y, fontsize=fs["base"])
    ax.tick_params(labelsize=fs["tick"], length=2)
    ax.annotate(f"tile {window_um:.0f} µm", xy=(0.995, 1.004),
                xycoords="axes fraction", ha="right", va="bottom",
                fontsize=fs["annot"], color="0.35")

    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")

    return dict(path=out_png, figure=fig, grid=(nrow, ncol),
                window_um=window_um, n_cells=nrow * ncol, n_filled=len(tiles),
                picked={f"r{r}c{c}": row["nucleus_uid"]
                        for (r, c), row in picks.items()})
