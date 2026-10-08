"""
An angular pseudo-time coordinate on a two-metric plane.

The cell cycle is closed, so if two metrics track it the nuclei should trace
a loop in the plane they span. Place a centre inside that loop and the angle
of each nucleus about it orders the cells -- the same construction `tricycle`
and Revelio use on cell-cycle gene expression, here on chromatin morphology.

The plane this was built for is `mid_dna_cv_corr` against
`mid_dna_persistence`: chromatin contrast and how much prominent internal
structure there is. Those two are essentially uncorrelated (r = 0.05), which
is what makes them worth spanning a plane with.

THE TRAP, and it cost an hour the first time
--------------------------------------------
An angle is a point on a circle, not a number on a line. A group of nuclei
sitting either side of theta = 0 has a perfectly well-defined centre -- at
zero -- but `numpy.median` of their angles returns something near pi,
halfway round the circle from where every one of them actually is. On the
validation data that mistake made theta look like it separated nothing; the
circular mean showed it separates condensed from interphase at p = 5e-9.

Use `circular_mean` and `circular_R` here. Never `mean`, `median` or `std`
on a column of angles.

Is it really a ring?
--------------------
Worth checking before trusting the coordinate, because an angle is defined
whether or not the geometry supports it: for a blob the angle is noise with
units. A ring has nuclei at a comparable radius all the way round; an arc
has them bunched over part of the range, and then outbound and return legs
share angles and theta is ambiguous between them.

On the validation data the condensed nuclei sit at median radius 1.14
against 2.09 for interphase -- closer to the centre rather than out at a
matching radius, which points to an arc rather than a closed ring. The
radius is therefore returned alongside theta: use it as a confidence, and
treat nuclei near the centre as unplaced rather than as having a small
angle.
"""

import numpy as np

from .style import sizes, HIGHLIGHT


def circular_mean(angles):
    """
    Mean direction of a set of angles, in [0, 2pi).

    Averages the unit vectors rather than the numbers, which is the only
    correct way: see the module docstring for what happens otherwise.
    """
    a = np.asarray(angles, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    return float(np.mod(np.arctan2(np.sin(a).mean(), np.cos(a).mean()), 2 * np.pi))


def circular_R(angles):
    """
    Resultant length: 1 when the angles coincide, 0 when spread uniformly.

    The angular analogue of 1 - variance, and the number to quote when
    claiming a group is concentrated at some phase.
    """
    a = np.asarray(angles, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    return float(np.hypot(np.cos(a).mean(), np.sin(a).mean()))


def angular_distance(angles, reference):
    """Shortest angular separation from `reference`, in [0, pi]."""
    a = np.asarray(angles, dtype=float)
    return np.abs(np.mod(a - reference + np.pi, 2 * np.pi) - np.pi)


def angular_pseudotime(table, centre, pi_toward, x="mid_dna_cv_corr",
                       y="mid_dna_persistence", direction="ccw",
                       standardize=True):
    """
    Add `theta` and `radius` columns: the angle of each nucleus about a centre.

    Parameters
    ----------
    table : DataFrame
        Per-nucleus measurements.
    centre : (cx, cy)
        Rotation centre, in the SAME units as the plane -- standardised
        units when `standardize` is true, which is the default and the
        recommended way to choose one from a plot.
    pi_toward : (px, py)
        A point such that theta = pi points from the centre toward it. This
        fixes the origin of the coordinate, which is otherwise arbitrary.
        Pick it so that a phase you can recognise lands somewhere memorable.
    x, y : str
        Columns spanning the plane.
    direction : {"ccw", "cw"}
        Which way theta increases. Neither is more correct; choose the one
        that puts the phases you can identify in their biological order,
        and say which you used.
    standardize : bool
        Z-score both columns first. Keep this ON. `atan2` mixes the two
        axes, so without it the angle is dominated by whichever column
        carries larger numbers -- on the validation data persistence spans
        0-20 and contrast spans 0.13-0.62, so the raw angle would be a
        persistence read-out with a contrast-shaped wobble.

    Returns
    -------
    (DataFrame, dict)
        A copy of `table` with `theta` (radians, [0, 2pi)) and `radius`
        added, and the parameters needed to reproduce or invert the
        mapping -- including the mean and sd used, so a centre chosen in
        standardised units can be reported in raw ones.
    """
    if direction not in ("ccw", "cw"):
        raise ValueError("direction must be 'ccw' or 'cw'")
    for c in (x, y):
        if c not in table.columns:
            raise KeyError(f"{c} not in table")

    out = table.copy()
    xv = out[x].to_numpy(dtype=float)
    yv = out[y].to_numpy(dtype=float)

    if standardize:
        mu = (float(np.nanmean(xv)), float(np.nanmean(yv)))
        sd = (float(np.nanstd(xv, ddof=1)), float(np.nanstd(yv, ddof=1)))
        xv = (xv - mu[0]) / sd[0]
        yv = (yv - mu[1]) / sd[1]
    else:
        mu, sd = (0.0, 0.0), (1.0, 1.0)

    cx, cy = float(centre[0]), float(centre[1])
    dx, dy = xv - cx, yv - cy

    raw = np.arctan2(dy, dx)
    ref = float(np.arctan2(float(pi_toward[1]) - cy, float(pi_toward[0]) - cx))
    turn = (raw - ref) if direction == "ccw" else -(raw - ref)

    out["theta"] = np.mod(turn + np.pi, 2 * np.pi)
    out["radius"] = np.hypot(dx, dy)

    params = dict(x=x, y=y, centre=(cx, cy),
                  pi_toward=(float(pi_toward[0]), float(pi_toward[1])),
                  direction=direction, standardize=bool(standardize),
                  mean=mu, sd=sd, reference_angle=ref,
                  centre_raw=(cx * sd[0] + mu[0], cy * sd[1] + mu[1]))
    return out, params


def pseudotime_geometry(table, out_png, params, highlight=None,
                        highlight_label="highlighted", base_fontsize=9,
                        dpi=200):
    """
    The plane, coloured by theta, with the centre and the theta = pi ray drawn.

    Look at this before using the coordinate. It is where an arc shows
    itself as an arc, and where a centre placed outside the cloud -- which
    makes every angle crowd into a narrow fan -- is obvious at a glance.

    `highlight` is an optional boolean mask (for example the nuclei you
    have independently scored) drawn as open rings on top.
    """
    import os

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fs = sizes(base_fontsize)
    x, y = params["x"], params["y"]
    mu, sd = params["mean"], params["sd"]
    cx, cy = params["centre"]
    px, py = params["pi_toward"]
    ref = params["reference_angle"]

    xv = (table[x].to_numpy(float) - mu[0]) / sd[0]
    yv = (table[y].to_numpy(float) - mu[1]) / sd[1]

    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    sc = ax.scatter(xv, yv, c=table["theta"], cmap="twilight", s=28,
                    vmin=0, vmax=2 * np.pi, edgecolor="none", zorder=2)
    if highlight is not None:
        m = np.asarray(highlight, dtype=bool)
        ax.scatter(xv[m], yv[m], s=92, facecolor="none", edgecolor=HIGHLIGHT,
                   linewidth=1.5, zorder=4, label=highlight_label)
        ax.legend(loc="upper left", frameon=False, fontsize=fs["annot"])

    ax.plot(cx, cy, "k+", ms=13, mew=1.8, zorder=5)
    ax.annotate("centre", (cx, cy), xytext=(9, 9), textcoords="offset points",
                fontsize=fs["annot"])
    ax.annotate("", xy=(px, py), xytext=(cx, cy), zorder=5,
                arrowprops=dict(arrowstyle="->", lw=1.5, color="k"))
    ax.annotate(r"$\theta=\pi$", (px, py), xytext=(8, -12),
                textcoords="offset points", fontsize=fs["base"])

    span = 1.15 * float(np.nanmax(np.hypot(xv - cx, yv - cy)))
    sign = 1.0 if params["direction"] == "ccw" else -1.0
    for t, lab in [(0.0, r"$\theta=0$"), (np.pi / 2, r"$\theta=\pi/2$"),
                   (3 * np.pi / 2, r"$\theta=3\pi/2$")]:
        a = ref + sign * (np.pi - t)
        ax.plot([cx, cx + span * np.cos(a)], [cy, cy + span * np.sin(a)],
                color="0.7", lw=0.7, ls=":", zorder=1)
        ax.annotate(lab, (cx + span * np.cos(a), cy + span * np.sin(a)),
                    fontsize=fs["annot"], color="0.45")

    ax.set_xlabel(f"{x}  (z)", fontsize=fs["base"])
    ax.set_ylabel(f"{y}  (z)", fontsize=fs["base"])
    ax.set_title("angular pseudo-time about the chosen centre", fontsize=fs["base"])
    ax.set_aspect("equal")
    ax.tick_params(labelsize=fs["tick"])
    fig.colorbar(sc, ax=ax, label=r"$\theta$ (rad)", shrink=0.85)

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=dpi)
    return dict(path=out_png, figure=fig, n=len(table))


def phase_order_check(table, phase_column, order, theta_column="theta"):
    """
    Circular mean theta per labelled phase, and whether it increases in
    `order`.

    This is the validation worth doing: if the labels were used only to
    ORIENT the coordinate (choosing `pi_toward` and `direction`), their
    ordering along theta is a free prediction and a violation is
    informative. Counts are returned because the test is only as good as
    they are -- two cells cannot establish an ordering.
    """
    rows = []
    for p in order:
        sel = table[table[phase_column] == p]
        if len(sel) == 0:
            continue
        rows.append(dict(phase=p, n=int(len(sel)),
                         theta_circular_mean=circular_mean(sel[theta_column]),
                         R=circular_R(sel[theta_column]),
                         radius_median=float(sel["radius"].median())
                         if "radius" in sel else np.nan))
    means = [r["theta_circular_mean"] for r in rows]
    return dict(per_phase=rows,
                monotonic=all(b > a for a, b in zip(means, means[1:])),
                n_total=int(sum(r["n"] for r in rows)))
