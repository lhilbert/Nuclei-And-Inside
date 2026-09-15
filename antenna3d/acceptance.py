"""The probe-absent acceptance test.

**What it is for.** If you have a control condition with no probe in it - not a vehicle control,
an actual *absent reporter* - then the correct length density there is **zero**, and that is the
only test in this package that can tell a working detector from a confidently broken one. Every
other number it produces is a measurement of something; this one is a measurement of the
measurement.

**How to read it.** `detect_neg` is the fraction of probe-absent nuclei in which ANY centreline
was traced. It should be 0. `ratio_med` and the rank-biserial say how well the two conditions
separate if it is not: rank-biserial is +1 when every probe-present nucleus exceeds every
control and 0 when the two are indistinguishable. A ratio of 2x with a rank-biserial of +0.18
means the medians differ and the distributions almost entirely overlap.

**It is reported as a CURVE over `min_branch_um`, never at one value** - and it re-traces from
the stored binary at each value rather than filtering edges afterwards, because pruning a spur
changes the topology. A shorter minimum branch is not the same graph with more edges kept.

**What it found on the data this package was built from.** It FAILS. 66 of 66 probe-absent
nuclei reported antennas; median density ratio 1.56, rank-biserial +0.178, p 0.024 at
`min_branch_um = 0.5`. See the README's "Things that will bite you" for the diagnosis. Run it on
your own control before quoting anything.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .reconnect import build_filaments
from .trace import trace

DEFAULT_GRID = (0.25, 0.5, 1.0, 2.0)


def mannwhitney(a, b):
    """(rank-biserial, two-sided p). NaN below n=3 rather than a number that means nothing."""
    from scipy import stats
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or len(b) < 3:
        return float("nan"), float("nan")
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(2 * u / (len(a) * len(b)) - 1), float(p)


def per_nucleus(work_dir, nuclei_df, params, grid_um=DEFAULT_GRID, log=print) -> pd.DataFrame:
    """Re-trace every stored nucleus at each `min_branch_um` and record its length density.

    Reads the small artefacts `run_folder(keep_work="small")` wrote. Where the flattened actin
    was kept as well (`keep_work="full"`), the branch-point columns are filled in too; without
    it they are NaN, because resolving a crossing needs the intensity along each arm.
    """
    work = Path(work_dir)
    by_uid = nuclei_df.set_index("nucleus_uid")
    rows = []
    for i, (uid, nr) in enumerate(by_uid.iterrows(), 1):
        p = work / f"{uid}.npz"
        if not p.exists():
            continue
        from .pipeline import load_work
        w = load_work(p)
        has_actin = "actin_flat" in w
        for mb in grid_um:
            tr = trace(w["fg"], w["grid"], params, min_branch_um=mb)
            # The SMOOTHED length, which needs no image - so the density reported here is the
            # same quantity whether or not `keep_work="full"` kept the flattened actin.
            length = float(tr["length_um"])
            n_branch = n_fil = float("nan")
            if has_actin:
                fil = build_filaments(tr["skeleton"], tr["grid"], w["actin_flat"], w["grid"],
                                      params, noise_sd=w["noise_sd"])
                length = float(sum(s["length_um"] for s in fil["segments"]))
                n_fil = len(fil["filaments"])
                deg: dict[int, int] = {}
                for s in fil["segments"]:
                    deg[s["src"]] = deg.get(s["src"], 0) + 1
                    deg[s["dst"]] = deg.get(s["dst"], 0) + 1
                resolved = {j["node"] for j in fil["junctions"] if j["resolved"]}
                n_branch = sum(1 for n, d in deg.items() if d >= 2 and n not in resolved)
            vol = float(nr["nucleus_volume_um3"])
            rows.append({
                "nucleus_uid": uid, "field_id": nr.get("field_id"),
                "condition": nr.get("condition"), "min_branch_um": mb,
                "total_length_um": length,
                "length_density_um_per_um3": length / max(vol, 1e-9),
                "n_antennas": n_fil, "n_branch_points": n_branch,
                "nucleus_volume_um3": vol, "detected": bool(length > 0)})
        if i % 20 == 0:
            log(f"  {i}/{len(by_uid)} nuclei re-traced")
    return pd.DataFrame(rows)


def summarise(per_nuc: pd.DataFrame, positive: str, negative: str) -> pd.DataFrame:
    """The curve: one row per `min_branch_um`."""
    out = []
    for mb, g in per_nuc.groupby("min_branch_um"):
        pos = g[g.condition == positive]
        neg = g[g.condition == negative]
        rb, p = mannwhitney(pos.length_density_um_per_um3, neg.length_density_um_per_um3)
        mp = float(pos.length_density_um_per_um3.median()) if len(pos) else float("nan")
        mn = float(neg.length_density_um_per_um3.median()) if len(neg) else float("nan")
        out.append({
            "min_branch_um": mb, "n_pos": len(pos), "n_neg": len(neg),
            "detect_pos": float(pos.detected.mean()) if len(pos) else float("nan"),
            # THE VERDICT LIVES HERE. It should be 0.0.
            "detect_neg": float(neg.detected.mean()) if len(neg) else float("nan"),
            "density_pos_med": mp, "density_neg_med": mn,
            "density_pos_p90": float(pos.length_density_um_per_um3.quantile(0.9)) if len(pos)
            else float("nan"),
            "density_neg_p90": float(neg.length_density_um_per_um3.quantile(0.9)) if len(neg)
            else float("nan"),
            "ratio_med": mp / mn if mn and np.isfinite(mn) and mn > 0 else float("nan"),
            "rank_biserial": rb, "p_value": p})
    return pd.DataFrame(out).sort_values("min_branch_um")


def run_acceptance(output_dir, positive: str, negative: str, params=None,
                   grid_um=DEFAULT_GRID, log=print) -> dict:
    """The whole test over an output tree that `run_folder` produced.

    `positive` and `negative` are values of the `condition` column - i.e. the input subfolder
    names. It REFUSES to run with no probe-absent nuclei rather than reporting a test it did not
    perform.
    """
    from .params import Params
    out = Path(output_dir).expanduser()
    params = params or Params()
    nuclei = pd.read_csv(out / "antenna_nuclei.csv")
    work = out / "work"
    if not work.exists():
        raise FileNotFoundError(
            f"{work} does not exist. The acceptance test re-traces from the stored binary, so "
            f"run_folder must have been called with keep_work='small' (the default) or 'full'.")

    have = set(nuclei.condition.dropna().unique())
    for name, role in ((positive, "positive"), (negative, "negative")):
        if name not in have:
            raise ValueError(f"condition {name!r} ({role}) is not in antenna_nuclei.csv; "
                             f"the conditions present are {sorted(have)}")
    n_neg = int((nuclei.condition == negative).sum())
    if n_neg == 0:
        raise ValueError(
            f"no nuclei in the probe-absent condition {negative!r}. This test is a "
            f"probe-present / probe-absent contrast and there is nothing to contrast against; "
            f"refusing to report a verdict it cannot support.")
    log(f"acceptance: {int((nuclei.condition == positive).sum())} {positive} / {n_neg} "
        f"{negative} nuclei, re-traced at {len(grid_um)} values of min_branch_um")

    pn = per_nucleus(work, nuclei, params, grid_um, log=log)
    s = summarise(pn, positive, negative)
    pn.to_csv(out / "acceptance_per_nucleus.csv", index=False)
    s.to_csv(out / "acceptance.csv", index=False)

    worst = float(s.detect_neg.max())
    verdict = "PASSES" if worst == 0.0 else "FAILS"
    log("\n" + s.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    log(f"\nVERDICT: {verdict}")
    if verdict == "FAILS":
        r = s.loc[s.min_branch_um.idxmax()]
        log(f"  Up to {100 * worst:.0f}% of {negative} nuclei report antennas where the correct "
            f"answer is zero.")
        log(f"  Even at min_branch_um = {r.min_branch_um:g} it is {100 * r.detect_neg:.0f}%.")
        log(f"  No per-nucleus number in this output tree is quotable as biology.")
        log(f"  Look at qc/traces/*{negative}* - that is where the cause is visible.")
    return {"curve": s, "per_nucleus": pn, "verdict": verdict}
