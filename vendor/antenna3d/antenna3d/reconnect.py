"""Crossing resolution and gap joining, by orientation and intensity continuity.

An explicit, testable step with its own parameters - **not** a side effect of morphological
cleanup. Skeletonising a network fuses filaments wherever they cross: a crossing becomes one
junction node and the identity of which arm continues into which is lost. Measured on real
data, reporter-positive nuclei carry 2.1-2.3 branch points per um of centreline, so this is the
common case rather than an edge case.

The decision at a junction is a matching problem over the incident branches. Two branches pair
as one filament passing through if their end tangents are close to anti-parallel AND their mean
intensities agree. Anything unpaired stays a genuine branch of the network.

**What this step is measured NOT to fix.** Widening `max_gap_um` does not reunite broken
filaments: measured, the number of joins that land on the same object is flat at every gap
width while the fraction landing on a DIFFERENT object rises from 83% to 97%. The gaps in a
real trace are mostly regions where nothing was detected at all, and nothing can be joined
across a region with no evidence in it. Leave the default alone.
"""
from __future__ import annotations

import itertools

import numpy as np
from skan import summarize

from .grid import Grid, sample_um
from .trace import safe_skeleton, smooth_polyline


def end_tangent(coords_um: np.ndarray, at_start: bool, window_um: float,
                skip_um: float = 0.0) -> np.ndarray:
    """Unit tangent at one end of a polyline, pointing AWAY from that end and INTO the segment.

    **THE CONVENTION MATTERS AND IS EASY TO INVERT.** At a free end the tangent points back
    along its own filament, i.e. away from any gap beyond it. Two broken halves of one filament
    therefore have tangents that point APART, not toward each other.

    A total-least-squares fit over an arc-length window, not a two-point difference. Two
    reasons, both measured: a one-voxel step on the work grid is quantised to 26 directions, and
    a skeleton junction lands up to ~0.2 um off the true crossing in z, so a short window reads
    the climb out of the junction rather than the filament. On an X phantom a 0.5 um two-point
    window put all four arms 55 degrees out of plane and no pair could be made.
    """
    p = np.asarray(coords_um, float)
    p = p if at_start else p[::-1]
    if len(p) < 2:
        return np.zeros(3)
    d = np.linalg.norm(np.diff(p, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(d)])
    # Skip the first `skip_um` of arc length: a skeleton climbs out of a junction over roughly
    # one work-grid z voxel, and a tangent fitted through that reads the junction.
    j = int(np.searchsorted(s, skip_um)) if skip_um > 0 else 0
    j = min(j, max(len(p) - 3, 0))
    k = int(np.searchsorted(s, s[j] + window_um))
    k = max(j + 1, min(k, len(p) - 1))
    seg = p[j:k + 1]
    anchor = p[0]
    if len(seg) < 3:
        v = seg[-1] - seg[0]
    else:
        c = seg - seg.mean(axis=0)
        _, _, vt = np.linalg.svd(c, full_matrices=False)
        v = vt[0]
        if np.dot(v, seg[-1] - seg[0]) < 0:
            v = -v
    n = np.linalg.norm(v)
    if n == 0:
        return np.zeros(3)
    v = v / n
    if len(seg) and np.dot(v, seg[-1] - anchor) < 0:     # still pointing away from the end
        v = -v
    return v


def turn_angle_deg(t_a: np.ndarray, t_b: np.ndarray) -> float:
    """Turn required for a filament arriving along -t_a to leave along t_b.

    Both tangents point away from the shared junction, so a straight pass-through has
    ``t_a . t_b = -1``, i.e. a turn of 0 degrees.
    """
    c = float(np.clip(-np.dot(t_a, t_b), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def intensity_ratio(ia: float, ib: float, floor: float) -> float:
    """Ratio of two mean intensities, with a noise floor.

    The intensities are background-SUBTRACTED, so they can be near zero or negative on a faint
    filament, and a bare ratio of two small numbers is unstable or undefined. Adding the noise
    sd as a floor makes the comparison "do these agree, given how well either is measured"
    rather than "is one a multiple of the other".
    """
    a, b = ia + floor, ib + floor
    if a <= 0 or b <= 0:
        return float("inf")
    return float(max(a / b, b / a))


def pair_at_junction(tangents, intensities, max_angle_deg, max_intensity_ratio,
                     intensity_floor=0.0):
    """Optimal pairing of the branches incident on one junction.

    Returns a list of ``(i, j)`` pairs. Branches left unpaired are real branch points. The
    search is exhaustive because a junction has 3-6 incident branches in practice.
    """
    n = len(tangents)
    best, best_cost = [], np.inf

    def admissible(i, j):
        a = turn_angle_deg(tangents[i], tangents[j])
        if a > max_angle_deg:
            return None
        r = intensity_ratio(intensities[i], intensities[j], intensity_floor)
        if not np.isfinite(r) or r > max_intensity_ratio:
            return None
        return a / max(max_angle_deg, 1e-9) + (r - 1.0) / max(max_intensity_ratio - 1.0, 1e-9)

    def search(remaining, pairs, cost):
        nonlocal best, best_cost
        if cost >= best_cost:
            return
        if len(remaining) < 2:
            if cost < best_cost:
                best_cost, best = cost, list(pairs)
            return
        i = remaining[0]
        search(remaining[1:], pairs, cost + 1.0)          # leave i unpaired
        for j in remaining[1:]:
            c = admissible(i, j)
            if c is None:
                continue
            search([k for k in remaining[1:] if k != j], pairs + [(i, j)], cost + c)

    search(list(range(n)), [], 0.0)
    return best


def join_gaps(ends, max_gap_um, max_angle_deg, max_intensity_ratio, intensity_floor=0.0):
    """Pair free endpoints across a gap.

    `ends` is a list of dicts with keys `id`, `point` (3,), `tangent` (3,), `intensity`, and
    optionally `group` - the segment the end belongs to.

    A join needs all four of: the two ends on DIFFERENT segments, the gap short enough, the two
    tangents anti-parallel enough, and each tangent pointing roughly ALONG the gap rather than
    across it.

    **The same-segment exclusion is not defensive tidiness.** Without it every segment shorter
    than `max_gap_um` joins to ITSELF: its own two ends are close, its two end tangents point
    into the segment and are therefore anti-parallel, each points along the gap, and its
    intensity ratio against itself is exactly 1. Every test passes. Measured before the fix: all
    19 "gap joins" in one nucleus were self-joins, so the step reported work it had not done.
    """
    out = []
    used: set = set()
    cands = []
    c = np.cos(np.radians(max_angle_deg))
    for a, b in itertools.combinations(range(len(ends)), 2):
        ga, gb = ends[a].get("group"), ends[b].get("group")
        if ga is not None and ga == gb:
            continue
        pa, pb = ends[a]["point"], ends[b]["point"]
        gap = float(np.linalg.norm(pb - pa))
        if gap > max_gap_um or gap == 0:
            continue
        ta, tb = ends[a]["tangent"], ends[b]["tangent"]
        turn = turn_angle_deg(ta, tb)
        if turn > max_angle_deg:
            continue
        u = (pb - pa) / gap
        # Each end tangent points INTO ITS OWN SEGMENT, so for a genuine break the two point
        # APART along the gap: dot(t_a, u) ~ -1 and dot(t_b, u) ~ +1. Testing the opposite sign
        # admits exactly one configuration - a segment joined to its own other end.
        if float(np.dot(ta, -u)) < c or float(np.dot(tb, u)) < c:
            continue
        ratio = intensity_ratio(ends[a]["intensity"], ends[b]["intensity"], intensity_floor)
        if not np.isfinite(ratio) or ratio > max_intensity_ratio:
            continue
        cands.append((turn / max_angle_deg + gap / max_gap_um, a, b, gap, turn, ratio))
    for cost, a, b, gap, turn, ratio in sorted(cands):
        if a in used or b in used:
            continue
        used.update((a, b))
        out.append({"a": ends[a]["id"], "b": ends[b]["id"], "gap_um": gap,
                    "turn_deg": turn, "intensity_ratio": ratio, "cost": cost})
    return out


def build_filaments(sk: np.ndarray, wgrid: Grid, actin: np.ndarray, agrid: Grid, params,
                    noise_sd: float = 0.0) -> dict:
    """Skeleton -> segments, junction decisions, gap joins and the filament chains they imply."""
    S = safe_skeleton(sk, wgrid.spacing)
    if S is None:
        return {"segments": [], "junctions": [], "gap_joins": [], "filaments": []}
    df = summarize(S, separator="_")
    t = params.trace
    s = params.scales()
    win, skip, merge_um = s["tangent_window_um"], s["tangent_skip_um"], s["junction_merge_um"]

    segs = []
    for pid in range(S.n_paths):
        c = S.path_coordinates(pid).astype(float)
        pts = c * np.array(wgrid.spacing) + np.array(wgrid.origin)
        pts = smooth_polyline(pts, t.polyline_smooth_passes)
        vals = sample_um(actin, agrid, pts)
        d = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        mean_i = float(np.mean(vals))
        segs.append({
            "id": int(pid),
            "src": int(df.node_id_src.iloc[pid]), "dst": int(df.node_id_dst.iloc[pid]),
            "branch_type": int(df.branch_type.iloc[pid]),
            "length_um": float(d.sum()),
            "chord_um": float(np.linalg.norm(pts[-1] - pts[0])),
            "mean_intensity": mean_i,
            "intensity_cv": float(np.std(vals) / abs(mean_i)) if mean_i else 0.0,
            "n_points": int(len(pts)),
            "points_um": np.round(pts, 4).tolist(),
            "t_src": end_tangent(pts, True, win, skip).tolist(),
            "t_dst": end_tangent(pts, False, win, skip).tolist(),
        })

    # --- merge split junction clusters -------------------------------------------------------
    # A skeleton junction routinely lands on two or three adjacent voxels, which skan reports as
    # separate nodes joined by a sub-voxel stub. Left alone, one Y becomes two branch points and
    # one X becomes a degree-5 node that cannot be paired at all.
    node_pt: dict[int, np.ndarray] = {}
    for sg in segs:
        node_pt.setdefault(sg["src"], np.asarray(sg["points_um"][0], float))
        node_pt.setdefault(sg["dst"], np.asarray(sg["points_um"][-1], float))
    parent = {n: n for n in node_pt}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    keys = list(node_pt)
    for ia in range(len(keys)):
        for ib in range(ia + 1, len(keys)):
            a, b = keys[ia], keys[ib]
            if np.linalg.norm(node_pt[a] - node_pt[b]) <= merge_um:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[rb] = ra
    if any(find(n) != n for n in node_pt):
        kept = []
        for sg in segs:
            sg["src"], sg["dst"] = find(sg["src"]), find(sg["dst"])
            if sg["src"] == sg["dst"] and sg["length_um"] <= 2 * merge_um:
                continue                        # the sub-voxel stub inside a junction cluster
            kept.append(sg)
        for i, sg in enumerate(kept):
            sg["id"] = i
        segs = kept

    # --- junction decisions ------------------------------------------------------------------
    incident: dict[int, list[tuple[int, str]]] = {}
    for sg in segs:
        incident.setdefault(sg["src"], []).append((sg["id"], "src"))
        incident.setdefault(sg["dst"], []).append((sg["id"], "dst"))
    junctions = []
    pair_of: dict[tuple[int, int], int] = {}
    for node, arms in incident.items():
        if len(arms) < 3:
            continue
        tang = [np.array(segs[i]["t_src" if e == "src" else "t_dst"]) for i, e in arms]
        inten = [segs[i]["mean_intensity"] for i, _ in arms]
        pairs = pair_at_junction(tang, inten, t.max_angle_deg, t.max_intensity_ratio,
                                 intensity_floor=noise_sd) if t.reconnect else []
        resolved = len(pairs) * 2 == len(arms)
        pt = segs[arms[0][0]]["points_um"][0 if arms[0][1] == "src" else -1]
        junctions.append({"node": int(node), "degree": len(arms),
                          "point_um": list(map(float, pt)),
                          "pairs": [[int(arms[a][0]), int(arms[b][0])] for a, b in pairs],
                          "resolved": bool(resolved)})
        if resolved:
            for a, b in pairs:
                pair_of[(node, arms[a][0])] = arms[b][0]
                pair_of[(node, arms[b][0])] = arms[a][0]

    # --- gap joins between free ends ---------------------------------------------------------
    ends = []
    for sg in segs:
        for node, key, tkey in ((sg["src"], "src", "t_src"), (sg["dst"], "dst", "t_dst")):
            if len(incident.get(node, [])) == 1:
                p = sg["points_um"][0 if key == "src" else -1]
                ends.append({"id": (sg["id"], key), "group": sg["id"],
                             "point": np.array(p, float),
                             "tangent": np.array(sg[tkey], float),
                             "intensity": sg["mean_intensity"]})
    gaps = join_gaps(ends, t.max_gap_um, t.max_angle_deg, t.max_intensity_ratio,
                     intensity_floor=noise_sd) if t.reconnect else []
    for g in gaps:
        g["a"] = [int(g["a"][0]), g["a"][1]]
        g["b"] = [int(g["b"][0]), g["b"][1]]

    # --- filament chains: walk pass-through pairings and gap joins ---------------------------
    adj: dict[int, set[int]] = {}
    for (node, a), b in pair_of.items():
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    for g in gaps:
        a, b = g["a"][0], g["b"][0]
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    seen, filaments = set(), []
    for sg in segs:
        if sg["id"] in seen:
            continue
        stack, comp = [sg["id"]], []
        while stack:
            k = stack.pop()
            if k in seen:
                continue
            seen.add(k)
            comp.append(k)
            stack.extend(adj.get(k, ()))
        filaments.append(sorted(comp))
    return {"segments": segs, "junctions": junctions, "gap_joins": gaps, "filaments": filaments}
