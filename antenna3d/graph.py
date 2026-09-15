"""The per-nucleus graph. This is the deliverable; everything else feeds it.

One `networkx` graph per nucleus, written as GraphML, plus the two tidy tables flattened out of
it. Nodes are ends, branch points and roots; edges are traced segments carrying their polyline,
length, tortuosity, intensity and apparent width.

**Every quantity is in microns and every density states its denominator.** A total is not
comparable between nuclei - a nucleus here is 250-900 um3 of imaged volume - so what is
reported is length per unit volume, length per unit envelope area, and counts per unit length.
"""
from __future__ import annotations

import networkx as nx
import numpy as np
from scipy.ndimage import distance_transform_edt

from .grid import Grid, sample_um


def signed_distance(mask: np.ndarray, grid: Grid) -> np.ndarray:
    """Signed distance to the nuclear surface in microns. Positive inside, negative outside.

    **Anisotropic sampling.** A native grid here is (0.2, 0.046, 0.046) um, so an unsampled EDT
    would report a z-step as 4.35x too short and every root classification would be wrong.
    """
    din = distance_transform_edt(mask, sampling=grid.spacing)
    dout = distance_transform_edt(~mask, sampling=grid.spacing)
    return np.where(mask, din, -dout).astype(np.float32)


def surface_normal_field(sdist: np.ndarray, grid: Grid) -> np.ndarray:
    """Unit OUTWARD normal, from the gradient of the signed distance field."""
    g = np.gradient(sdist, *grid.spacing)
    v = -np.stack(g, axis=-1)                       # sdist increases inward
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return np.divide(v, n, out=np.zeros_like(v), where=n > 1e-9)


def boundary_arclength(mask_plane: np.ndarray, grid: Grid, y_um: float, x_um: float):
    """Arc-length position of the nearest boundary point on this z-plane, and the perimeter."""
    from skimage.measure import find_contours
    cs = find_contours(mask_plane.astype(float), 0.5)
    if not cs:
        return float("nan"), float("nan")
    best = None
    for c in cs:
        p = np.column_stack([c[:, 0] * grid.dy_um + grid.oy_um,
                             c[:, 1] * grid.dx_um + grid.ox_um])
        d = np.hypot(p[:, 0] - y_um, p[:, 1] - x_um)
        i = int(np.argmin(d))
        if best is None or d[i] < best[0]:
            seg = np.linalg.norm(np.diff(p, axis=0, append=p[:1]), axis=1)
            best = (float(d[i]), float(np.concatenate([[0.0], np.cumsum(seg)])[i]),
                    float(seg.sum()))
    return best[1], best[2]


def _fwhm(prof, step):
    p = prof - prof.min()
    i = int(np.argmax(p))
    if p[i] <= 0:
        return float("nan")
    h = p[i] / 2.0
    l = i
    while l > 0 and p[l] > h:
        l -= 1
    r = i
    while r < len(p) - 1 and p[r] > h:
        r += 1
    if p[l] > h or p[r] > h:
        return float("nan")
    xl = l + (h - p[l]) / max(p[l + 1] - p[l], 1e-9)
    xr = r - (h - p[r]) / max(p[r - 1] - p[r], 1e-9)
    return (xr - xl) * step


def perpendicular_fwhm(actin, agrid, pts_um, tangents, half_um=0.46, n=41):
    """Apparent cross-section FWHM along a polyline, over two perpendicular directions.

    Measured on the DATA; it is not the filter's scale. On this optics it is PSF-dominated -
    0.251 um measured against a 0.180 um PSF - so the quadrature-deconvolved intrinsic width is
    ~0.175 um and the number should be read as an UPPER BOUND, not a width. The graph carries
    `width_is_psf_limited` to say so.

    Three independent estimators of filament width on the same data disagree by 40% (0.251,
    0.301, 0.31-0.35 um) depending on how profiles are selected. What they agree on is the
    thing worth carrying: filaments are wider than the PSF, everywhere, by 1.35-1.51x.
    """
    out = []
    t = np.linspace(-half_um, half_um, n)
    for p, tan in zip(pts_um, tangents):
        tan = np.asarray(tan, float)
        if np.linalg.norm(tan) < 1e-9:
            continue
        a = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(a, tan)) > 0.9:
            a = np.array([0.0, 1.0, 0.0])
        u = np.cross(tan, a)
        u /= np.linalg.norm(u)
        v = np.cross(tan, u)
        for d in (u, v):
            prof = sample_um(actin, agrid, p[None, :] + t[:, None] * d[None, :])
            w = _fwhm(prof, t[1] - t[0])
            if np.isfinite(w):
                out.append(w)
    return float(np.median(out)) if out else float("nan")


def _poly_str(pts):
    return ";".join(f"{z:.4f},{y:.4f},{x:.4f}" for z, y, x in pts)


def build_graph(nucleus, crop, fil, params, field=None) -> nx.Graph:
    """Assemble the per-nucleus graph from the resolved filament objects."""
    segs = fil["segments"]
    mask, mgrid = crop.mask, crop.grid
    actin, agrid = crop.actin_flat, crop.grid
    gp = params.graph
    psf = params.optics.psf()
    sdist = signed_distance(mask, mgrid)
    normals = surface_normal_field(sdist, mgrid)
    tol = float(gp.root_tolerance_um)
    psf_fwhm = float(psf.lateral_fwhm_um)

    resolved_nodes = {j["node"] for j in fil["junctions"] if j["resolved"]}
    fil_of = {sid: k for k, comp in enumerate(fil["filaments"]) for sid in comp}

    G = nx.Graph()
    node_pts: dict[int, np.ndarray] = {}
    degree: dict[int, int] = {}
    for s in segs:
        pts = np.asarray(s["points_um"], float)
        node_pts.setdefault(s["src"], pts[0])
        node_pts.setdefault(s["dst"], pts[-1])
        degree[s["src"]] = degree.get(s["src"], 0) + 1
        degree[s["dst"]] = degree.get(s["dst"], 0) + 1

    # Where the mask is cut in z, the surface there is MANUFACTURED - the boundary mesh is
    # floored off with a flat lid. `sdist` is distance to that surface, so a filament ending
    # against the lid sits within `root_tolerance_um` of it and would be classified **root**,
    # carrying the exit angle, boundary azimuth and the rest of the sub-envelope observables.
    # It is not a root; it is a cut end, and its exit angle points along z.
    #
    # This is the slab on a thin acquisition AND the segmentation on a thick one. One test
    # covers both: an end is at a cut face when the mask's cross-section THERE is still a
    # substantial fraction of its peak, which is what distinguishes a cut from a nucleus
    # tapering to its pole. Tested against the MASK's own axial extent, never the crop's - the
    # crop carries a 1.5 um margin, so a mask cut by the SEGMENTATION ends inside the crop and
    # a `mask[0].any()` test would miss it entirely.
    areas = mask.reshape(mask.shape[0], -1).sum(axis=1)
    occ = np.flatnonzero(areas)
    if occ.size:
        peak = float(areas[occ].max())
        cap_lo = bool(areas[occ[0]] / peak > gp.cut_area_fraction)
        cap_hi = bool(areas[occ[-1]] / peak > gp.cut_area_fraction)
        cap_z_lo = mgrid.origin[0] + occ[0] * mgrid.dz_um
        cap_z_hi = mgrid.origin[0] + occ[-1] * mgrid.dz_um
    else:
        cap_lo = cap_hi = False
        cap_z_lo = cap_z_hi = float("nan")
    # One axial PSF FWHM: the distance over which a cut face is optically indistinct, so an end
    # nearer than this cannot be told from a filament continuing through it.
    cut_tol = float(psf.axial_fwhm_um)

    def at_cut_face(pz: float) -> bool:
        return ((cap_lo and abs(pz - cap_z_lo) <= cut_tol)
                or (cap_hi and abs(pz - cap_z_hi) <= cut_tol))

    for n, p in node_pts.items():
        if n in resolved_nodes:
            continue                              # a dissolved crossing is not a graph node
        d = float(sample_um(sdist, mgrid, p[None, :])[0])
        deg = degree.get(n, 0)
        # `cut` takes precedence over BOTH root and tip: it is neither, and counting it as
        # either inflates a number that is then read as biology.
        if deg <= 1 and at_cut_face(float(p[0])):
            kind = "cut"
        else:
            kind = "root" if abs(d) <= tol else ("tip" if deg <= 1 else "branch")
        idx = tuple(int(round((p[a] - mgrid.origin[a]) / mgrid.spacing[a])) for a in range(3))
        idx = tuple(min(max(i, 0), s - 1) for i, s in zip(idx, mask.shape))
        nrm = normals[idx]
        attrs = {"kind": kind, "degree": int(deg),
                 "z": float(p[0]), "y": float(p[1]), "x": float(p[2]),
                 "boundary_distance": d,
                 "normal_z": float(nrm[0]), "normal_y": float(nrm[1]), "normal_x": float(nrm[2]),
                 "at_cut_face": bool(at_cut_face(float(p[0])))}
        if kind == "root":
            arclen, perim = boundary_arclength(mask[idx[0]], mgrid, p[1], p[2])
            r = p - np.array(nucleus.centroid_um)
            attrs.update({"boundary_arclength_um": arclen, "boundary_perimeter_um": perim,
                          "boundary_azimuth_rad": float(np.arctan2(r[2], r[1])),
                          "boundary_elevation_rad": float(np.arctan2(
                              r[0], np.hypot(r[1], r[2]) + 1e-12))})
        G.add_node(int(n), **attrs)

    n_edges_kept = 0
    for s in segs:
        for n in (s["src"], s["dst"]):
            if n in resolved_nodes and not G.has_node(int(n)):
                # incident on a dissolved crossing: keep the edge, but anchor it to the
                # crossing point as a degree-2 pass-through node so the polyline is not lost
                p = node_pts[n]
                G.add_node(int(n), kind="crossing", degree=int(degree.get(n, 0)),
                           z=float(p[0]), y=float(p[1]), x=float(p[2]),
                           boundary_distance=float(sample_um(sdist, mgrid, p[None, :])[0]),
                           normal_z=0.0, normal_y=0.0, normal_x=0.0, at_cut_face=False)
        pts = np.asarray(s["points_um"], float)
        tang = np.gradient(pts, axis=0)
        tang /= np.maximum(np.linalg.norm(tang, axis=1, keepdims=True), 1e-12)
        w = perpendicular_fwhm(actin, agrid, pts, tang) if gp.measure_width else float("nan")
        chord = s["chord_um"]
        e = {"length": s["length_um"], "chord_length": chord,
             "tortuosity": float(s["length_um"] / chord) if chord > 1e-9 else float("nan"),
             "mean_intensity": s["mean_intensity"], "intensity_cv": s["intensity_cv"],
             "mean_width": w,
             "mean_width_deconvolved": float(np.sqrt(max(w ** 2 - psf_fwhm ** 2, 0.0)))
             if np.isfinite(w) else float("nan"),
             "n_points": s["n_points"], "filament_id": int(fil_of.get(s["id"], -1)),
             "segment_id": int(s["id"]), "polyline": _poly_str(pts),
             "orientation_src_z": float(tang[0][0]), "orientation_src_y": float(tang[0][1]),
             "orientation_src_x": float(tang[0][2]),
             "orientation_dst_z": float(tang[-1][0]), "orientation_dst_y": float(tang[-1][1]),
             "orientation_dst_x": float(tang[-1][2])}
        # `tvec` is the EXIT direction: the way the filament would continue if extended out
        # through the envelope. For a radial filament that is the outward normal, so a radial
        # antenna has exit_angle 0 and one running along the envelope has 90.
        for n, tvec, key in ((s["src"], -tang[0], "src"), (s["dst"], tang[-1], "dst")):
            if G.nodes.get(int(n), {}).get("kind") == "root":
                nd = G.nodes[int(n)]
                nrm = np.array([nd["normal_z"], nd["normal_y"], nd["normal_x"]])
                if np.linalg.norm(nrm) > 1e-9:
                    e[f"exit_angle_{key}"] = float(
                        np.degrees(np.arccos(np.clip(np.dot(tvec, nrm), -1, 1))))
        G.add_edge(int(s["src"]), int(s["dst"]), **e)
        n_edges_kept += 1

    total_len = float(sum(s["length_um"] for s in segs))
    n_branch = sum(1 for _, d in G.nodes(data=True) if d["kind"] == "branch")
    n_fil = len(fil["filaments"])
    vol = float(nucleus.volume_um3)
    surf_true = float(nucleus.surface_um2_true)
    imaged_vol = float(mask.sum()) * mgrid.voxel_volume_um3

    G.graph.update({
        "nucleus_uid": nucleus.nucleus_uid, "label": int(nucleus.label_id),
        "field_id": getattr(field, "field_id", ""),
        "file": getattr(field, "stem", ""), "position": int(getattr(field, "position", 0)),
        "source_path": str(getattr(field, "path", "")),
        # From the input subfolder. Carried verbatim - nothing here parses it, because a
        # trailing token in a folder name is not automatically a drug.
        "condition": str(getattr(field, "condition", "") or ""),

        # --- denominators. Every density below names the one it used. ---------------------
        "nucleus_volume_um3": vol,
        # Counted on the NATIVE crop grid rather than the coarse label grid. The two measure the
        # same nucleus at two resolutions and differ by a percent or two; a density must not mix
        # them, so both are named.
        "nucleus_imaged_volume_um3": imaged_vol,
        "nucleus_surface_um2": float(nucleus.surface_um2),
        "nucleus_surface_um2_true": surf_true,
        "nucleus_surface_cap_fraction": float(nucleus.cap_fraction),
        "nucleus_equivalent_diameter_um": float(nucleus.equivalent_diameter_um),
        "axial_edge_fraction": float(nucleus.axial_edge_fraction),
        "nucleus_centroid_z_um": float(nucleus.centroid_um[0]),
        "nucleus_centroid_y_um": float(nucleus.centroid_um[1]),
        "nucleus_centroid_x_um": float(nucleus.centroid_um[2]),

        "n_antennas": int(n_fil), "n_segments": int(n_edges_kept),
        "total_antenna_length_um": total_len,
        "length_density_um_per_um3": total_len / max(vol, 1e-9),
        # Divided by the envelope with the manufactured lid REMOVED. The pipeline this was
        # ported from divided by the lid-inflated surface, which on a cut nucleus is up to 42.5%
        # too large; `nucleus_surface_um2` is still emitted so the two can be compared.
        "length_density_um_per_um2_surface": total_len / max(surf_true, 1e-9),
        "antennas_per_um3": n_fil / max(vol, 1e-9),

        "n_branch_points": int(n_branch),
        "branching_ratio": float(n_branch / n_fil) if n_fil else 0.0,
        "branch_points_per_um": float(n_branch / total_len) if total_len > 0 else 0.0,
        "n_crossings_resolved": int(len(resolved_nodes)),
        "n_gap_joins": int(len(fil["gap_joins"])),

        # Raw and cut-excluded, IN MATCHED PAIRS. Every component- and tip-derived observable is
        # reported both ways, because the difference between them is the artefact and it is not
        # separable from the biology by inspection.
        "n_roots": int(sum(1 for _, d in G.nodes(data=True) if d["kind"] == "root")),
        "n_tips": int(sum(1 for _, d in G.nodes(data=True) if d["kind"] == "tip")),
        "n_cut_ends": int(sum(1 for _, d in G.nodes(data=True) if d["kind"] == "cut")),
        "n_ends_total": int(sum(1 for _, d in G.nodes(data=True) if d["degree"] <= 1)),
        # Components of THIS graph, raw and cut-excluded. Deliberately NOT called
        # `n_antennas_uncut`: `n_antennas` is the filament count assembled by junction pairing
        # and gap joining, and a connected component of G is not the same object. A ratio
        # between the two would look meaningful and mean nothing.
        "n_components": int(sum(1 for c in nx.connected_components(G)
                                if any(G.degree(m) > 0 for m in c))),
        "n_components_uncut": int(sum(
            1 for c in nx.connected_components(G)
            if not any(G.nodes[m].get("kind") == "cut" for m in c)
            and any(G.degree(m) > 0 for m in c))),
        "length_in_uncut_components_um": float(sum(
            sum(G[u][v]["length"] for u, v in G.subgraph(c).edges())
            for c in nx.connected_components(G)
            if not any(G.nodes[m].get("kind") == "cut" for m in c))),

        "mask_cut_low": cap_lo, "mask_cut_high": cap_hi, "cut_face_tolerance_um": cut_tol,
        "min_branch_um": float(fil.get("min_branch_um", params.trace.min_branch_um)),
        "width_is_psf_limited": True, "psf_fwhm_um": psf_fwhm,
        "units": "micron", "coordinate_order": "zyx",
        "scales_digest": params.optics.digest(),
    })
    return G


# --------------------------------------------------------------------------------- tidy tables

#: Graph attributes that identify the nucleus; repeated on every edge row so the edge table
#: stands alone.
ID_KEYS = ("nucleus_uid", "file", "position", "field_id", "condition", "label", "source_path")

#: The four denominators carried onto every edge row, so a per-edge density never has to be
#: joined back to get one.
DENOMINATORS = ("nucleus_volume_um3", "nucleus_imaged_volume_um3",
                "nucleus_surface_um2_true", "nucleus_equivalent_diameter_um")


def graph_attrs(G) -> dict:
    """Scalar graph attributes only.

    networkx injects `node_default` and `edge_default` dicts into `G.graph` on a GraphML read,
    and pyarrow refuses a column of structs. Dropping non-scalars here is what lets the same
    code path serve a freshly built graph and one read back from disk.
    """
    return {k: v for k, v in G.graph.items() if not isinstance(v, (dict, list, set))}


def edge_rows(G) -> list[dict]:
    """One row per edge, carrying its nucleus's ids and denominators."""
    a = graph_attrs(G)
    ids = {k: a.get(k) for k in ID_KEYS}
    dens = {k: a.get(k) for k in DENOMINATORS}
    rows = []
    for u, v, d in G.edges(data=True):
        rows.append({**ids, "u": int(u), "v": int(v),
                     "kind_u": G.nodes[u].get("kind"), "kind_v": G.nodes[v].get("kind"),
                     "boundary_distance_u": G.nodes[u].get("boundary_distance"),
                     "boundary_distance_v": G.nodes[v].get("boundary_distance"),
                     **d, **dens})
    return rows


def nucleus_row(G) -> dict:
    """One row per nucleus: every scalar graph attribute, plus the graph's size."""
    a = graph_attrs(G)
    a["n_nodes"] = int(G.number_of_nodes())
    a["n_edges"] = int(G.number_of_edges())
    a["detected"] = bool(a.get("total_antenna_length_um", 0.0) > 0)
    return a
