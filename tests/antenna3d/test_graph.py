"""The deliverable: node kinds, edge attributes, denominators, and the GraphML round trip."""
import numpy as np
import networkx as nx
import pytest

from antenna3d.graph import DENOMINATORS, build_graph, edge_rows, graph_attrs, nucleus_row
from antenna3d.grid import Grid
from antenna3d.params import Params
from antenna3d.reconnect import build_filaments
from antenna3d.segment import nuclei_from_labels

import phantom
from test_trace import run_detect


@pytest.fixture(scope="module")
def built():
    """One phantom nucleus, all the way to a graph."""
    p = Params()
    vol, mask, grid, _, true_um = phantom.one_filament_nucleus(noise_sd=300.0, amplitude=3000.0)
    crop, enh, fg, tr = run_detect(vol, mask, grid, p)

    lgrid = Grid(grid.dz_um, grid.dy_um * 5, grid.dx_um * 5)
    from scipy.ndimage import zoom
    coarse = zoom(mask.astype(np.float32), (1.0, 1 / 5, 1 / 5), order=1) > 0.5
    labels = coarse.astype(np.uint16)
    n = nuclei_from_labels(labels, lgrid, p)[0]
    n.nucleus_uid = "phantom_p00_nuc001"

    fil = build_filaments(tr["skeleton"], tr["grid"], crop.actin_flat, crop.grid, p,
                          noise_sd=crop.noise_sd)
    fil["min_branch_um"] = tr["min_branch_um"]
    return build_graph(n, crop, fil, p), true_um, p


def test_the_graph_has_edges_and_they_carry_a_polyline(built):
    G, true_um, _ = built
    assert G.number_of_edges() >= 1
    for _, _, e in G.edges(data=True):
        pts = [q.split(",") for q in e["polyline"].split(";")]
        assert len(pts) == e["n_points"] and all(len(q) == 3 for q in pts)
        assert e["length"] > 0 and e["chord_length"] > 0


def test_total_length_matches_the_phantom(built):
    G, true_um, _ = built
    assert G.graph["total_antenna_length_um"] == pytest.approx(true_um, rel=0.10)


def test_every_node_has_a_kind_from_the_declared_set(built):
    G, _, _ = built
    kinds = {d["kind"] for _, d in G.nodes(data=True)}
    assert kinds <= {"root", "tip", "branch", "cut", "crossing"}
    for _, d in G.nodes(data=True):
        for k in ("degree", "z", "y", "x", "boundary_distance", "at_cut_face"):
            assert k in d


def test_a_straight_filament_inside_a_closed_nucleus_has_tips_not_cut_ends(built):
    """The phantom nucleus is a closed sphere, so nothing is against a manufactured face."""
    G, _, _ = built
    assert G.graph["n_cut_ends"] == 0
    assert not G.graph["mask_cut_low"] and not G.graph["mask_cut_high"]
    assert G.graph["n_ends_total"] >= 2


def test_cut_ends_are_found_when_the_mask_really_is_cut():
    """Against a cut face an end is neither root nor tip, and calling it either inflates a
    number that is then read as biology."""
    # The phantom nucleus is ~107 um3 whole, so half of one is below the shipped 100 um3 gate.
    # Lowered here only so there is a nucleus to test the cut-face logic on.
    from dataclasses import replace
    p = Params()
    p = p.with_(segment=replace(p.segment, min_volume_um3=20.0))
    vol, mask, grid, _, _ = phantom.one_filament_nucleus(noise_sd=300.0, amplitude=3000.0)
    z0 = mask.any(axis=(1, 2)).argmax() + 6
    vol, mask = vol[z0:], mask[z0:]                     # slice the sphere through its middle
    grid = Grid(grid.dz_um, grid.dy_um, grid.dx_um, oz_um=z0 * grid.dz_um)
    crop, _, fg, tr = run_detect(vol, mask, grid, p)

    lgrid = Grid(grid.dz_um, grid.dy_um * 5, grid.dx_um * 5)
    from scipy.ndimage import zoom
    labels = (zoom(mask.astype(np.float32), (1.0, 1 / 5, 1 / 5), order=1) > 0.5).astype(np.uint16)
    n = nuclei_from_labels(labels, lgrid, p)[0]
    n.nucleus_uid = "cut_p00_nuc001"
    assert n.axial_edge_fraction > 0.3, "the test mask must actually be cut"

    fil = build_filaments(tr["skeleton"], tr["grid"], crop.actin_flat, crop.grid, p,
                          noise_sd=crop.noise_sd)
    G = build_graph(n, crop, fil, p)
    assert G.graph["mask_cut_low"], "the mask reaches its own first plane at full cross-section"
    # Every component/tip observable is reported raw AND cut-excluded, in matched pairs.
    for raw, uncut in (("n_components", "n_components_uncut"),):
        assert raw in G.graph and uncut in G.graph
        assert G.graph[uncut] <= G.graph[raw]


def test_every_density_states_its_denominator(built):
    G, _, _ = built
    a = G.graph
    for k in DENOMINATORS:
        assert k in a and a[k] > 0, k
    assert a["length_density_um_per_um3"] == pytest.approx(
        a["total_antenna_length_um"] / a["nucleus_volume_um3"])
    # The SURFACE density divides by the lid-free envelope, not the marching-cubes total.
    assert a["length_density_um_per_um2_surface"] == pytest.approx(
        a["total_antenna_length_um"] / a["nucleus_surface_um2_true"])
    assert a["nucleus_surface_um2_true"] <= a["nucleus_surface_um2"]


def test_no_total_is_reported_without_a_density_beside_it(built):
    G, _, _ = built
    assert "total_antenna_length_um" in G.graph
    assert "length_density_um_per_um3" in G.graph
    assert "n_antennas" in G.graph and "antennas_per_um3" in G.graph


def test_the_width_is_flagged_as_psf_limited(built):
    G, _, _ = built
    assert G.graph["width_is_psf_limited"] is True
    assert G.graph["psf_fwhm_um"] == pytest.approx(0.2302)
    for _, _, e in G.edges(data=True):
        if np.isfinite(e["mean_width"]):
            assert e["mean_width_deconvolved"] <= e["mean_width"]


def test_units_and_order_are_declared(built):
    G, _, _ = built
    assert G.graph["units"] == "micron" and G.graph["coordinate_order"] == "zyx"


def test_graphml_round_trip(built, tmp_path):
    G, _, _ = built
    p = tmp_path / "n.graphml"
    nx.write_graphml(G, p)
    H = nx.read_graphml(p)
    assert H.number_of_edges() == G.number_of_edges()
    # networkx injects node_default/edge_default on read; pyarrow refuses a struct column.
    assert "node_default" not in graph_attrs(H)
    assert float(graph_attrs(H)["total_antenna_length_um"]) == pytest.approx(
        G.graph["total_antenna_length_um"])


def test_the_tables_stand_alone(built):
    G, _, _ = built
    rows = edge_rows(G)
    assert len(rows) == G.number_of_edges()
    for r in rows:
        assert r["nucleus_uid"] == G.graph["nucleus_uid"]
        for k in DENOMINATORS:
            assert k in r, "a per-edge density must not need a join to get its denominator"
    nr = nucleus_row(G)
    assert nr["n_edges"] == G.number_of_edges() and nr["detected"] is True
