"""Nuclei: bring your own labels, or let cellpose find them.

**Bringing your own is the intended path.** Nucleus segmentation is a solved, separately
maintained problem - if you already have label volumes from `nucleus3d` or anything else, pass
them and skip this module's dependency on cellpose and torch entirely. What this module always
does is the part that is specific to measuring things INSIDE a nucleus: the boundary mesh with
its manufactured axial lid separated out, the truncation covariate, and the size gates.

Two measured facts govern what a nucleus record carries.

**The marching-cubes surface includes a lid that is not envelope.** Closing the mesh requires
padding, and where a mask reaches its own first or last plane the pad floors the object off
with a flat cap. On a cut nucleus that cap was measured at 42.5% of the reported surface, and
0.00 um2 on an untruncated one. `surface_um2` keeps it (for comparison with other tools);
`surface_um2_true` removes it, and **`surface_um2_true` is the one a density divides by.**

**A mask that stops mid-nucleus is common and does not look wrong.** `axial_edge_fraction` is
0 when the mask tapers away to nothing and 1 when it is cut at its widest. Measured on data
that is NOT axially truncated by its acquisition, 2D-per-plane segmentation gave a median of
0.8764 - pure segmentation defect - against 0.0548 for a 3D flow field on the same fields. The
covariate is what let that be seen at all, so it is carried into every table.
"""
from __future__ import annotations

from dataclasses import dataclass, field as _field
from pathlib import Path

import numpy as np
from scipy.ndimage import find_objects

from .grid import Grid


@dataclass
class Nucleus:
    """One nucleus: where it is, how big, how cut, and the mask itself."""
    label_id: int
    bbox_vox: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]
    mask: np.ndarray                      # bool, on the label grid, cropped to bbox
    grid: Grid                            # the label grid (coarse), with the bbox origin
    centroid_um: tuple[float, float, float]
    volume_um3: float
    equivalent_diameter_um: float
    surface_um2: float                    # marching cubes, INCLUDING any manufactured lid
    surface_um2_true: float               # lid removed. THE denominator for a surface density.
    surface_um2_cap: float
    cap_fraction: float
    axial_edge_fraction: float
    axial_peak_area_um2: float
    axial_extent_um: float
    axial_peak_at_edge: bool
    n_planes_spanned: int
    xy_border_touching: bool
    z_truncated: bool
    used: bool
    #: Set by the pipeline once the field is known: `<stem>_p<NN>_nuc<NNN>`.
    nucleus_uid: str = ""
    #: The envelope mesh, in microns. EMPTY unless `nuclei_from_labels(..., keep_mesh=True)`:
    #: `surface_um2_true` above needs marching cubes, but the vertex and face lists themselves
    #: are read by nothing in this package and are ~70 kB per nucleus. Ask for them if you want
    #: to render an envelope.
    vertices_um: list = _field(default_factory=list, repr=False)
    faces: list = _field(default_factory=list, repr=False)

    def to_row(self) -> dict:
        d = {k: v for k, v in self.__dict__.items()
             if k not in ("mask", "grid", "bbox_vox", "vertices_um", "faces", "centroid_um")}
        d["centroid_z_um"], d["centroid_y_um"], d["centroid_x_um"] = self.centroid_um
        return d


def axial_completeness(sub_mask: np.ndarray, grid: Grid) -> dict:
    """How cut is this nucleus in z? **Model-free** - no shape is assumed and none is fitted.

    A closed object's cross-section tapers to nothing at its poles. A cut one stops while its
    cross-section is still large. So the area at the mask's first and last occupied plane,
    relative to its peak, measures truncation directly:

    * 0.0 - the mask tapers away to nothing; nothing was cut off;
    * 1.0 - the mask is cut at its widest point; at least half the nucleus is missing.

    It cannot say how much VOLUME is missing - that needs a shape model, and a parabolic
    spheroid fit was measured failing at a median 20% relative error precisely because the
    masks it fits are themselves cut. It says WHETHER volume is missing and orders nuclei by
    how much, which is what a covariate has to do.
    """
    per_plane = sub_mask.reshape(sub_mask.shape[0], -1).sum(axis=1)
    A = (per_plane * grid.dy_um * grid.dx_um).astype(float)
    occ = np.flatnonzero(A > 0)
    if occ.size == 0:
        return {"axial_edge_fraction": float("nan"), "axial_peak_area_um2": float("nan"),
                "axial_extent_um": 0.0, "axial_peak_at_edge": False}
    A = A[occ[0]:occ[-1] + 1]
    peak = float(A.max())
    return {"axial_edge_fraction": float(max(A[0], A[-1]) / peak) if peak > 0 else float("nan"),
            "axial_peak_area_um2": peak,
            "axial_extent_um": float(len(A) * grid.dz_um),
            "axial_peak_at_edge": bool(int(np.argmax(A)) in (0, len(A) - 1))}


def boundary_surface(mask: np.ndarray, grid: Grid, step: int = 2,
                     cut_area_fraction: float = 0.30) -> dict:
    """Marching-cubes surface, with the MANUFACTURED axial caps separated out.

    The pad is what makes the mesh closed, and it is also what invents a flat lid: where the
    mask is cut in z, the pad puts a zero plane against it and marching cubes floors the object
    off there. That lid is not nuclear envelope, and counting it inflates a DENOMINATOR.

    A lid is found by the mesh's OWN extreme z, not by an assumed plane: `step_size` subsamples
    the volume before marching, so the isosurface does not sit half a voxel outside the
    boundary as it does at step 1. Testing a fixed plane instead made the whole correction a
    silent no-op at step 2 - cap 0.00 um2 where step 1 gives 69.88 against an analytic 75.40.

    A face is a lid when it lies at that extreme, faces along z, AND the mask's cross-section
    there is still above `cut_area_fraction` of its peak.

    **That last test is not in the pipeline this was ported from, and it is needed.** `mask` here
    is cropped to the nucleus's bounding box, so its first and last plane ALWAYS have voxels and
    a bare `mask[0].any()` guard is always true. The normal test alone then still counts the one
    flat facet at a closed pole: measured, 7.19 um2 or 2.6% of the surface of an ellipsoid that
    is not cut at all. The cross-section test is the same criterion `graph.at_cut_face` uses,
    and it distinguishes what it is supposed to - 0.00 um2 on a closed nucleus against 18.5% on
    a cut one.
    """
    from skimage.measure import marching_cubes, mesh_surface_area
    empty = {"n_vertices": 0, "n_faces": 0, "surface_um2": 0.0, "surface_um2_true": 0.0,
             "surface_um2_cap": 0.0, "cap_fraction": 0.0, "vertices": [], "faces": []}
    pad = np.pad(mask, 1)
    try:
        verts, faces, _, _ = marching_cubes(pad.astype(np.float32), 0.5,
                                            spacing=grid.spacing, step_size=int(step))
    except (RuntimeError, ValueError):
        return empty
    verts = verts - np.array(grid.spacing)          # undo the pad
    total = float(mesh_surface_area(verts, faces))

    tri = verts[faces]
    zs = tri[:, :, 0]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1)
    axial = np.divide(np.abs(n[:, 0]), ln, out=np.zeros(len(n)), where=ln > 0) > 0.9
    tol = max(0.75 * int(step) * grid.dz_um, 0.75 * grid.dz_um)
    cap = np.zeros(len(faces), bool)
    areas = mask.reshape(mask.shape[0], -1).sum(axis=1)
    occ = np.flatnonzero(areas)
    peak = float(areas[occ].max()) if occ.size else 0.0
    if peak > 0 and areas[occ[0]] / peak > cut_area_fraction:
        cap |= np.all(np.abs(zs - verts[:, 0].min()) <= tol, axis=1) & axial
    if peak > 0 and areas[occ[-1]] / peak > cut_area_fraction:
        cap |= np.all(np.abs(zs - verts[:, 0].max()) <= tol, axis=1) & axial
    cap_um2 = float(mesh_surface_area(verts, faces[cap])) if cap.any() else 0.0
    return {"n_vertices": int(len(verts)), "n_faces": int(len(faces)),
            "surface_um2": total, "surface_um2_true": total - cap_um2,
            "surface_um2_cap": cap_um2, "cap_fraction": (cap_um2 / total) if total > 0 else 0.0,
            "vertices": np.round(verts, 4).tolist(), "faces": faces.tolist()}


def nuclei_from_labels(labels: np.ndarray, grid: Grid, params,
                       keep_mesh: bool = False) -> list[Nucleus]:
    """Turn a label volume into `Nucleus` records, applying the size and border gates.

    `grid` is the spacing of `labels`, in microns, as (z, y, x). **Get it right.** A label
    volume from another tool is usually at the spacing that tool worked at, and assuming native
    sampling for a binned label volume scales every volume by the bin factor cubed.

    `keep_mesh=True` also returns each envelope's vertices and faces; see `Nucleus`.
    """
    seg = params.segment
    out: list[Nucleus] = []
    Z, Y, X = labels.shape
    vv = grid.voxel_volume_um3
    for li, sl in enumerate(find_objects(labels.astype(np.int32)), start=1):
        if sl is None:
            continue
        sub = labels[sl] == li
        vol_um3 = float(sub.sum()) * vv
        if not (seg.min_volume_um3 <= vol_um3 <= seg.max_volume_um3):
            continue
        zz, yy, xx = np.nonzero(sub)
        z0, y0, x0 = sl[0].start, sl[1].start, sl[2].start
        centroid = ((zz.mean() + z0) * grid.dz_um, (yy.mean() + y0) * grid.dy_um,
                    (xx.mean() + x0) * grid.dx_um)
        border = bool(sl[1].start == 0 or sl[2].start == 0 or sl[1].stop >= Y or sl[2].stop >= X)
        ztrunc = bool(sl[0].start == 0 or sl[0].stop >= Z)
        surf = boundary_surface(sub, grid, seg.boundary_surface_step,
                                params.graph.cut_area_fraction)
        axial = axial_completeness(sub, grid)
        out.append(Nucleus(
            label_id=int(li), bbox_vox=tuple((s.start, s.stop) for s in sl), mask=sub,
            grid=Grid(grid.dz_um, grid.dy_um, grid.dx_um,
                      z0 * grid.dz_um, y0 * grid.dy_um, x0 * grid.dx_um),
            centroid_um=centroid, volume_um3=vol_um3,
            equivalent_diameter_um=float((6 * vol_um3 / np.pi) ** (1 / 3)),
            surface_um2=surf["surface_um2"], surface_um2_true=surf["surface_um2_true"],
            surface_um2_cap=surf["surface_um2_cap"], cap_fraction=surf["cap_fraction"],
            n_planes_spanned=int(sl[0].stop - sl[0].start),
            xy_border_touching=border, z_truncated=ztrunc,
            used=bool(not (border and seg.drop_xy_border_touching)
                      and not (ztrunc and seg.drop_z_truncated)),
            vertices_um=surf["vertices"] if keep_mesh else [],
            faces=surf["faces"] if keep_mesh else [], **axial))
    return out


def labels_from_tiff(path: str | Path, voxel_um=None) -> tuple[np.ndarray, Grid]:
    """A (Z, Y, X) label volume from a TIFF, with its spacing.

    The spacing is read from the OME metadata when present. **If it is not, you must pass
    `voxel_um`** - a label volume with no spacing is not a measurement, and guessing native
    sampling for a binned volume scales every reported volume by the bin factor cubed.
    """
    import tifffile
    with tifffile.TiffFile(str(path)) as tf:
        labels = tf.asarray()
        dz = dy = dx = None
        if voxel_um is None and tf.ome_metadata:
            import re
            m = tf.ome_metadata
            def _g(k):
                r = re.search(rf'PhysicalSize{k}="([0-9.eE+-]+)"', m)
                return float(r.group(1)) if r else None
            dz, dy, dx = _g("Z"), _g("Y"), _g("X")
    if voxel_um is not None:
        dz, dy, dx = voxel_um
    if None in (dz, dy, dx):
        raise ValueError(
            f"{path}: no physical voxel size in the file and none passed. Pass "
            f"voxel_um=(dz, dy, dx) in MICRONS - a label volume without its spacing cannot be "
            f"measured, and assuming one silently rescales every volume it produces.")
    if labels.ndim != 3:
        raise ValueError(f"{path}: expected a 3D (Z, Y, X) label volume, got shape "
                         f"{labels.shape}")
    return labels, Grid(float(dz), float(dy), float(dx))


def normalise_for_cellpose(vol: np.ndarray, percentiles=(1.0, 99.0)) -> np.ndarray:
    """Percentile-normalise the WHOLE volume, once.

    Cellpose's own `normalize=True` runs with `norm3D=True`, so it normalises whatever array it
    is handed. Hand it anything less than the whole field and the same plane enters the network
    as a different image depending on what it was handed with. Measured, when this was done per
    z-chunk: of 512 nuclei, **103 masks ended exactly at the last plane of the first chunk and
    85 began exactly at the next one.**

    Matches `cellpose.transforms.normalize99` term for term, including its guard for a
    collapsed range.
    """
    x = vol.astype(np.float32, copy=True)
    lo, hi = np.percentile(x, percentiles)
    if hi - lo > 1e-3:
        x -= lo
        x /= (hi - lo)
    else:
        x[:] = 0.0
    return x


def block_reduce_mean(vol: np.ndarray, b: int) -> np.ndarray:
    """Bin by `b` in y and x. A 10 um object does not need 0.046 um sampling."""
    if b == 1:
        return vol.astype(np.float32)
    Z, Y, X = vol.shape
    return (vol[:, : Y // b * b, : X // b * b]
            .reshape(Z, Y // b, b, X // b, b).mean(axis=(2, 4)).astype(np.float32))


def label_grid_for(native: Grid, params) -> tuple[Grid, int]:
    """(label grid, bin factor) for segmenting at `SegParams.nuclei_grid_um`."""
    b = int(round(params.segment.nuclei_grid_um / native.dx_um))
    return Grid(dz_um=native.dz_um, dy_um=native.dy_um * b, dx_um=native.dx_um * b), b


def bin_planes(planes, b: int, n_planes: int) -> np.ndarray:
    """Bin an iterable of (Y, X) planes by `b` in y and x, one plane at a time.

    Plane by plane rather than volume-then-bin because the unbinned volume is the largest array
    in the whole pipeline: a 169 x 2280 x 2588 uint16 field is 2.0 GB, against 160 MB binned.
    Nothing needs the full-resolution DNA volume, so nothing should hold one.
    """
    out = None
    for i, pl in enumerate(planes):
        if i >= n_planes:
            break
        row = block_reduce_mean(np.asarray(pl)[None, ...], b)[0]
        if out is None:
            out = np.empty((n_planes,) + row.shape, np.float32)
        out[i] = row
    return out if out is not None else np.zeros((0, 0, 0), np.float32)


def segment_binned(vol_binned: np.ndarray, lgrid: Grid, params, log=print) -> np.ndarray:
    """Cellpose on an already-binned DNA volume. Returns a label volume.

    **`anisotropy` is required and is computed from the grid.** Without it cellpose treats a
    0.2 um z step as a 0.23 um lateral one, computes the flow field on the wrong geometry, and
    is wrong in a way that does not fail.
    """
    try:
        from cellpose import models
    except ImportError as e:                                       # pragma: no cover
        raise ImportError(
            "cellpose is not installed. Either `pip install antenna3d[cellpose]`, or - better - "
            "pass labels you already have to nuclei_from_labels().") from e
    seg = params.segment
    vol = normalise_for_cellpose(vol_binned)
    log(f"      {vol.shape} at {lgrid.dx_um:.4f} um; segmenting")
    model = models.CellposeModel(gpu=bool(seg.gpu))
    extra = dict(do_3D=True, anisotropy=lgrid.dz_um / lgrid.dy_um) if seg.do_3d else {}
    masks, _, _ = model.eval(
        vol, z_axis=0, diameter=seg.diameter_um / lgrid.dy_um,
        flow_threshold=float(seg.flow_threshold),
        cellprob_threshold=float(seg.cellprob_threshold),
        batch_size=8, normalize=False,      # done once above, over the whole volume
        **extra)
    return np.asarray(masks).astype(np.uint16)


def segment_nuclei_cellpose(dna: np.ndarray, native: Grid, params, log=print
                            ) -> tuple[np.ndarray, Grid]:
    """Cellpose on a full-resolution DNA volume. Returns (labels, label grid).

    Optional dependency: `pip install antenna3d[cellpose]`. If you already have masks, use
    `nuclei_from_labels` and none of this runs.
    """
    lgrid, b = label_grid_for(native, params)
    return segment_binned(block_reduce_mean(dna, b), lgrid, params, log), lgrid
