"""The batch driver: a folder of `.nd2` in, an output tree out.

One function, `run_folder`, walks files -> fields -> nuclei and writes everything documented in
the README. It holds no parameters of its own; everything comes from `Params`.

**A field that errors is recorded and the run continues.** One unreadable file, one field whose
segmentation returns nothing, must not cost the other twelve. `failures.csv` is where they go,
and an empty `failures.csv` is the thing to check before reading any table.
"""
from __future__ import annotations

import json
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from . import validate as _validate
from .enhance import enhance
from .graph import build_graph, edge_rows, nucleus_row
from .io import FieldRef, list_fields, nucleus_uid, open_field
from .params import Params
from .preprocess import PlaneCache, Crop, crop_bounds, make_crop
from .reconnect import build_filaments
from .segment import labels_from_tiff, nuclei_from_labels, segment_nuclei_cellpose
from .trace import binarize, trace

WORK_SMALL, WORK_FULL, WORK_NONE = "small", "full", False


def _write_table(df: pd.DataFrame, out: Path, stem: str) -> None:
    """CSV always; parquet as well when pyarrow is available.

    CSV is the one every collaborator can open. Parquet keeps dtypes and is what a large edge
    table should be read from - `polyline` alone is most of the bytes.
    """
    df.to_csv(out / f"{stem}.csv", index=False)
    try:
        df.to_parquet(out / f"{stem}.parquet", index=False)
    except Exception:
        pass


def _save_work(path: Path, crop: Crop, fg: np.ndarray, null_level: float, keep: str) -> None:
    """The small per-nucleus artefact `run_acceptance.py` re-traces from.

    Bit-packed, so a nucleus costs ~30 kB rather than ~9 MB. `keep="full"` adds the flattened
    actin as int16, which is what the branch-point columns of the acceptance table need - and
    which is ~300x larger.
    """
    d = dict(fg=np.packbits(fg), mask=np.packbits(crop.mask),
             shape=np.array(fg.shape), grid=json.dumps(crop.grid.to_dict()),
             noise_sd=np.float32(crop.noise_sd), null_level=np.float32(null_level),
             label_id=int(crop.label_id), nucleus_uid=crop.nucleus_uid)
    if keep == WORK_FULL:
        d["actin_flat"] = np.clip(np.round(crop.actin_flat), -32768, 32767).astype(np.int16)
    np.savez_compressed(path, **d)


def load_work(path: Path):
    """Read back what `_save_work` wrote. Used by `run_acceptance.py`."""
    from .grid import Grid
    z = np.load(path)
    shape = tuple(int(v) for v in z["shape"])
    n = int(np.prod(shape))
    out = {"fg": np.unpackbits(z["fg"], count=n).astype(bool).reshape(shape),
           "mask": np.unpackbits(z["mask"], count=n).astype(bool).reshape(shape),
           "grid": Grid.from_dict(json.loads(str(z["grid"]))),
           "noise_sd": float(z["noise_sd"]), "null_level": float(z["null_level"]),
           "label_id": int(z["label_id"]), "nucleus_uid": str(z["nucleus_uid"])}
    if "actin_flat" in z.files:
        out["actin_flat"] = z["actin_flat"].astype(np.float32)
    return out


def _labels_for(ref: FieldRef, labels_dir, params, field, log):
    """Nuclear labels for one field: yours if you have them, cellpose otherwise."""
    if labels_dir is not None:
        d = Path(labels_dir).expanduser()
        cands = [d / f"{ref.field_id}.tif", d / f"{ref.field_id}.tiff",
                 d / f"{ref.stem}.tif", d / f"{ref.stem}.tiff"]
        hit = next((c for c in cands if c.exists()), None)
        if hit is None:
            raise FileNotFoundError(
                f"no label volume for {ref.field_id} in {d}. Looked for: "
                + ", ".join(c.name for c in cands))
        log(f"      labels from {hit.name}")
        return labels_from_tiff(hit)
    log("      no LABELS_DIR: segmenting with cellpose")
    with_dna = field.volume("dna")
    return segment_nuclei_cellpose(with_dna, field.grid, params, log=log)


def run_field(ref: FieldRef, params: Params, out: Path, labels_dir=None,
              keep_work: str = WORK_SMALL, qc: bool = True, limit: int | None = None,
              multi_series: bool = False, log=print) -> dict:
    """One field, end to end. Returns the QC summary row; writes graphs and work artefacts."""
    t0 = time.time()
    with open_field(ref, params.optics, params.channels, multi_series=multi_series) as field:
        labels, lgrid = _labels_for(ref, labels_dir, params, field, log)
        if labels.shape[0] != field.shape[0]:
            raise ValueError(
                f"{ref.field_id}: the label volume has {labels.shape[0]} planes and the image "
                f"has {field.shape[0]}. They must be the same stack.")
        bin_factor = int(round(lgrid.dy_um / field.grid.dy_um))
        nuclei = nuclei_from_labels(labels, lgrid, params)
        usable = [n for n in nuclei if n.used]
        used = usable[:limit] if limit else usable
        for n in used:
            n.nucleus_uid = nucleus_uid(ref.stem, ref.position, n.label_id)
        log(f"      {len(nuclei)} nuclei pass the volume gate, {len(usable)} usable"
            + (f", {len(used)} kept by limit_nuclei" if len(used) != len(usable) else ""))

        # Plan every crop's bounds BEFORE reading a pixel, then work through them in increasing
        # z. Nuclei in one field overlap heavily in z, so z order turns the plane cache's job
        # into a sweep: each plane is read once and dropped when the window has passed it. Label
        # order is spatially arbitrary and makes the same planes be re-read as the window jumps.
        planned = [(n, crop_bounds(n.bbox_vox, labels.shape, bin_factor,
                                   params.detect.margin_um, lgrid, field.shape)) for n in used]
        planned.sort(key=lambda t: (t[1][3].start, t[1][3].stop))

        (out / "graphs").mkdir(parents=True, exist_ok=True)
        if keep_work:
            (out / "work").mkdir(parents=True, exist_ok=True)

        cache = PlaneCache(field, "actin", params.detect.destripe, params.detect.destripe_axis,
                           int(params.detect.plane_cache_mb) * 1024 ** 2)
        graphs, per_nucleus = [], []
        for n, bounds in planned:
            crop = make_crop(cache, labels, n.label_id, bounds, bin_factor, field.grid,
                             params, n.nucleus_uid)
            enh = enhance(crop, params)
            fg, thr = binarize(enh["response"], enh["null_level"], crop.mask, params,
                               vol=crop.actin_flat, sigma_vox=enh["sigma_vox"],
                               spacing=crop.grid.spacing)
            tr = trace(fg, crop.grid, params)
            fil = build_filaments(tr["skeleton"], tr["grid"], crop.actin_flat, crop.grid,
                                  params, noise_sd=crop.noise_sd)
            fil["min_branch_um"] = tr["min_branch_um"]
            G = build_graph(n, crop, fil, params, field)
            G.graph.update({
                "noise_sd": crop.noise_sd, "null_level": enh["null_level"],
                "separation": enh["separation"], "raw_median": crop.raw_median,
                "n_puncta_rejected": tr["n_puncta_rejected"],
                "n_short_components_dropped": tr["n_short_components_dropped"],
                "fg_voxels": int(fg.sum()),
                "fg_fraction": float(fg.sum() / max(crop.mask.sum(), 1)),
                "threshold_low": thr["low"], "threshold_high": thr["high"],
                "thinned_by_nms": thr["nms"],
            })
            import networkx as nx
            nx.write_graphml(G, out / "graphs" / f"{n.nucleus_uid}.graphml")
            graphs.append(G)
            per_nucleus.append(nucleus_row(G))
            if keep_work:
                _save_work(out / "work" / f"{n.nucleus_uid}.npz", crop, fg,
                           enh["null_level"], keep_work)
            if qc:
                _validate.trace_overlay(crop, tr, out / "qc" / "traces" /
                                        f"{n.nucleus_uid}.png")
        log(f"      {len(graphs)} graphs in {time.time() - t0:.0f}s; {cache.summary}")
        if qc and graphs:
            _validate.field_figure(field, labels, lgrid, graphs, params,
                                   out / "qc" / f"{ref.field_id}.png")
    return {"rows": per_nucleus,
            "edges": [r for G in graphs for r in edge_rows(G)],
            "summary": _validate.field_summary(ref, per_nucleus, time.time() - t0)}


def run_folder(input_dir, output_dir, params: Params | None = None, labels_dir=None,
               fields=None, positions=None, limit_nuclei=None, keep_work: str = WORK_SMALL,
               qc: bool = True, multi_series: bool = False, exclude=(), log=print) -> dict:
    """Every field in `input_dir` -> the output tree documented in the README."""
    params = params or Params()
    out = Path(output_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    (out / "qc").mkdir(exist_ok=True)
    if qc:
        (out / "qc" / "traces").mkdir(exist_ok=True)

    # Written FIRST, so a run that dies half way still says what it was doing.
    (out / "run_parameters.json").write_text(json.dumps(
        {**params.to_dict(), "input_dir": str(input_dir), "labels_dir": str(labels_dir),
         "keep_work": keep_work, "multi_series": bool(multi_series)}, indent=2, sort_keys=True))

    refs, skipped = list_fields(input_dir, params.optics, params.channels,
                                exclude=exclude, multi_series=multi_series)
    if fields:
        refs = [r for r in refs if r.field_id in set(fields) or r.stem in set(fields)]
    if positions is not None:
        refs = [r for r in refs if r.position in set(positions)]
    if not refs:
        raise RuntimeError(
            f"no usable fields in {input_dir}. {len(skipped)} file(s) were skipped; the first "
            f"reason was: {skipped[0]['reason'] if skipped else 'no .nd2 files found'}")

    notes = params.scales(strict=False)["notes"]
    log(f"antenna3d: {len(refs)} field(s), {len(skipped)} file(s) skipped")
    log(f"  scales digest {params.optics.digest()}"
        + ("".join(f"\n  note: {n}" for n in notes) if notes else ""))

    rows, edges, summaries, failures = [], [], [], []
    for i, ref in enumerate(refs, 1):
        log(f"[{i}/{len(refs)}] {ref.field_id}"
            + (f"  ({ref.condition})" if ref.condition else ""))
        try:
            r = run_field(ref, params, out, labels_dir, keep_work, qc, limit_nuclei,
                          multi_series, log)
        except Exception as e:                       # recorded, never fatal to the other fields
            log(f"      FAILED: {type(e).__name__}: {e}")
            failures.append({"field_id": ref.field_id, "file": str(ref.path),
                             "position": ref.position, "error": f"{type(e).__name__}: {e}",
                             "traceback": traceback.format_exc()})
            continue
        rows.extend(r["rows"])
        edges.extend(r["edges"])
        summaries.append(r["summary"])

    if rows:
        _write_table(pd.DataFrame(rows), out, "antenna_nuclei")
    if edges:
        _write_table(pd.DataFrame(edges), out, "antenna_edges")
    if summaries:
        pd.DataFrame(summaries).to_csv(out / "field_summary.csv", index=False)
    if skipped:
        pd.DataFrame(skipped).to_csv(out / "skipped_files.csv", index=False)
    if failures:
        pd.DataFrame(failures).to_csv(out / "failures.csv", index=False)
        log(f"\n{len(failures)} field(s) FAILED - see {out / 'failures.csv'}")

    nd = pd.DataFrame(rows)
    if len(nd):
        log(f"\n{len(nd)} nuclei over {nd.field_id.nunique()} field(s), {len(edges)} edges")
        log(f"  detected in {100 * nd.detected.mean():.0f}%; "
            f"median length density {nd.length_density_um_per_um3.median():.4f} um/um3")
        log(f"  QC: {out / 'qc'} - look at these before you use the table")
    return {"nuclei": nd, "edges": pd.DataFrame(edges),
            "field_summary": pd.DataFrame(summaries), "failures": pd.DataFrame(failures)}
