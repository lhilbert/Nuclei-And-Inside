"""The driver: .nd2 + nucleus3d label volume in, spot tables out.

Nuclei are never segmented here. Each field's mask comes from nucleus3d's
`labels/<stem>_pNN.tif`; a field without one is recorded in failures.csv and
skipped.

Output tree
    <outdir>/
        run_parameters.json      written FIRST
        spots.csv                one row per spot (channel x nucleus x spot)
        spots_per_nucleus.csv    one row per nucleus x channel
        field_summary.csv        one row per field x channel
        failures.csv             fields that errored or had no labels
        spot_labels/             <stem>_pNN_<channel>.tif, spot label volumes (OME, um)
        qc/                      field overview + per-nucleus panels
"""
from __future__ import annotations

import json
import os
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from nucleus3d.core.pipeline import condition_from_path, find_nd2
from nucleus3d.core.quantify import nucleus_uid

from .detect import detect_in_nucleus, erode_mask_union, erosion_radius_px
from .features import spot_features, spot_intensities
from .io import ND2Field, labels_path_for, read_labels
from .params import SpotParams


def detect_field(images: dict, labels: np.ndarray, voxel_um, params: SpotParams,
                 source_path: str, position: int = 0, keep_crops: bool = False) -> dict:
    """Detect and measure spots in every nucleus of one field.

    `images` maps each ChannelSpec.name to its raw (Z, Y, X) volume; `labels`
    is the nucleus label volume on the same grid. Array-level, so it runs on
    synthetic data without an .nd2.
    """
    dz, dy, dx = voxel_um
    vox_um3 = dx * dx * dz
    dtype = np.dtype(params.image_dtype)
    specs = params.channels
    owner = SimpleNamespace(source_path=source_path, position=position)
    base = dict(file=os.path.basename(source_path), position=position)

    r_px = erosion_radius_px(params.erosion_um, dx) if params.erosion_um > 0 else 0
    used = erode_mask_union(labels > 0, r_px)

    spot_vols = {s.name: np.zeros(labels.shape, np.uint16) for s in specs}
    next_id = {s.name: 1 for s in specs}
    spot_rows, nuc_rows, crops = [], [], []

    for lab_id, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is None:
            continue
        uid = nucleus_uid(owner, lab_id)
        whole = labels[sl] == lab_id
        m_full = whole & used[sl]
        nuc_base = dict(base, nucleus_uid=uid, label=lab_id,
                        nucleus_voxels=int(whole.sum()),
                        nucleus_volume_um3=float(whole.sum()) * vox_um3,
                        mask_voxels=int(m_full.sum()),
                        mask_volume_um3=float(m_full.sum()) * vox_um3,
                        erosion_um=params.erosion_um, erosion_px=r_px)
        if not m_full.any():
            for s in specs:
                nuc_rows.append(dict(nuc_base, channel=s.name, nd2_index=s.index,
                                     status="empty_after_erosion", n_spots=0))
            continue

        # MATLAB crops to the bounding box of the ERODED nucleus (l.656-701)
        sub = ndi.find_objects(m_full.astype(np.uint8))[0]
        csl = tuple(slice(a.start + b.start, a.start + b.stop) for a, b in zip(sl, sub))
        mask = m_full[sub]
        origin = tuple(s.start for s in csl)
        qcrops = {s.name: images[s.name][csl].astype(dtype) for s in specs}
        nuc_mean = {s.name: float(qcrops[s.name][mask].mean(dtype=np.float64))
                    for s in specs}

        for s in specs:
            det = detect_in_nucleus(qcrops[s.name], mask, voxel_um, s,
                                    params.small_sigma_cutoff_px)
            lab = det["labels"]
            feats = spot_features(lab, voxel_um, origin)
            ratios = {}
            for q in specs:
                ratios[q.name], nuc_med = spot_intensities(
                    lab, qcrops[q.name], mask, q.quant_blur_um, dx)
                ratios[q.name + "_nuc_median"] = nuc_med

            n = len(feats)
            for k, f in enumerate(feats):
                row = dict(base, nucleus_uid=uid, label=lab_id, channel=s.name,
                           nd2_index=s.index, **f,
                           n_merged_components=int(det["n_merged"][k]),
                           field_spot_id=next_id[s.name] + k)
                for q in specs:
                    row[f"int_ratio_{q.name}"] = float(ratios[q.name][k])
                row.update(threshold=det["threshold"], erosion_um=params.erosion_um,
                           nucleus_mask_volume_um3=nuc_base["mask_volume_um3"])
                spot_rows.append(row)

            if n:
                region = spot_vols[s.name][csl]
                region[lab > 0] = (lab[lab > 0] + next_id[s.name] - 1).astype(np.uint16)
                next_id[s.name] += n

            vols = [f["volume_um3"] for f in feats]
            nrow = dict(nuc_base, channel=s.name, nd2_index=s.index, status="ok",
                        n_spots=n,
                        median_spot_volume_um3=float(np.median(vols)) if n else np.nan,
                        total_spot_volume_um3=float(np.sum(vols)),
                        thr_mean=det["thr_mean"], thr_sd=det["thr_sd"],
                        threshold=det["threshold"],
                        n_voxels_above_threshold=det["n_above_threshold"],
                        n_components_raw=det["n_components_raw"],
                        n_components_after_volume_gate=det["n_components_gated"],
                        qc_more_than_max_spots=bool(n > params.qc_max_spots))
            for q in specs:
                nrow[f"nuc_mean_{q.name}"] = nuc_mean[q.name]
            nuc_rows.append(nrow)
            if keep_crops:
                # projections only: full crops of every nucleus would be GBs
                crops.append(dict(nucleus_uid=uid, channel=s.name, origin=origin,
                                  raw=qcrops[s.name].max(axis=0),
                                  dog=det["dog"].max(axis=0), mask=mask.any(axis=0),
                                  spots=lab.max(axis=0), feats=feats,
                                  threshold=det["threshold"]))

    return dict(spots=pd.DataFrame(spot_rows), nuclei=pd.DataFrame(nuc_rows),
                spot_volumes=spot_vols, used_mask=used, crops=crops)


def _field_summary(nuclei: pd.DataFrame, field_info: dict, params: SpotParams) -> list:
    rows = []
    for s in params.channels:
        d = nuclei[nuclei["channel"] == s.name] if len(nuclei) else nuclei
        ok = d[d["status"] == "ok"] if len(d) else d
        rows.append(dict(field_info, channel=s.name, nd2_index=s.index,
                         n_nuclei=int(len(ok)),
                         n_nuclei_skipped=int(len(d) - len(ok)),
                         n_spots=int(ok["n_spots"].sum()) if len(ok) else 0,
                         mean_spots_per_nucleus=float(ok["n_spots"].mean()) if len(ok) else np.nan,
                         n_nuclei_with_0=int((ok["n_spots"] == 0).sum()) if len(ok) else 0,
                         n_nuclei_more_than_max=int(ok["qc_more_than_max_spots"].sum()) if len(ok) else 0,
                         # MATLAB fractionAboveFour > 0 marks the whole file invalid
                         # (l.2046-2050). Reported here, never acted on.
                         qc_any_nucleus_more_than_max=bool(ok["qc_more_than_max_spots"].any())
                         if len(ok) else False))
    return rows


def _write_spot_volume(path, vol, voxel_um):
    import tifffile
    dz, dy, dx = voxel_um
    tifffile.imwrite(path, vol, ome=True, compression="zlib", photometric="minisblack",
                     metadata={"axes": "ZYX", "PhysicalSizeX": dx, "PhysicalSizeXUnit": "µm",
                               "PhysicalSizeY": dy, "PhysicalSizeYUnit": "µm",
                               "PhysicalSizeZ": dz, "PhysicalSizeZUnit": "µm"})


def run_field(nd2_path, labels_dir, outdir, params: SpotParams, position: int = 0,
              qc: bool = True, save_spot_labels: bool = True, log=print) -> dict:
    """One field: read, detect, write the per-field artefacts, return the tables."""
    from . import qc as _qc
    out = Path(outdir)
    with ND2Field(nd2_path, position) as fld:
        lpath = labels_path_for(labels_dir, fld.stem, position)
        if not lpath.exists():
            raise FileNotFoundError(f"no nucleus3d label volume {lpath.name} in {labels_dir} "
                                    f"(segmentation is not re-run here)")
        labels = read_labels(lpath, fld.shape, fld.voxel_um)
        images = {s.name: fld.channel(s.index, s.expected_nd2_name) for s in params.channels}
        voxel = fld.voxel_um
        stem = fld.stem
    log(f"    labels {lpath.name}: {int(labels.max())} max label, voxel {voxel} um")

    res = detect_field(images, labels, voxel, params, str(nd2_path), position, keep_crops=qc)
    tag = f"{stem}_p{position:02d}"
    if save_spot_labels:
        (out / "spot_labels").mkdir(parents=True, exist_ok=True)
        for name, vol in res["spot_volumes"].items():
            _write_spot_volume(out / "spot_labels" / f"{tag}_{name}.tif", vol, voxel)
    if qc:
        _qc.field_overview(images, labels, res["spots"], params, voxel,
                           out / "qc" / f"{tag}_overview.png")
        _qc.nucleus_panels(res["crops"], params, out / "qc" / "nuclei")
    return res


def run_folder(input_path, labels_dir, outdir, params: SpotParams | None = None,
               positions=(0,), files=None, qc: bool = True, save_spot_labels: bool = True,
               log=print) -> dict:
    """Every .nd2 under `input_path` (or that one file) that has a label volume."""
    params = params or SpotParams()
    out = Path(outdir).expanduser()
    if (out / "nuclei_measurements.csv").exists():
        raise FileExistsError(f"{out} is a nucleus3d output folder; write spots elsewhere")
    out.mkdir(parents=True, exist_ok=True)
    (out / "run_parameters.json").write_text(json.dumps(
        dict(params=params.to_dict(), input_path=str(input_path), labels_dir=str(labels_dir),
             positions=list(positions)), indent=2))

    paths = find_nd2(str(input_path))
    if files:
        paths = [p for p in paths if os.path.splitext(os.path.basename(p))[0] in set(files)]
    root = input_path if os.path.isdir(input_path) else os.path.dirname(input_path)

    spots, nuclei, summary, failures = [], [], [], []
    t0 = time.time()
    for k, p in enumerate(paths, 1):
        for pos in positions:
            t1 = time.time()
            log(f"[{k}/{len(paths)}] {os.path.basename(p)} p{pos:02d}")
            try:
                res = run_field(p, labels_dir, out, params, pos, qc=qc,
                                save_spot_labels=save_spot_labels, log=log)
                cond = condition_from_path(p, root)
                for df in (res["spots"], res["nuclei"]):
                    if len(df):
                        df.insert(2, "condition", cond)
                spots.append(res["spots"])
                nuclei.append(res["nuclei"])
                summary += _field_summary(res["nuclei"], dict(
                    file=os.path.basename(p), position=pos, condition=cond,
                    seconds=round(time.time() - t1, 1)), params)
                log(f"    " + ", ".join(
                    f"{r['channel']}: {r['n_spots']} spots in {r['n_nuclei']} nuclei"
                    for r in summary[-len(params.channels):]))
            except Exception as exc:                          # noqa: BLE001
                failures.append(dict(file=os.path.basename(p), position=pos,
                                     error=f"{type(exc).__name__}: {exc}",
                                     traceback=traceback.format_exc()))
                log(f"    FAILED: {type(exc).__name__}: {exc}")

    tables = dict(
        spots=pd.concat(spots, ignore_index=True) if spots else pd.DataFrame(),
        nuclei=pd.concat(nuclei, ignore_index=True) if nuclei else pd.DataFrame(),
        field_summary=pd.DataFrame(summary), failures=pd.DataFrame(failures))
    tables["spots"].to_csv(out / "spots.csv", index=False)
    tables["nuclei"].to_csv(out / "spots_per_nucleus.csv", index=False)
    tables["field_summary"].to_csv(out / "field_summary.csv", index=False)
    tables["failures"].to_csv(out / "failures.csv", index=False)
    log(f"done in {time.time() - t0:.0f} s; {len(failures)} failure(s) -> {out}")
    return tables
