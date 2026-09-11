"""
The run driver: walks files and fields, calls the core functions in order,
and writes the output tree.

Everything here is bookkeeping. The science lives in segment.py,
quantify.py, export.py and validate.py; this module only sequences them
and handles per-field failures without losing the rest of the run.

Output tree
    <outdir>/
        nuclei_measurements.csv   one row per nucleus, all channels
        field_summary.csv         one row per field, QC numbers
        run_parameters.json       every parameter used, for reproducibility
        qc/                       one validation figure per field
        nucleus_boxes/            one OME-TIFF per nucleus  (if enabled)
        nucleus_boxes_index.csv   nucleus -> file map        (if enabled)
"""

import glob
import json
import os
import time
import traceback
from dataclasses import asdict

import numpy as np
import pandas as pd

from . import io as n3io
from .segment import SegParams, segment_nuclei
from .quantify import quantify_nuclei
from .export import export_nucleus_boxes
from .validate import validation_figure


def find_nd2(root, recursive=True):
    """All .nd2 files under `root` (or `root` itself if it is a file)."""
    if os.path.isfile(root):
        return [root]
    pattern = "**/*.nd2" if recursive else "*.nd2"
    return sorted(glob.glob(os.path.join(root, pattern), recursive=recursive))


def condition_from_path(path, root):
    """
    Condition label for a file: the name of its immediate parent folder
    relative to `root`, or 'root' for files sitting directly in it.

    Matches the usual layout of one folder per treatment.
    """
    rel = os.path.relpath(os.path.dirname(os.path.abspath(path)),
                          os.path.abspath(root))
    return "root" if rel in (".", "") else rel.replace(os.sep, "/")


def process_field(stack, dna_channel, params, outdir, save_qc=True,
                  save_boxes=False, box_pad_um=1.0, box_include_mask=True,
                  extra_columns=None, min_blob_um3=5.0):
    """
    Segment, quantify and (optionally) export one field.

    Returns (measurements, box_index, field_qc) -- any of which may be None
    or empty when the field contains no nuclei.
    """
    vol = stack.channel(dna_channel)
    labels, props = segment_nuclei(vol, stack.voxel_um, params)

    meas = pd.DataFrame()
    if len(props):
        quant = quantify_nuclei(stack, labels, dna_channel,
                                extra_columns=extra_columns)
        # props and quant overlap on geometry columns; keep one copy
        shared = [c for c in quant.columns if c in props.columns and c != "label"]
        meas = props.merge(quant.drop(columns=shared), on="label", how="left")

        front = [c for c in ("nucleus_uid", "file", "position", "field", "label")
                 if c in meas.columns]
        meas = meas[front + [c for c in meas.columns if c not in front]]

    boxes = pd.DataFrame()
    if save_boxes and len(props):
        boxes = export_nucleus_boxes(
            stack, labels, os.path.join(outdir, "nucleus_boxes"),
            pad_um=box_pad_um, include_mask=box_include_mask,
            extra_columns=extra_columns,
        )

    qc = None
    if save_qc:
        fig_name = stack.name.replace("#", "_").replace(".nd2", "") + "_qc.png"
        qc = validation_figure(stack, labels, props, dna_channel,
                               os.path.join(outdir, "qc", fig_name),
                               min_blob_um3=min_blob_um3)
        if extra_columns:
            qc.update(extra_columns)

    return meas, boxes, qc


def run(input_dir, outdir, dna_channel, params=None, positions=None,
        recursive=True, save_qc=True, save_boxes=False, box_pad_um=1.0,
        box_include_mask=True, label_conditions=True, min_blob_um3=5.0,
        verbose=True):
    """
    Run the full pipeline over a folder of .nd2 files.

    Parameters
    ----------
    input_dir : folder of .nd2 files (searched recursively by default), or
                a single .nd2 file
    outdir : output folder; created if absent
    dna_channel : channel name (e.g. 'DAPI') or index. CHECK THIS FIRST with
                  io.describe_file -- a wrong channel yields a plausible,
                  entirely meaningless table.
    params : SegParams or None
    positions : list of stage-position indices, or None for all of them
    save_boxes : write one 3D OME-TIFF per nucleus
    label_conditions : add a `condition` column from the parent folder name

    Returns
    -------
    dict with the measurement table, field summary and box index.
    """
    p = params or SegParams()
    os.makedirs(outdir, exist_ok=True)

    paths = find_nd2(input_dir, recursive)
    if not paths:
        raise FileNotFoundError(f"no .nd2 files under {input_dir}")

    jobs = []
    for fp in paths:
        npos = n3io.n_positions(fp)
        wanted = range(npos) if positions is None else [
            q for q in positions if q < npos]
        jobs.extend((fp, int(q)) for q in wanted)

    if verbose:
        print(f"{len(jobs)} field(s) from {len(paths)} file(s) -> {outdir}",
              flush=True)

    meas_all, box_all, qc_all, failures = [], [], [], []
    t0 = time.time()

    for k, (fp, pos) in enumerate(jobs, 1):
        extra = ({"condition": condition_from_path(fp, input_dir)}
                 if label_conditions and os.path.isdir(input_dir) else None)
        try:
            stack = n3io.load_field(fp, pos)
            if dna_channel not in stack.channels and isinstance(dna_channel, str):
                raise KeyError(f"channel {dna_channel!r} not in {stack.channels}")

            meas, boxes, qc = process_field(
                stack, dna_channel, p, outdir, save_qc=save_qc,
                save_boxes=save_boxes, box_pad_um=box_pad_um,
                box_include_mask=box_include_mask, extra_columns=extra,
                min_blob_um3=min_blob_um3)

            if len(meas):
                meas_all.append(meas)
            if len(boxes):
                box_all.append(boxes)
            if qc:
                qc_all.append(qc)

            if verbose:
                print(f"  [{k:3d}/{len(jobs)}] {stack.name:42s} "
                      f"{len(meas):3d} nuclei   ({time.time() - t0:5.0f}s)",
                      flush=True)

        except Exception as exc:                      # noqa: BLE001
            failures.append(dict(file=os.path.basename(fp), position=pos,
                                 error=f"{type(exc).__name__}: {exc}",
                                 traceback=traceback.format_exc()))
            if verbose:
                print(f"  [{k:3d}/{len(jobs)}] FAILED {os.path.basename(fp)} "
                      f"p{pos}: {type(exc).__name__}: {exc}", flush=True)

    table = (pd.concat(meas_all, ignore_index=True) if meas_all
             else pd.DataFrame())
    summary = pd.DataFrame(qc_all)
    box_index = (pd.concat(box_all, ignore_index=True) if box_all
                 else pd.DataFrame())

    table.to_csv(os.path.join(outdir, "nuclei_measurements.csv"), index=False)
    if len(summary):
        summary.to_csv(os.path.join(outdir, "field_summary.csv"), index=False)
    if len(box_index):
        box_index.to_csv(os.path.join(outdir, "nucleus_boxes_index.csv"),
                         index=False)
    if failures:
        pd.DataFrame(failures).to_csv(os.path.join(outdir, "failures.csv"),
                                      index=False)

    with open(os.path.join(outdir, "run_parameters.json"), "w") as fh:
        json.dump(dict(input_dir=os.path.abspath(input_dir),
                       outdir=os.path.abspath(outdir),
                       dna_channel=dna_channel,
                       positions=positions,
                       save_boxes=save_boxes, box_pad_um=box_pad_um,
                       box_include_mask=box_include_mask,
                       n_files=len(paths), n_fields=len(jobs),
                       n_failed=len(failures),
                       n_nuclei=int(len(table)),
                       elapsed_s=round(time.time() - t0, 1),
                       seg_params=asdict(p)), fh, indent=2)

    if verbose:
        print(f"\n{len(table)} nuclei from {len(jobs) - len(failures)} field(s) "
              f"in {time.time() - t0:.0f}s", flush=True)
        if failures:
            print(f"{len(failures)} field(s) failed -- see failures.csv",
                  flush=True)
        if len(summary):
            worst = summary.nlargest(3, "unsegmented_fraction")
            print("\nfields with most unsegmented DNA (check these figures):",
                  flush=True)
            print(worst[["field", "n_nuclei", "unsegmented_fraction"]]
                  .to_string(index=False), flush=True)

    return dict(measurements=table, field_summary=summary,
                box_index=box_index, failures=pd.DataFrame(failures))
