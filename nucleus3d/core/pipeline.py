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
import hashlib
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
from .export import export_nucleus_boxes, params_fingerprint
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


def code_fingerprint():
    """
    Short hash of the source that determines a measurement table.

    Cached results must not survive an edit to the code that produced them,
    and this project has no release cadence to hang a version number on, so
    the fingerprint is taken over the modules themselves. Editing a metric
    in `quantify.py` or a threshold default in `segment.py` therefore
    invalidates every cached field automatically.
    """
    h = hashlib.sha1()
    here = os.path.dirname(os.path.abspath(__file__))
    for name in ("segment.py", "quantify.py", "io.py"):
        with open(os.path.join(here, name), "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()[:16]


def field_cache_key(stack, dna_channel, params, min_blob_um3):
    """Identity of a field's measurement table: inputs, settings, code."""
    try:
        st = os.stat(stack.source_path)
        src = (int(st.st_size), int(st.st_mtime))
    except OSError:
        src = (-1, -1)
    payload = (os.path.abspath(stack.source_path), src, int(stack.position),
               dna_channel, float(min_blob_um3),
               params_fingerprint(params), code_fingerprint())
    return hashlib.sha1(repr(payload).encode()).hexdigest()[:16]


def _write_cache(path, df):
    """
    Write a cached measurement table that reads back bit-identical.

    Two details are needed for that, and neither is the default.

    `%.17g` is the shortest decimal form guaranteed to round-trip float64;
    pandas' default repr is shorter and loses the last bits.

    The column dtypes go in a sidecar, because CSV carries no schema: a
    column of whole-valued floats (a median of integer pixel values, say)
    comes back as int64, and concatenating such a cached field with a
    freshly computed one would then produce object columns. The sidecar is
    JSON so a stale cache stays inspectable by eye.
    """
    df.to_csv(path, index=False, float_format="%.17g")
    with open(path + ".schema.json", "w") as fh:
        json.dump({c: str(t) for c, t in df.dtypes.items()}, fh, indent=1)


def _read_cache(path):
    """
    Read a cached table back, or None if it cannot be trusted.

    `float_precision="round_trip"` selects the correctly-rounded parser;
    pandas' default C parser is fast but not exact, which silently
    introduces differences of order 1e-14 -- small, but enough that a
    cached run stops being reproducible.

    A cache is disposable by construction, so any problem here is a miss
    rather than an error: the field is simply re-segmented.
    """
    try:
        df = pd.read_csv(path, float_precision="round_trip")
        schema_path = path + ".schema.json"
        if os.path.exists(schema_path):
            with open(schema_path) as fh:
                df = df.astype(json.load(fh))
        return df
    except Exception:                                         # noqa: BLE001
        return None


def process_field(stack, dna_channel, params, outdir, save_qc=True,
                  save_boxes=False, box_pad_um=1.0, box_include_mask=True,
                  extra_columns=None, min_blob_um3=5.0, reuse_boxes=True,
                  cache_dir=None):
    """
    Segment, quantify and (optionally) export one field.

    Returns (measurements, box_index, field_qc) -- any of which may be None
    or empty when the field contains no nuclei.

    With `cache_dir`, the measurement table is written there under a key
    covering the source file (path, size, mtime), the stage position, the
    DNA channel, the segmentation parameters AND a hash of the module
    sources, and is read back instead of re-segmenting on a later run. This
    is what makes a resumed or repeated run cheap: segmentation is ~38 s per
    field against ~0.4 s for the box export, so reusing boxes alone saves
    almost nothing.

    The cache is bypassed whenever `save_qc` or `save_boxes` is set, because
    both need the label image, which is not cached -- only the table is.
    """
    cache_path = None
    if cache_dir and not (save_qc or save_boxes):
        key = field_cache_key(stack, dna_channel, params, min_blob_um3)
        os.makedirs(cache_dir, exist_ok=True)
        cache_path = os.path.join(
            cache_dir,
            f"{stack.name.replace('#', '_').replace('.nd2', '')}__{key}.csv")
        if os.path.exists(cache_path):
            cached = _read_cache(cache_path)
            if cached is not None:
                return cached, pd.DataFrame(), None

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
            extra_columns=extra_columns, params=params, reuse=reuse_boxes,
        )

    qc = None
    if save_qc:
        fig_name = stack.name.replace("#", "_").replace(".nd2", "") + "_qc.png"
        qc = validation_figure(stack, labels, props, dna_channel,
                               os.path.join(outdir, "qc", fig_name),
                               min_blob_um3=min_blob_um3)
        if extra_columns:
            qc.update(extra_columns)

    if cache_path is not None:
        # written last, so an interrupted field leaves no cache entry
        _write_cache(cache_path, meas)

    return meas, boxes, qc


def run(input_dir, outdir, dna_channel, params=None, positions=None,
        recursive=True, save_qc=True, save_boxes=False, box_pad_um=1.0,
        box_include_mask=True, label_conditions=True, min_blob_um3=5.0,
        n_workers=1, max_workers=None, memory_fraction=0.75,
        reuse_boxes=True, use_cache=True, verbose=True):
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
    n_workers : 1 for serial (default), an integer for a fixed pool, or
                "auto" to size the pool from measured memory use. Peak RAM
                per field is several GB, so "auto" is usually far below the
                core count -- see nucleus3d.parallel.
    max_workers : hard cap on the pool when n_workers="auto"
    memory_fraction : share of available RAM "auto" is allowed to spend

    Returns
    -------
    dict with the measurement table, field summary and box index.
    """
    p = params or SegParams()
    os.makedirs(outdir, exist_ok=True)
    cache_dir = (os.path.join(outdir, "field_cache")
                 if use_cache and not (save_qc or save_boxes) else None)

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
    parallel_plan = None

    def _extra_for(fp):
        return ({"condition": condition_from_path(fp, input_dir)}
                if label_conditions and os.path.isdir(input_dir) else None)

    if n_workers != 1:
        from .parallel import run_parallel

        results, parallel_plan = run_parallel(
            jobs, dna_channel, p, outdir, save_qc=save_qc,
            save_boxes=save_boxes, box_pad_um=box_pad_um,
            box_include_mask=box_include_mask, min_blob_um3=min_blob_um3,
            extra_for=_extra_for, n_workers=n_workers,
            max_workers=max_workers, memory_fraction=memory_fraction,
            reuse_boxes=reuse_boxes, cache_dir=cache_dir, verbose=verbose)

        for res in results:
            if not res["ok"]:
                failures.append(res["failure"])
                continue
            if res["meas"] is not None and len(res["meas"]):
                meas_all.append(res["meas"])
            if res["boxes"] is not None and len(res["boxes"]):
                box_all.append(res["boxes"])
            if res["qc"]:
                qc_all.append(res["qc"])

        jobs_iter = []                    # the loop below has nothing left to do
    else:
        jobs_iter = list(enumerate(jobs, 1))

    for k, (fp, pos) in jobs_iter:
        extra = _extra_for(fp)
        try:
            stack = n3io.load_field(fp, pos)
            if dna_channel not in stack.channels and isinstance(dna_channel, str):
                raise KeyError(f"channel {dna_channel!r} not in {stack.channels}")

            meas, boxes, qc = process_field(
                stack, dna_channel, p, outdir, save_qc=save_qc,
                save_boxes=save_boxes, box_pad_um=box_pad_um,
                box_include_mask=box_include_mask, extra_columns=extra,
                min_blob_um3=min_blob_um3, reuse_boxes=reuse_boxes,
                cache_dir=cache_dir)

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
                       n_workers=n_workers, parallel_plan=parallel_plan,
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
