"""
Parallel field processing under an explicit memory budget.

One field of this data is ~100 MB on disk but peaks near 3 GB in RAM while
it is being segmented, so the number of workers that fits in memory is
smaller than the number of cores. Running `n_workers = cpu_count()` swaps
the machine to a standstill; running one worker wastes thirteen cores.

This module decides the number for you, in two stages:

    1. `estimate_peak_bytes`  predicts peak RAM per field from the file
       geometry alone (metadata read, no pixels).
    2. `run_parallel`         then MEASURES the true peak on the first
       field and re-sizes the pool from the measurement before starting
       the rest. The calibration field is a real field: its results are
       kept, nothing is processed twice.

Where the memory goes
    It used to go almost entirely into the distance transform in
    `segment._distance_map`. That step replicates the first and last slices
    by `pad = clip(round(12 um / dz), 1, 64)` planes so nuclei clipped at the
    slab faces are treated as interior, and `distance_transform_edt` then
    works in float64: at dz = 0.1 um the pad saturates at 64 planes each
    side, so a 31-plane stack was transformed as 159 planes, 5.1x the
    original volume, for a measured 7.1-7.5 GB peak.

    `SegParams.blockwise_distance` (default True) now computes that
    transform one connected component at a time on small lateral crops,
    which is provably the same array -- see `segment._distance_map_blockwise`
    -- and brings the measured peak to 2.8-3.0 GB on the same fields. What
    remains is spread across the DoG, the boundary refinement and the
    convex hulls in `regionprops`, none of them dominant, so peak now
    scales with the plain voxel count.
"""

import os
import platform
import resource
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

import numpy as np
import pandas as pd

from . import io as n3io
from .segment import SegParams


# ======================================================================
# memory accounting
# ======================================================================

PEAK_BYTES_PER_VOXEL = 170.0
"""
Calibration constant for `estimate_peak_bytes`, default (blockwise) path.

Measured in clean subprocesses on two 31 x 716 x 794 fields (17.6e6 voxels,
dz = 0.1 um): 3.03 GB and 2.82 GB peak resident, i.e. 166 and 154 bytes per
voxel once the raw stack is accounted separately. Taken at the top of that
range.
"""

PEAK_BYTES_PER_PADDED_VOXEL = 82.0
"""
The same, for `SegParams.blockwise_distance = False`.

There the distance transform runs on the whole z-padded volume, so peak
scales with the PADDED voxel count instead: the same two fields peaked at
7.07 GB and 7.50 GB over 90.4e6 padded voxels, i.e. 77 and 82 bytes each.

Both constants are properties of the algorithm rather than of the machine,
but they are empirical numbers from one dataset. `run_parallel` measures
the real value on the first field and overrides them, so they only have to
be in the right ballpark.
"""


def _maxrss_bytes():
    """Peak resident set size of this process, in bytes."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is bytes on macOS, kilobytes on Linux
    return float(rss) if platform.system() == "Darwin" else float(rss) * 1024.0


def total_memory_bytes():
    """Physical RAM installed."""
    try:
        return float(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES"))
    except (ValueError, OSError, AttributeError):
        return 8e9                       # conservative fallback


def available_memory_bytes():
    """
    RAM that can be claimed right now, without evicting the user's work.

    Falls back to a fraction of the installed total when the platform
    cannot be queried -- this only needs to be roughly right, and being
    wrong low costs throughput while being wrong high costs the machine.
    """
    system = platform.system()

    if system == "Linux":
        try:
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemAvailable:"):
                        return float(line.split()[1]) * 1024.0
        except OSError:
            pass

    elif system == "Darwin":
        try:
            out = subprocess.run(["vm_stat"], capture_output=True, text=True,
                                 timeout=10).stdout
            page = 4096
            if "page size of" in out:
                page = int(out.split("page size of")[1].split()[0])
            counts = {}
            for line in out.splitlines()[1:]:
                if ":" in line:
                    k, v = line.split(":", 1)
                    v = v.strip().rstrip(".")
                    if v.isdigit():
                        counts[k.strip()] = int(v)
            # free + inactive + speculative + purgeable are all reclaimable
            free = sum(counts.get(k, 0) for k in (
                "Pages free", "Pages inactive", "Pages speculative",
                "Pages purgeable"))
            if free:
                return float(free) * page
        except (OSError, ValueError, subprocess.SubprocessError):
            pass

    return 0.6 * total_memory_bytes()


def field_geometry(path):
    """Shape and sampling of one .nd2, from metadata only (no pixels read)."""
    d = n3io.describe_file(path)
    ny, nx = d["shape_yx"]
    return dict(path=path, nz=int(d["n_z"]), nc=len(d["channels"]),
                ny=int(ny), nx=int(nx), dz=float(d["voxel_um"][0]),
                n_positions=int(d["n_positions"]))


def estimate_peak_bytes(geom, params=None):
    """
    Predicted peak RAM for segmenting one field of this geometry.

    Two terms: the raw multichannel stack held in memory as uint16, and the
    segmentation working set, which dominates.

    Which voxel count drives the second term depends on how the distance
    transform is computed. With `blockwise_distance` (the default) it runs
    per connected component on small crops, so peak scales with the plain
    voxel count. With the whole-volume transform it scales with the
    z-PADDED count, which at fine dz is several times larger -- see the
    module docstring.
    """
    p = params or SegParams()
    raw = geom["nz"] * geom["nc"] * geom["ny"] * geom["nx"] * 2
    voxels = geom["nz"] * geom["ny"] * geom["nx"]

    if not p.watershed_split:
        return float(raw) + PEAK_BYTES_PER_VOXEL * voxels

    if p.blockwise_distance:
        return float(raw) + PEAK_BYTES_PER_VOXEL * voxels

    pad = int(np.clip(round(12.0 / geom["dz"]), 1, 64)) if p.ignore_z_border else 0
    padded_voxels = (geom["nz"] + 2 * pad) * geom["ny"] * geom["nx"]
    return float(raw) + PEAK_BYTES_PER_PADDED_VOXEL * padded_voxels


def plan_workers(peak_bytes_per_field, n_fields, max_workers=None,
                 memory_fraction=0.75, reserve_bytes=2e9):
    """
    How many fields can run at once without exhausting memory.

    Parameters
    ----------
    peak_bytes_per_field : predicted or measured peak RSS of one worker
    n_fields : total fields to process (never spawn more workers than work)
    max_workers : hard user cap, or None
    memory_fraction : share of currently-available RAM to spend
    reserve_bytes : held back for the parent process and the OS

    Returns a dict describing the decision, so callers can print or log it.
    """
    avail = available_memory_bytes()
    budget = max(avail * memory_fraction - reserve_bytes, peak_bytes_per_field)

    by_memory = max(1, int(budget // max(peak_bytes_per_field, 1.0)))
    by_cpu = max(1, (os.cpu_count() or 2) - 1)

    n = min(by_memory, by_cpu, max(1, n_fields))
    if max_workers:
        n = min(n, int(max_workers))

    return dict(n_workers=int(n), by_memory=int(by_memory), by_cpu=int(by_cpu),
                peak_per_field_gb=peak_bytes_per_field / 1e9,
                available_gb=avail / 1e9, budget_gb=budget / 1e9,
                total_gb=total_memory_bytes() / 1e9,
                limited_by=("memory" if by_memory < by_cpu else "cores"))


# ======================================================================
# worker
# ======================================================================

def _thread_limit_env():
    """
    Env vars that keep each worker single-threaded.

    Without this every worker starts its own BLAS thread pool and the
    processes fight for the same cores, which on an oversubscribed machine
    is slower than running fewer workers. Must be set BEFORE the child
    imports numpy, so they go into the parent's environment and are
    inherited -- setting them inside the child is too late.
    """
    return {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                             "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                             "VECLIB_MAXIMUM_THREADS")}


def _ping(x):
    """Trivial picklable payload for `pool_available`."""
    return x


_POOL_OK = None


def pool_available(timeout=20.0, verbose=False):
    """
    Can this platform actually RUN a process pool? Cached after the first call.

    Constructing the executor is not a sufficient test. There are two
    distinct failure modes, and only the first announces itself:

      * construction raises (a sandbox denying the POSIX semaphore probe
        `sysconf("SC_SEM_NSEMS_MAX")` gives PermissionError);
      * construction succeeds, but the workers never come up -- the manager
        thread marks the pool broken while the foreground call sits in
        `fut.result()` forever.

    The second mode would hang a run indefinitely, so this probe submits a
    trivial task with a deadline instead of trusting construction. Twenty
    seconds is generous for starting one interpreter and far cheaper than
    discovering the problem partway through a long batch.
    """
    global _POOL_OK
    if _POOL_OK is not None:
        return _POOL_OK

    ex = None
    try:
        ex = ProcessPoolExecutor(max_workers=1)
        _POOL_OK = (ex.submit(_ping, 1).result(timeout=timeout) == 1)
    except (PermissionError, OSError, BrokenProcessPool, TimeoutError) as exc:
        if verbose:
            print(f"  process pools unusable here ({type(exc).__name__}): "
                  f"running serially", flush=True)
        _POOL_OK = False
    except Exception as exc:                                  # noqa: BLE001
        if verbose:
            print(f"  process pool probe failed ({type(exc).__name__}: "
                  f"{exc}); running serially", flush=True)
        _POOL_OK = False
    finally:
        if ex is not None:
            try:
                ex.shutdown(wait=False, cancel_futures=True)
            except Exception:                                 # noqa: BLE001
                pass

    return _POOL_OK


def _map_jobs(job_tuples, n_workers, verbose=True, label="", isolate=False):
    """
    Run job tuples, in a process pool when the platform allows one.

    Whether a pool is usable is decided by `pool_available`, which probes
    with a deadline -- some hosts construct an executor happily and then
    never deliver a result. When pools are unusable the same jobs run in
    this process, one at a time; when a pool breaks partway through, the
    unfinished jobs finish here too.

    `isolate=True` forces a one-worker pool even for a single job. The
    calibration step needs that: `_maxrss_bytes` reads RUSAGE_SELF, a
    whole-process high-water mark, so measuring a field inside the
    long-lived driver would report the driver's history rather than the
    field's own peak, and the several GB would stay resident afterwards.

    Returns (results, used_pool). `used_pool` is False when the work ran
    in this process, which the caller needs in order to know whether a
    measured peak is trustworthy.
    """
    job_tuples = list(job_tuples)
    if not job_tuples:
        return [], False

    done, t0 = [], time.time()
    remaining = list(job_tuples)

    if (n_workers > 1 or isolate) and pool_available(verbose=verbose):
        futures = {}
        try:
            # max_tasks_per_child=1 returns every byte to the OS between
            # fields: at ~40 s per field the ~1 s respawn is under 3%, and it
            # removes any chance of memory creep pushing the pool past budget.
            with ProcessPoolExecutor(max_workers=max(1, n_workers),
                                     max_tasks_per_child=1) as ex:
                futures = {ex.submit(_process_one, j): j for j in remaining}
                for fut in as_completed(futures):
                    done.append(fut.result())
                    futures[fut] = None          # this job is accounted for
                    if verbose:
                        _log(done[-1], len(done), len(job_tuples), t0, label)
            return done, True

        except BrokenProcessPool as exc:
            # A worker died mid-run, most often to the OS memory killer --
            # i.e. the budget was too optimistic. Keep the fields that
            # finished and grind out the rest here rather than losing a
            # multi-hour batch.
            remaining = [j for j in futures.values() if j is not None]
            if verbose:
                print(f"  pool broke ({type(exc).__name__}); {len(done)} field(s) "
                      f"done, {len(remaining)} continuing serially", flush=True)

    for j in remaining:
        done.append(_process_one(j))
        if verbose:
            _log(done[-1], len(done), len(job_tuples), t0, label)
    return done, False


def _log(res, done, total, t0, label=""):
    tag = (f"{len(res['meas']):3d} nuclei" if res["ok"]
           else f"FAILED {res['failure']['error'][:48]}")
    print(f"  {label}[{done:3d}/{total}] {res['name']:42s} {tag}"
          f"   ({time.time() - t0:5.0f}s)", flush=True)


def _process_one(job):
    """
    Segment, quantify and export one field. Runs in a worker process.

    Returns plain data (DataFrames and dicts), never arrays: the label
    image and the stack stay in the worker and are freed when it exits.
    Failures are captured, not raised -- one bad field must not abort a
    run of two hundred.
    """
    from .pipeline import process_field          # imported here: keeps the
    from . import io as _io                      # parent import graph light

    (path, position, dna_channel, params, outdir, save_qc, save_boxes,
     box_pad_um, box_include_mask, extra_columns, min_blob_um3,
     reuse_boxes, cache_dir, save_labels) = job

    t0 = time.time()
    try:
        stack = _io.load_field(path, position)
        if isinstance(dna_channel, str) and dna_channel not in stack.channels:
            raise KeyError(f"channel {dna_channel!r} not in {stack.channels}")

        meas, boxes, qc = process_field(
            stack, dna_channel, params, outdir, save_qc=save_qc,
            save_boxes=save_boxes, box_pad_um=box_pad_um,
            box_include_mask=box_include_mask, extra_columns=extra_columns,
            min_blob_um3=min_blob_um3, reuse_boxes=reuse_boxes,
            cache_dir=cache_dir, save_labels=save_labels)

        return dict(ok=True, name=stack.name, meas=meas, boxes=boxes, qc=qc,
                    peak_bytes=_maxrss_bytes(), elapsed=time.time() - t0,
                    failure=None)

    except Exception as exc:                                  # noqa: BLE001
        return dict(ok=False, name=f"{os.path.basename(path)}#p{position:02d}",
                    meas=None, boxes=None, qc=None,
                    peak_bytes=_maxrss_bytes(), elapsed=time.time() - t0,
                    failure=dict(file=os.path.basename(path), position=position,
                                 error=f"{type(exc).__name__}: {exc}",
                                 traceback=traceback.format_exc()))


# ======================================================================
# driver
# ======================================================================

def run_parallel(jobs, dna_channel, params, outdir, save_qc=True,
                 save_boxes=False, box_pad_um=1.0, box_include_mask=True,
                 min_blob_um3=5.0, extra_for=None, n_workers="auto",
                 max_workers=None, memory_fraction=0.75, calibrate=True,
                 reuse_boxes=True, cache_dir=None, save_labels=True,
                 verbose=True):
    """
    Process a list of (path, position) jobs across processes.

    `n_workers="auto"` sizes the pool from measured memory use; pass an
    integer to override. `extra_for` is a callable (path -> dict) supplying
    the per-file constant columns, e.g. the condition label.

    Returns (results, plan) where `results` is the list of per-field dicts
    from `_process_one` and `plan` records how the pool was sized.
    """
    jobs = list(jobs)
    if not jobs:
        return [], {}

    def _job_tuple(path, pos):
        extra = extra_for(path) if extra_for else None
        return (path, pos, dna_channel, params, outdir, save_qc, save_boxes,
                box_pad_um, box_include_mask, extra, min_blob_um3,
                reuse_boxes, cache_dir, save_labels)

    # --- predict from geometry -----------------------------------------
    geom = field_geometry(jobs[0][0])
    predicted = estimate_peak_bytes(geom, params)
    plan = plan_workers(predicted, len(jobs), max_workers, memory_fraction)
    plan["source"] = "predicted"

    if verbose:
        print(f"memory: {plan['total_gb']:.0f} GB installed, "
              f"{plan['available_gb']:.0f} GB available, "
              f"budget {plan['budget_gb']:.0f} GB", flush=True)
        print(f"predicted peak {plan['peak_per_field_gb']:.1f} GB/field "
              f"-> {plan['n_workers']} worker(s) "
              f"(memory allows {plan['by_memory']}, cores allow {plan['by_cpu']})",
              flush=True)

    # keep BLAS single-threaded in the children (inherited at spawn)
    saved_env = {k: os.environ.get(k) for k in _thread_limit_env()}
    os.environ.update(_thread_limit_env())

    results = []
    try:
        # --- calibrate on one real field --------------------------------
        if calibrate and n_workers == "auto" and len(jobs) > 1:
            if verbose:
                print("calibrating on field 1 ...", flush=True)
            cal, isolated = _map_jobs([_job_tuple(*jobs[0])], 1,
                                      verbose=False, isolate=True)
            first = cal[0]
            results.append(first)

            measured = first["peak_bytes"]
            plan = plan_workers(measured, len(jobs) - 1, max_workers,
                                memory_fraction)
            plan.update(source="measured" if isolated else "measured-in-process",
                        calibration_isolated=isolated,
                        predicted_gb=predicted / 1e9,
                        calibration_s=first["elapsed"])
            jobs = jobs[1:]

            if verbose:
                # An in-process measurement is the driver's whole-process
                # high-water mark, so it can only over-state one field's
                # peak. That errs toward fewer workers, which is the safe
                # direction, but say so rather than passing it off as clean.
                note = "" if isolated else " [in-process: upper bound]"
                print(f"measured peak {measured / 1e9:.1f} GB/field{note} "
                      f"(predicted {predicted / 1e9:.1f}) in "
                      f"{first['elapsed']:.0f}s -> {plan['n_workers']} worker(s), "
                      f"limited by {plan['limited_by']}", flush=True)

        n = plan["n_workers"] if n_workers == "auto" else int(n_workers)
        n = max(1, min(n, len(jobs))) if jobs else 1
        plan["n_workers"] = n

        if not jobs:
            return results, plan

        # --- the rest ----------------------------------------------------
        rest, used_pool = _map_jobs([_job_tuple(p, q) for p, q in jobs],
                                    n, verbose=verbose)
        results.extend(rest)
        plan["used_pool"] = used_pool
    finally:
        for k, v in saved_env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    return results, plan
