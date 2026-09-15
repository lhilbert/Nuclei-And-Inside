"""
Export of individual nuclei as 3D OME-TIFF boxes.

One file per nucleus: the bounding box padded by a physical margin, with
every channel and every z-slice retained, and the segmentation mask
appended as a final channel so it travels with the pixels.
"""

import dataclasses
import hashlib
import json
import os
import re
from xml.sax.saxutils import unescape

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops

from .quantify import nucleus_uid


# Provenance fields that decide whether an existing box may be reused.
# Everything here changes the PIXELS in the file; anything not here (the
# free-text note, condition labels) does not.
_REUSE_KEYS = ("source_path", "source_bytes", "source_mtime", "position",
               "nucleus_label", "bbox_yx_in_source", "pad_um", "channels",
               "voxel_um_zyx", "mask_sha1", "seg_params_sha1")


def _sha1(obj):
    return hashlib.sha1(repr(obj).encode()).hexdigest()[:16]


def params_fingerprint(params):
    """
    Short hash of the segmentation parameters that produced a mask.

    Stored in each box so a later run can tell whether an existing file was
    extracted under the same settings. `None` hashes to a sentinel rather
    than to the default parameters: not knowing which settings were used is
    different from knowing they were the defaults, and must not silently
    compare equal to anything.
    """
    if params is None:
        return "unknown"
    d = dataclasses.asdict(params) if dataclasses.is_dataclass(params) else dict(params)
    return _sha1(sorted(d.items()))


def box_is_current(path, expected):
    """
    Is the box at `path` already exactly what we would write now?

    Compares the embedded provenance against `expected` on the fields in
    `_REUSE_KEYS`. Returns (bool, reason) -- the reason names the first
    field that differs, which is what you want in a log when a run you
    expected to be a no-op starts rewriting files.
    """
    if not os.path.exists(path):
        return False, "absent"
    try:
        prov = read_box_provenance(path)
    except Exception as exc:                                  # noqa: BLE001
        return False, f"unreadable ({type(exc).__name__})"
    if not prov:
        return False, "no provenance"

    for k in _REUSE_KEYS:
        if k not in prov:
            return False, f"provenance predates this check ({k} missing)"
        if prov[k] != expected[k]:
            return False, f"{k} differs"
    return True, "match"


def export_nucleus_boxes(stack, labels, outdir, pad_um=1.0, include_mask=True,
                         compression="zlib", extra_columns=None, params=None,
                         reuse=True):
    """
    Write one OME-TIFF per nucleus.

    The crop is the nucleus bounding box expanded by `pad_um` laterally and
    clipped at the image edge. Z is NEVER cropped: these stacks are already
    thin slabs through the nuclei, so every slice carries signal.

    Written as ZCYX with PhysicalSize metadata, so Fiji and BioFormats open
    them correctly scaled. With `include_mask`, a binary "NucleusMask"
    channel is appended.

    Provenance (source file, position, label, bounding box in source
    coordinates) is embedded as JSON in the OME Description -- read it back
    with `read_box_provenance`.

    Reuse
    -----
    With `reuse=True` (default), a box already on disk is left alone when
    its provenance shows it would be written identically -- same source
    file (path, size, mtime), same crop, same padding, same channels, same
    voxel size, same segmentation parameters, and the same mask content by
    SHA-1. Re-running an interrupted export then costs a header read per
    nucleus instead of a re-encode.

    The mask hash is what makes this safe. Parameters alone are not enough:
    two runs can share a `SegParams` and still produce different masks if
    the input changed, and a box whose pixels no longer match the current
    segmentation is worse than no box at all. Hashing the mask that would
    be written makes the check content-addressed, so the parameter
    fingerprint only has to explain WHY something differs.

    Note this skips re-EXTRACTION, not re-segmentation: the labels must
    already exist to compute the hash. What it saves is the encode and
    write of unchanged crops, which is the slow, disk-hungry part.

    Returns
    -------
    DataFrame mapping each nucleus to its file, sharing `nucleus_uid` with
    the measurement table so the two join directly. The `reused` column
    flags rows served from disk, and `reuse_reason` says why a box was
    rewritten.
    """
    os.makedirs(outdir, exist_ok=True)

    dz, dy, dx = stack.voxel_um
    pad_y, pad_x = int(round(pad_um / dy)), int(round(pad_um / dx))
    nz, _, ny, nx = stack.data.shape

    chan_names = list(stack.channels) + (["NucleusMask"] if include_mask else [])

    try:
        src_stat = os.stat(stack.source_path)
        src_bytes, src_mtime = int(src_stat.st_size), int(src_stat.st_mtime)
    except OSError:
        src_bytes, src_mtime = -1, -1

    seg_fp = params_fingerprint(params)

    records = []
    for r in regionprops(labels):
        _, y0, x0, _, y1, x1 = r.bbox
        y0, y1 = max(0, y0 - pad_y), min(ny, y1 + pad_y)
        x0, x1 = max(0, x0 - pad_x), min(nx, x1 + pad_x)

        mask_crop = labels[:, y0:y1, x0:x1] == r.label

        uid = nucleus_uid(stack, r.label)
        fname = f"{uid}.ome.tif"
        path = os.path.join(outdir, fname)

        provenance = {
            "nucleus_uid": uid,
            "source_file": os.path.basename(stack.source_path),
            "source_path": stack.source_path,
            "source_bytes": src_bytes,
            "source_mtime": src_mtime,
            "position": int(stack.position),
            "nucleus_label": int(r.label),
            "bbox_yx_in_source": [int(y0), int(x0), int(y1), int(x1)],
            "pad_um": pad_um,
            "voxel_um_zyx": [dz, dy, dx],
            "channels": chan_names,
            "mask_sha1": hashlib.sha1(np.packbits(mask_crop).tobytes()).hexdigest(),
            "seg_params_sha1": seg_fp,
            "note": ("z not cropped; source stack is a thin slab and nuclei "
                     "are truncated top and bottom"),
        }
        if extra_columns:
            provenance.update({k: str(v) for k, v in extra_columns.items()})

        current, reason = (box_is_current(path, provenance) if reuse
                           else (False, "reuse disabled"))

        if not current:
            sub = stack.data[:, :, y0:y1, x0:x1]            # (Z, C, y, x)
            if include_mask:
                m = mask_crop.astype(stack.data.dtype) * np.iinfo(stack.data.dtype).max
                sub = np.concatenate([sub, m[:, None]], axis=1)

            tifffile.imwrite(
                path, sub, ome=True, compression=compression,
                metadata={
                    "axes": "ZCYX",
                    "PhysicalSizeX": dx, "PhysicalSizeXUnit": "µm",
                    "PhysicalSizeY": dy, "PhysicalSizeYUnit": "µm",
                    "PhysicalSizeZ": dz, "PhysicalSizeZUnit": "µm",
                    "Channel": {"Name": chan_names},
                    "Description": json.dumps(provenance),
                },
            )

        rec = dict(nucleus_uid=uid, reused=bool(current), reuse_reason=reason,
                   file=os.path.basename(stack.source_path),
                   position=stack.position, label=int(r.label),
                   ome_tiff=fname,
                   box_y0=int(y0), box_x0=int(x0), box_y1=int(y1), box_x1=int(x1),
                   box_height_px=int(y1 - y0), box_width_px=int(x1 - x0),
                   n_z=int(nz), n_channels=len(chan_names))
        if extra_columns:
            rec.update(extra_columns)
        records.append(rec)

    return pd.DataFrame(records)


def read_box_provenance(path):
    """
    Read back the provenance dict embedded in an exported nucleus OME-TIFF.

    The JSON lives in the OME <Description> element and therefore comes back
    XML-escaped; this unescapes and parses it.
    """
    with tifffile.TiffFile(path) as tf:
        xml = tf.ome_metadata

    m = re.search(r"<Description>(.*?)</Description>", xml, re.S)
    if not m:
        return {}
    return json.loads(unescape(m.group(1), {"&quot;": '"', "&apos;": "'"}))
