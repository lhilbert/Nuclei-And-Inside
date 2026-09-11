"""
Export of individual nuclei as 3D OME-TIFF boxes.

One file per nucleus: the bounding box padded by a physical margin, with
every channel and every z-slice retained, and the segmentation mask
appended as a final channel so it travels with the pixels.
"""

import json
import os
import re
from xml.sax.saxutils import unescape

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops

from .quantify import nucleus_uid


def export_nucleus_boxes(stack, labels, outdir, pad_um=1.0, include_mask=True,
                         compression="zlib", extra_columns=None):
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

    Returns
    -------
    DataFrame mapping each nucleus to its file, sharing `nucleus_uid` with
    the measurement table so the two join directly.
    """
    os.makedirs(outdir, exist_ok=True)

    dz, dy, dx = stack.voxel_um
    pad_y, pad_x = int(round(pad_um / dy)), int(round(pad_um / dx))
    nz, _, ny, nx = stack.data.shape

    chan_names = list(stack.channels) + (["NucleusMask"] if include_mask else [])

    records = []
    for r in regionprops(labels):
        _, y0, x0, _, y1, x1 = r.bbox
        y0, y1 = max(0, y0 - pad_y), min(ny, y1 + pad_y)
        x0, x1 = max(0, x0 - pad_x), min(nx, x1 + pad_x)

        sub = stack.data[:, :, y0:y1, x0:x1]               # (Z, C, y, x)

        if include_mask:
            m = ((labels[:, y0:y1, x0:x1] == r.label)
                 .astype(stack.data.dtype) * np.iinfo(stack.data.dtype).max)
            sub = np.concatenate([sub, m[:, None]], axis=1)

        uid = nucleus_uid(stack, r.label)
        fname = f"{uid}.ome.tif".replace("#", "_")
        path = os.path.join(outdir, fname)

        provenance = {
            "nucleus_uid": uid,
            "source_file": os.path.basename(stack.source_path),
            "source_path": stack.source_path,
            "position": int(stack.position),
            "nucleus_label": int(r.label),
            "bbox_yx_in_source": [int(y0), int(x0), int(y1), int(x1)],
            "pad_um": pad_um,
            "voxel_um_zyx": [dz, dy, dx],
            "channels": chan_names,
            "note": ("z not cropped; source stack is a thin slab and nuclei "
                     "are truncated top and bottom"),
        }
        if extra_columns:
            provenance.update({k: str(v) for k, v in extra_columns.items()})

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

        rec = dict(nucleus_uid=uid, file=os.path.basename(stack.source_path),
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
