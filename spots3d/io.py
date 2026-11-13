"""Reading: one .nd2 channel at a time, and the nucleus3d label volume.

A TripleOP field is 27 x 8 x 2048 x 2048 uint16, 1.8 GB. Only the channels
that are detected in are read, each through the dask interface, so the other
five never enter memory.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np


class ND2Field:
    """One stage position of an .nd2, opened lazily. Use as a context manager."""

    def __init__(self, path, position: int = 0):
        import nd2
        self.path = str(path)
        self.position = int(position)
        self._fh = nd2.ND2File(self.path)
        sizes = dict(self._fh.sizes)
        if sizes.get("T", 1) != 1:
            raise ValueError(f"{self.path}: time series are not supported ({sizes})")
        n_pos = sizes.get("P", 1)
        if self.position >= n_pos:
            raise IndexError(f"{self.path} has {n_pos} position(s); asked for {position}")
        self.axes = list(sizes)
        self.shape = (sizes.get("Z", 1), sizes["Y"], sizes["X"])
        vs = self._fh.voxel_size()
        self.voxel_um = (float(vs.z), float(vs.y), float(vs.x))
        self.channel_names = [c.channel.name for c in self._fh.metadata.channels]

    @property
    def stem(self) -> str:
        return os.path.splitext(os.path.basename(self.path))[0]

    def channel(self, index: int, expected_name: str | None = None) -> np.ndarray:
        """The (Z, Y, X) volume of channel `index`, checked against `expected_name`."""
        if not 0 <= index < len(self.channel_names):
            raise IndexError(f"{self.path}: channel {index} does not exist "
                             f"({self.channel_names})")
        if expected_name is not None and self.channel_names[index] != expected_name:
            raise ValueError(
                f"{self.path}: channel {index} is {self.channel_names[index]!r}, "
                f"expected {expected_name!r}. Channels: {self.channel_names}")
        idx = []
        for ax in self.axes:
            idx.append({"P": self.position, "T": 0, "C": index}.get(ax, slice(None)))
        arr = np.asarray(self._fh.to_dask()[tuple(idx)].compute())
        if "Z" not in self.axes:
            arr = arr[None]
        if arr.shape != self.shape:
            raise ValueError(f"{self.path}: channel {index} has shape {arr.shape}, "
                             f"expected {self.shape}")
        return arr

    def close(self):
        self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def labels_path_for(labels_dir, stem: str, position: int) -> Path:
    """nucleus3d's handoff name, `<stem>_p<NN>.tif` (export_field_labels)."""
    return Path(labels_dir).expanduser() / f"{stem}_p{position:02d}.tif"


def read_labels(path, image_shape, image_voxel_um, tol_um: float = 1e-6) -> np.ndarray:
    """A nucleus3d label volume, accepted only if it is on the image's own grid.

    The label TIFF must carry its voxel size (nucleus3d writes it into the OME
    metadata) and both that and its shape must equal the image's. A label
    volume that is binned, cropped or uncalibrated is refused rather than
    rescaled: spot detection runs on the native grid, and a mask on any other
    grid puts every nucleus in the wrong place without raising anything.
    """
    import tifffile
    with tifffile.TiffFile(str(path)) as tf:
        labels = tf.asarray()
        meta = tf.ome_metadata or ""

    def _size(axis):
        m = re.search(rf'PhysicalSize{axis}="([0-9.eE+-]+)"', meta)
        return float(m.group(1)) if m else None

    voxel = (_size("Z"), _size("Y"), _size("X"))
    if None in voxel:
        raise ValueError(
            f"{path}: no OME PhysicalSize in the label file. spots3d only accepts "
            f"calibrated label volumes on the image grid, such as nucleus3d's "
            f"labels/<stem>_pNN.tif.")
    if labels.ndim != 3 or tuple(labels.shape) != tuple(image_shape):
        raise ValueError(f"{path}: label shape {labels.shape} != image shape "
                         f"{tuple(image_shape)}. They must be the same stack, unbinned.")
    if any(abs(a - b) > tol_um for a, b in zip(voxel, image_voxel_um)):
        raise ValueError(f"{path}: label voxel size {voxel} um != image voxel size "
                         f"{tuple(image_voxel_um)} um.")
    return labels
