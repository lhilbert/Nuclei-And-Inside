"""
Reading Nikon .nd2 stacks.

A file may hold several stage positions ("fields"). Each field is returned
as a `Stack`: one 3D multichannel volume plus the metadata needed to turn
pixel counts into microns.
"""

import os
from dataclasses import dataclass

import numpy as np
import nd2


@dataclass
class Stack:
    """One 3D multichannel field of view."""

    data: np.ndarray      # (Z, C, Y, X)
    channels: list        # channel names, e.g. ['DAPI', 'Cy5', 'mCherry']
    voxel_um: tuple       # (dz, dy, dx) in microns
    name: str             # identifier, e.g. 'SetC_Control_004_crop.nd2#p03'
    source_path: str      # full path of the .nd2 it came from
    position: int         # stage position index within that file

    def channel(self, which):
        """Return the (Z, Y, X) volume for a channel name or index."""
        idx = self.channels.index(which) if isinstance(which, str) else which
        return self.data[:, idx]

    @property
    def voxel_volume_um3(self):
        return float(np.prod(self.voxel_um))

    @property
    def z_extent_um(self):
        return self.data.shape[0] * self.voxel_um[0]


def describe_file(path):
    """
    Channel names, voxel size and position count, without reading pixels.

    Worth calling before any batch run: it is the cheapest way to confirm
    that the channel you are about to segment is the one you think it is.
    """
    with nd2.ND2File(path) as fh:
        return dict(
            path=path,
            name=os.path.basename(path),
            channels=[c.channel.name for c in fh.metadata.channels],
            emission_nm=[getattr(c.channel, "emissionLambdaNm", None)
                         for c in fh.metadata.channels],
            n_positions=int(fh.sizes.get("P", 1)),
            n_z=int(fh.sizes.get("Z", 1)),
            shape_yx=(int(fh.sizes["Y"]), int(fh.sizes["X"])),
            voxel_um=(fh.voxel_size().z, fh.voxel_size().y, fh.voxel_size().x),
        )


def n_positions(path):
    """Number of stage positions in an .nd2 file."""
    with nd2.ND2File(path) as fh:
        return int(fh.sizes.get("P", 1))


def load_field(path, position=0):
    """
    Read a single stage position into a Stack.

    Multi-position files are read through the dask interface so that only
    the requested position is pulled into memory. Reading the whole file
    and slicing afterwards costs the full file size in RAM, which for a
    24-position stack is the difference between ~0.1 GB and several GB.
    """
    with nd2.ND2File(path) as fh:
        names = [c.channel.name for c in fh.metadata.channels]
        vs = fh.voxel_size()
        npos = int(fh.sizes.get("P", 1))

        if position >= npos:
            raise IndexError(f"{os.path.basename(path)} has {npos} position(s); "
                             f"asked for index {position}")

        arr = (np.asarray(fh.to_dask()[position].compute()) if npos > 1
               else np.asarray(fh.asarray()))

    arr = np.squeeze(arr)
    if arr.ndim != 4:
        raise ValueError(f"expected 4D (Z,C,Y,X), got {arr.shape} from {path}")

    base = os.path.basename(path)
    name = f"{base}#p{position:02d}" if npos > 1 else base

    return Stack(data=arr, channels=names, voxel_um=(vs.z, vs.y, vs.x),
                 name=name, source_path=path, position=position)


def iter_fields(path, positions=None):
    """
    Yield every field in a file (or a chosen subset).

    `positions` may be a list of indices, or None for all of them.
    """
    npos = n_positions(path)
    for pos in (range(npos) if positions is None else positions):
        yield load_field(path, pos)
