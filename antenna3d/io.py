"""ND2 access.

Three invariants live here, and each of them produces wrong numbers that do not crash:

1. **Channel identity comes from the OME name, never a stored index.** One dataset in this
   project stores actin at index 2 and another at index 0 - same instrument, five days apart.
2. **The optical regime comes from the voxel size, never the filename.** One folder here holds
   100x and 40x stacks together. A filename marker is corroborating evidence only.
3. **One `.nd2` is many fields, and the position loop encloses the z loop.** A P axis read as a
   plane axis turns 34 fields into one stack 34x too deep, with unrelated nuclei on top of each
   other - and a plausible number for every one.
"""
from __future__ import annotations

import fnmatch
import re
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import nd2
import numpy as np

from .grid import Grid
from .params import Channels, Optics


class IngestError(RuntimeError):
    """Something about the data does not match what was declared. Not a bug."""


# --------------------------------------------------------------------------------------- headers

def read_header(path: str | Path) -> dict:
    """Everything about a file that does not require reading pixels."""
    path = Path(path)
    with nd2.ND2File(str(path)) as f:
        v = f.voxel_size()
        chans = [c.channel.name for c in f.metadata.channels]
        emis = [c.channel.emissionLambdaNm for c in f.metadata.channels]
        mic = f.metadata.channels[0].microscope
        sizes = dict(f.sizes)
        return {
            "file": str(path),
            "stem": path.stem,
            "axes": "".join(sizes.keys()),
            "sizes": {k: int(n) for k, n in sizes.items()},
            "dtype": str(f.dtype),
            "dxy_um": float(v.x),
            "dy_um": float(v.y),
            "dz_um": float(v.z),
            "channels_stored": chans,
            "emission_nm": [None if e is None else float(e) for e in emis],
            "na": float(mic.objectiveNumericalAperture or 0.0),
            "magnification": float(mic.objectiveMagnification or 0.0),
            "immersion_ri": float(mic.immersionRefractiveIndex or 0.0),
            "modality": str(mic.modalityFlags),
        }


def describe_file(path: str | Path) -> str:
    """What is actually in this file. **Run this on every new dataset before anything else.**

    It takes two seconds and prevents the most expensive mistake in the pipeline: a wrong
    channel produces a full, confident, meaningless table.
    """
    h = read_header(path)
    nz = h["sizes"].get("Z", 1)
    npos = h["sizes"].get("P", 1)
    lines = [
        f"{h['stem']}",
        f"  axes          {h['axes']}  {h['sizes']}",
        f"  dtype         {h['dtype']}",
        f"  voxel (z,y,x) {h['dz_um']:.6g} x {h['dy_um']:.8g} x {h['dxy_um']:.8g} um",
        f"  slab          {nz * h['dz_um']:.2f} um over {nz} planes",
        f"  positions     {npos}" + ("  <- one file, many fields" if npos > 1 else ""),
        f"  objective     {h['magnification']:g}x  NA {h['na']:g}  n {h['immersion_ri']:g}"
        f"  [{h['modality']}]",
        "  channels      (index -> OME name, emission)",
    ]
    for i, (name, em) in enumerate(zip(h["channels_stored"], h["emission_nm"])):
        lines.append(f"      {i}  {name!r}" + (f"  {em:g} nm" if em else ""))
    lines.append("")
    lines.append("  Set Channels(dna=..., actin=...) to the NAMES above, never the indices.")
    return "\n".join(lines)


# ------------------------------------------------------------------------------ what is declared

def resolve_channels(header: dict, channels: Channels) -> dict[str, int]:
    """Map role -> stored index by OME NAME. Never by index."""
    stored = header["channels_stored"]
    out: dict[str, int] = {}
    for role, name in (("dna", channels.dna), ("actin", channels.actin)):
        if name not in stored:
            raise IngestError(
                f"{header['stem']}: channel role {role!r} is declared as OME name {name!r}, "
                f"which this file does not have; it stores {stored}. "
                f"Run describe_file() on it.")
        if stored.count(name) > 1:
            raise IngestError(f"{header['stem']}: channel name {name!r} is ambiguous ({stored})")
        out[role] = stored.index(name)
    return out


def check_sampling(header: dict, optics: Optics, tol: float = 1e-6) -> None:
    """THE REGIME IS THE VOXEL SIZE. The filename marker is checked against it, not trusted."""
    stem = header["stem"]
    if abs(header["dxy_um"] - optics.dxy_um) > tol:
        raise IngestError(
            f"{stem}: lateral sampling is {header['dxy_um']!r} um but these Optics declare "
            f"{optics.dxy_um!r} um. Analysing it anyway would put every length, the PSF and all "
            f"ten derived scales at the wrong physical scale. Either this file belongs to a "
            f"different optical regime, or Optics needs updating.")
    if abs(header["dz_um"] - optics.dz_um) > tol:
        raise IngestError(
            f"{stem}: axial step is {header['dz_um']!r} um but these Optics declare "
            f"{optics.dz_um!r} um.")
    if optics.stem_marker is not None and optics.stem_marker not in stem:
        raise IngestError(
            f"{stem}: voxel size matches these Optics but the stem does not carry the declared "
            f"marker {optics.stem_marker!r}. A file misnamed at acquisition must fail loudly "
            f"rather than be analysed under the wrong label. Set Optics.stem_marker=None if "
            f"your filenames do not carry one.")
    if optics.na and header["na"] and abs(header["na"] - float(optics.na)) > 0.01:
        raise IngestError(
            f"{stem}: Optics declare NA {optics.na} but the file reports {header['na']}")


def n_positions(header: dict, multi_series: bool = False) -> int:
    """The number of fields in this file, refusing an undeclared position axis.

    A P axis read as a plane axis turns 34 fields into one stack 34x too deep - a slab 34x too
    thick, unrelated nuclei stacked on each other, and a plausible number for every one.
    """
    n_p = header["sizes"].get("P", 1)
    if n_p > 1 and not multi_series:
        raise IngestError(
            f"{header['stem']}: file has a position axis of length {n_p}, but multi_series was "
            f"not declared. Refusing to index frames. Pass multi_series=True if that is really "
            f"a stage-position loop.")
    if n_p > 1:
        axes = header["axes"]
        if "P" not in axes or "Z" not in axes or axes.index("P") > axes.index("Z"):
            raise IngestError(
                f"{header['stem']}: axis order {axes!r} does not put the position loop outside "
                f"the z loop, so `position * n_planes + z` would address the wrong frame.")
    return n_p


# ------------------------------------------------------------------------------------ the fields

@dataclass
class Field:
    """One field of view: pixels, the grid that maps them to microns, and channels by name."""
    path: Path
    stem: str
    position: int
    n_positions: int
    grid: Grid
    channel_index: dict[str, int]
    header: dict
    condition: str | None = None
    _f: object | None = None

    @property
    def field_id(self) -> str:
        return f"{self.stem}_p{self.position:02d}"

    @property
    def shape(self) -> tuple[int, int, int]:
        s = self.header["sizes"]
        return (int(s.get("Z", 1)), int(s["Y"]), int(s["X"]))

    @property
    def slab_um(self) -> float:
        return self.shape[0] * self.grid.dz_um

    @property
    def voxel_um(self) -> tuple[float, float, float]:
        return self.grid.spacing

    def plane(self, z: int, role: str = "actin") -> np.ndarray:
        """One (Y, X) plane of one channel."""
        f = self._require_open()
        n_z = f.sizes.get("Z", 1)
        fr = f.read_frame(self.position * n_z + int(z))
        ci = self.channel_index[role]
        return fr[ci] if f.sizes.get("C", 1) > 1 else fr

    def volume(self, role: str = "actin", z_slice=None, y_slice=None, x_slice=None) -> np.ndarray:
        """A (Z, Y, X) subvolume of one channel, read plane by plane."""
        f = self._require_open()
        n_z = f.sizes.get("Z", 1)
        zs = list(range(*(z_slice or slice(None)).indices(n_z)))
        ys = y_slice or slice(None)
        xs = x_slice or slice(None)
        first = self.plane(zs[0], role)[ys, xs]
        out = np.empty((len(zs),) + first.shape, first.dtype)
        out[0] = first
        for i, z in enumerate(zs[1:], start=1):
            out[i] = self.plane(z, role)[ys, xs]
        return out

    def _require_open(self):
        if self._f is None:
            raise IngestError(f"{self.field_id}: pixels requested outside `with open_field(...)`")
        return self._f


@dataclass(frozen=True)
class FieldRef:
    """A field named but not opened. What `list_fields` returns."""
    path: Path
    position: int
    n_positions: int
    condition: str | None

    @property
    def stem(self) -> str:
        return self.path.stem

    @property
    def field_id(self) -> str:
        return f"{self.stem}_p{self.position:02d}"


def list_files(root: str | Path, include=("*.nd2",), exclude=()) -> list[Path]:
    root = Path(root).expanduser()
    if not root.exists():
        raise IngestError(f"input folder does not exist: {root}")
    if root.is_file():
        return [root]
    out: list[Path] = []
    for pat in include:
        out.extend(sorted(root.rglob(pat)))
    keep = [p for p in out if not any(fnmatch.fnmatch(p.name, e) for e in exclude)]
    return sorted(set(keep))


def condition_of(path: Path, root: Path) -> str | None:
    """The immediate subfolder under the input root, or None when files sit at the top.

    Same convention as `nucleus3d`: put files in per-treatment subfolders and the column
    appears. **The trailing token of a folder name is not automatically a drug** - in this
    project's own data `LatB_009` and `LatA_004` are the SAME compound in two embryos, and
    treating them as two conditions manufactures a drug comparison out of one condition imaged
    twice. Nothing here parses the name; it is carried through verbatim for you to group on.
    """
    root = Path(root).expanduser()
    for candidate_path, candidate_root in ((path, root),
                                           (path.resolve(), root.resolve())):
        # The unresolved pair FIRST, and that ordering is the point: a symlinked input tree is
        # the normal way to stage a subset, and resolving the file jumps out of the tree to
        # wherever the data really lives - which silently returned None for every condition.
        try:
            rel = candidate_path.relative_to(candidate_root)
        except ValueError:
            continue
        return rel.parts[0] if len(rel.parts) > 1 else None
    return None


def list_fields(root: str | Path, optics: Optics, channels: Channels,
                include=("*.nd2",), exclude=(), multi_series: bool = False,
                ) -> tuple[list[FieldRef], list[dict]]:
    """Every field in the input folder, and one row per file that was skipped and why.

    A file is skipped rather than fatal when it simply belongs to another optical regime - one
    folder holding two magnifications is normal. Anything else raises.
    """
    refs: list[FieldRef] = []
    skipped: list[dict] = []
    root = Path(root).expanduser()
    for p in list_files(root, include, exclude):
        try:
            h = read_header(p)
        except Exception as e:                      # a file that is not readable is not silent
            skipped.append({"file": str(p), "reason": f"unreadable: {e}"})
            continue
        try:
            check_sampling(h, optics)
        except IngestError as e:
            skipped.append({"file": str(p), "reason": str(e)})
            continue
        resolve_channels(h, channels)               # raises: a wrong channel is never skippable
        n_p = n_positions(h, multi_series)
        cond = condition_of(p, root)
        refs.extend(FieldRef(p, i, n_p, cond) for i in range(n_p))
    return refs, skipped


@contextmanager
def open_field(ref: FieldRef | str | Path, optics: Optics, channels: Channels,
               position: int = 0, multi_series: bool = False):
    """Open one field for reading. Validates the header before any pixel is touched."""
    if isinstance(ref, (str, Path)):
        ref = FieldRef(Path(ref), position, 1, None)
    h = read_header(ref.path)
    check_sampling(h, optics)
    ci = resolve_channels(h, channels)
    n_p = n_positions(h, multi_series)
    if ref.position >= n_p:
        raise IngestError(f"{ref.path.name}: position {ref.position} of {n_p}")
    grid = Grid(dz_um=h["dz_um"], dy_um=h["dy_um"], dx_um=h["dxy_um"])
    f = nd2.ND2File(str(ref.path))
    try:
        yield Field(ref.path, ref.path.stem, ref.position, n_p, grid, ci, h, ref.condition, f)
    finally:
        f.close()


# ----------------------------------------------------------------------------------- identifiers

_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def nucleus_uid(stem: str, position: int, label: int) -> str:
    """`<stem>_p<NN>_nuc<NNN>` - the join key, and the same format `nucleus3d` uses.

    So `antenna_nuclei.csv.merge(nuclei_measurements.csv, on="nucleus_uid")` works once both
    pipelines have run over the same data.

    It is built FORWARD from (stem, position, label) and is never parsed back out of a path. An
    id that has to be decoded from a filename is ambiguous the moment a stem contains the
    separator - which happened here, and silently dropped a fifth of a control arm with no
    error. Every table carries `file`, `position` and `label` as their own columns instead.
    """
    return f"{_SAFE.sub('-', stem)}_p{int(position):02d}_nuc{int(label):03d}"
