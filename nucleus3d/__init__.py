"""
nucleus3d -- 3D nucleus segmentation and quantification for Nikon .nd2 stacks.

Typical use is through the template script (`run_analysis.py`); import the
package directly when you want to work field by field:

    from nucleus3d import load_field, describe_file, SegParams, segment_nuclei

    print(describe_file("field.nd2")["channels"])     # confirm the DNA channel
    stack = load_field("field.nd2", position=0)
    labels, props = segment_nuclei(stack.channel("DAPI"), stack.voxel_um)

Module map
    io        reading .nd2 files             -> Stack
    segment   segmentation                   -> labels, props
    quantify  per-nucleus intensities        -> DataFrame
    export    3D OME-TIFF crops per nucleus  -> files + index
    validate  per-field QC figures           -> figure + QC numbers
    pipeline  the batch driver tying it together
"""

__version__ = "1.0.0"

from .io import Stack, load_field, iter_fields, describe_file, n_positions
from .segment import SegParams, segment_nuclei, otsu_limit
from .quantify import quantify_nuclei, background_level, nucleus_uid
from .export import export_nucleus_boxes, read_box_provenance
from .validate import validation_figure, unsegmented_fraction
from .pipeline import run, process_field, find_nd2

__all__ = [
    "Stack", "load_field", "iter_fields", "describe_file", "n_positions",
    "SegParams", "segment_nuclei", "otsu_limit",
    "quantify_nuclei", "background_level", "nucleus_uid",
    "export_nucleus_boxes", "read_box_provenance",
    "validation_figure", "unsegmented_fraction",
    "run", "process_field", "find_nd2",
]
