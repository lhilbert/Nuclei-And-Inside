"""
Step 1: images to tables.

Reads Nikon .nd2 z-stacks, segments nuclei in 3D, quantifies each one, and
writes the measurement table -- optionally with a per-nucleus substack per
segmented object. Nothing in here imports from `nucleus3d.analysis`, and
that direction is the point: the measurement table is the only interface
between the two halves, so this half can run on a headless cluster node
with no plotting concerns, and the analysis half can run anywhere the table
is, without the image data.

    io          .nd2 access: Stack, load_field, iter_fields, describe_file
    segment     3D segmentation: SegParams, segment_nuclei
    quantify    per-nucleus measurements: quantify_nuclei, midplane_metrics
    export      per-nucleus OME-TIFF substacks, with reuse detection,
                and the field label volume antenna3d reads
    validate    per-field QC figures
    parallel    memory-aware multiprocessing of fields
    pipeline    the driver: run, process_field
"""

from .io import Stack, load_field, iter_fields, describe_file, n_positions
from .segment import SegParams, segment_nuclei, otsu_limit
from .quantify import (quantify_nuclei, background_level, nucleus_uid,
                       midplane_metrics, weighted_radius)
from .export import (export_nucleus_boxes, export_field_labels,
                     read_box_provenance, box_is_current, params_fingerprint)
from .validate import validation_figure, unsegmented_fraction
from .parallel import (run_parallel, plan_workers, estimate_peak_bytes,
                       field_geometry, available_memory_bytes, pool_available)
from .pipeline import run, process_field, find_nd2, code_fingerprint, field_cache_key

__all__ = [
    "Stack", "load_field", "iter_fields", "describe_file", "n_positions",
    "SegParams", "segment_nuclei", "otsu_limit",
    "quantify_nuclei", "background_level", "nucleus_uid",
    "midplane_metrics", "weighted_radius",
    "export_nucleus_boxes", "export_field_labels", "read_box_provenance",
    "box_is_current", "params_fingerprint",
    "validation_figure", "unsegmented_fraction",
    "run_parallel", "plan_workers", "estimate_peak_bytes", "field_geometry",
    "available_memory_bytes", "pool_available",
    "run", "process_field", "find_nd2", "code_fingerprint", "field_cache_key",
]
