"""antenna3d - nuclear F-actin antenna networks as one graph per nucleus.

    from antenna3d import describe_file
    print(describe_file("field.nd2"))          # ALWAYS do this first on a new dataset

    from antenna3d import run_folder, Params
    run_folder("data/", "results/", Params(), labels_dir="masks/")

The output is not an image. It is `graphs/<nucleus_uid>.graphml`, `antenna_edges.csv` and
`antenna_nuclei.csv`.

**Read `run_acceptance.py` before you quote a number.** See the README.
"""
from .acceptance import run_acceptance, summarise as summarise_acceptance
from .enhance import centerline_nms, frangi_3d, noise_null
from .graph import build_graph, edge_rows, nucleus_row
from .grid import Grid, sample_um
from .io import (Field, FieldRef, IngestError, describe_file, list_fields, nucleus_uid,
                 open_field, read_header)
from .optics import PSF, E1_100X_FROZEN, SamplingError, scales_for, theoretical
from .params import (Channels, DetectParams, GraphParams, Optics, Params, SegParams,
                     TraceParams, WIDEFIELD_100X)
from .pipeline import load_work, run_field, run_folder
from .preprocess import Crop, destripe_plane, flatten_in_mask, make_crop, noise_sd
from .reconnect import build_filaments, end_tangent, join_gaps, pair_at_junction
from .segment import (Nucleus, boundary_surface, labels_from_tiff, nuclei_from_labels,
                      segment_binned, segment_nuclei_cellpose)
from .trace import binarize, path_length_um, polyline_length_um, resample_to, trace
from .validate import field_figure, field_summary, trace_overlay

__version__ = "0.1.0"

__all__ = [
    "Channels", "Crop", "DetectParams", "E1_100X_FROZEN", "Field", "FieldRef", "Grid",
    "GraphParams", "IngestError", "Nucleus", "Optics", "PSF", "Params", "SamplingError",
    "SegParams", "TraceParams", "WIDEFIELD_100X", "binarize", "boundary_surface",
    "build_filaments", "build_graph", "centerline_nms", "describe_file", "destripe_plane",
    "edge_rows", "end_tangent", "field_figure", "field_summary", "flatten_in_mask", "frangi_3d",
    "join_gaps", "labels_from_tiff", "list_fields", "load_work", "make_crop", "noise_null",
    "noise_sd", "nuclei_from_labels", "nucleus_row", "nucleus_uid", "open_field",
    "pair_at_junction", "path_length_um", "polyline_length_um", "read_header", "resample_to", "run_acceptance",
    "run_field", "run_folder", "sample_um", "scales_for", "segment_binned",
    "segment_nuclei_cellpose",
    "summarise_acceptance", "theoretical", "trace", "trace_overlay",
]
