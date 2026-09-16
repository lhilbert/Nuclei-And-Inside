# antenna3d

> This is the antenna pipeline's own reference, kept whole from the standalone package it was
> extracted from. It is the long form: the measurements behind every default, the bake-offs that
> chose each tool, and the caveats in full. The [repository README](../README.md) is the short
> form and the place to start.

Nuclear F-actin **antenna** networks from Nikon `.nd2` z-stacks, as **one graph per nucleus**.

Per-nucleus 3D ridge detection, centreline tracing, explicit crossing resolution and graph
construction, for HA-K-actin-stained (or equivalent F-actin-reporting) whole-mount nuclei.

**The output is not an image.** It is `graphs/<nucleus_uid>.graphml`, plus a tidy per-edge and
per-nucleus table that goes straight into downstream statistics.

Nucleus segmentation is deliberately not this package's job — bring masks from
[`nucleus3d`](../README.md#step-1--from-images-to-tables), the sibling package in this
repository, or from anywhere else. `nucleus_uid` is built from the same three parts on both
sides, so the two tables join on it; `scripts/run_combined.py` runs the two in order and does
the join for you.

---

## New to Git and GitHub?

See the repository root README — everything there applies here.

---

## Install

One install covers the whole toolbox — see [the repository README](../README.md#install).
`uv sync`, or `pip install -e .`, and both packages are importable.

Optional extras:

```bash
uv sync --extra parquet      # parquet beside the CSVs; the edge table is ~10x smaller
uv sync --extra cellpose     # ONLY if you want this package to segment nuclei itself
```

### Dependencies

`nd2`, `numpy`, `scipy`, `scikit-image`, `networkx`, `skan`, `pandas`, `tifffile`, `matplotlib`.

`cellpose` and `torch` are **optional** and are imported only if you ask this package to
segment. In this toolbox the masks come from `nucleus3d`, so neither is ever installed.

---

## Example data and code validation

Validated two ways.

**Synthetic phantoms with exactly known geometry** (`tests/`, no data needed): a straight
oblique filament of known length inside a spherical nucleus, at the measured widefield PSF.
`pytest` recovers its length to **1.03×** truth and traces it as **one** component.

**Real widefield 100× HA-K-actin whole-mount sphere data.** Two fields — one reporter-positive
and one **reporter-negative**, 37 nuclei, 11 750 edges — were run end to end through graphs,
tables and QC while this package was written, against masks produced elsewhere. Independently of
the research pipeline, on different masks and a different subsample, it lands on the same
per-edge geometry:

| | research pipeline, 363 nuclei | this package, 37 nuclei |
|---|---|---|
| median edge length, pos / neg | 0.757 / 0.739 µm | 0.763 / 0.746 µm |
| p90 edge length, pos / neg | 1.728 / 1.500 µm | 1.691 / 1.556 µm |
| median tortuosity, pos / neg | 1.035 / 1.022 | 1.035 / 1.009 |
| median intensity CV, pos / neg | 0.471 / 0.642 | 0.489 / 0.645 |
| control `separation` | 280 | 277 |
| acceptance `detect_neg` | 1.000 | 1.000 |

Be clear about what that does and does not establish. It shows the code path works on real
`.nd2` data and reproduces the control behaviour described under "Things that will bite you".
**It is not an independent validation of the method.** Every measured number quoted in this
README — the enhancer bake-off, the 4.35:1 detected cross-section, the acceptance-test verdict
on 363 nuclei across 13 fields — comes from the research pipeline this package was extracted
from, on that same dataset. Treat them as the provenance of the defaults, not as a claim that
the method has been validated on yours.

There is no public example dataset for this assay yet. `pytest` runs with no data at all; set
`ANTENNA3D_TEST_ND2` to one of your own `.nd2` files to exercise the reader as well.

---

## Use

**1. Check your channels and your voxel size first.** Two seconds, and it prevents the two most
expensive mistakes in the pipeline:

```bash
python -c "from antenna3d import describe_file; print(describe_file('yourfile.nd2'))"
```

```
130726_100xHAKactin-488_S5P-647_sphere_postfix_002_crop
  axes          ZCYX  {'Z': 123, 'C': 3, 'Y': 2280, 'X': 2588}
  voxel (z,y,x) 0.2 x 0.045998535 x 0.045998535 um
  slab          24.60 um over 123 planes
  positions     1
  objective     100x  NA 1.49  n 1.515  [['fluorescence', 'camera']]
  channels      (index -> OME name, emission)
      0  'DAPI'  460 nm
      1  'GFP'  525 nm
      2  'RFP'  632 nm
```

**2. Copy `scripts/run_antennas.py`**, edit the `SETTINGS` block — input folder, output folder, the two
channel **names**, the voxel size, the PSF, and `LABELS_DIR` — and run it:

```bash
python scripts/run_antennas.py
```

Put files in per-treatment subfolders of the input folder and a `condition` column appears,
named after the subfolder.

**3. Look at the QC figures** in `<output>/qc/` before you use the table. Not optional here —
see "Reading the QC figure".

**4. If you have a probe-absent control, run the acceptance test** — it is the only check here
that distinguishes a working detector from a confidently broken one, and on the data these
defaults come from it **fails**. See *Things that will bite you*.

```bash
python run_acceptance.py
```

To script something custom, import the package and skip the template:

```python
from antenna3d import Params, open_field, nuclei_from_labels, labels_from_tiff

params = Params()
labels, lgrid = labels_from_tiff("field_labels.tif")
nuclei = nuclei_from_labels(labels, lgrid, params)
with open_field("field.nd2", params.optics, params.channels) as field:
    ...
```

---

## Layout

```
antenna3d/
    grid.py         Grid: THE micron <-> voxel converter. Nothing else converts.
    optics.py       the PSF, and the ten scales derived from it; SamplingError
    params.py       every parameter, in five frozen dataclasses
    io.py           .nd2 -> Field (pixels + Grid + channels BY NAME); describe_file
    segment.py      bring-your-own labels | cellpose; boundary mesh; truncation covariate
    preprocess.py   destripe, flatten inside the mask, estimate the noise
    enhance.py      3D Frangi with gamma in noise units; the pure-noise null; across-ridge NMS
    trace.py        hysteresis -> NMS -> work grid -> skeletonize -> prune -> reject puncta
    reconnect.py    junction merging, crossing resolution, gap joining -> filaments
    graph.py        the per-nucleus networkx graph, and the two tidy tables
    validate.py     per-field QC figure, per-nucleus trace overlay, QC numbers
    acceptance.py   the probe-absent control test
    pipeline.py     batch driver: walks files/fields/nuclei, writes the output tree
scripts/run_antennas.py    TEMPLATE -- the only file you normally edit
scripts/run_acceptance.py  TEMPLATE -- the control test
```

The split is deliberate and is the same one `nucleus3d` uses: **`run_antennas.py` holds
parameters and no logic; the modules hold logic and no parameters.** One addition here —
`params.py` holds the defaults and the reason each one has the value it does, so a parameter's
justification lives next to the parameter rather than in a paper you no longer have.

---

## Output

```
<output>/
    antenna_nuclei.csv        one row per nucleus       (+ .parquet if pyarrow)
    antenna_edges.csv         one row per traced edge   (+ .parquet if pyarrow)
    field_summary.csv         one row per field: QC numbers
    run_parameters.json       every parameter used, plus the resolved scales
    graphs/                   one GraphML per nucleus            <- the deliverable
    qc/                       one validation figure per field
    qc/traces/                one overlay per nucleus
    work/                     bit-packed binary + mask, ~30 kB per nucleus
    skipped_files.csv         files not in this optical regime   (if any)
    failures.csv              fields that errored                (if any)
```

`work/` is what `run_acceptance.py` re-traces from. `KEEP_WORK = "full"` additionally stores the
flattened actin (~9 MB per nucleus), which the acceptance table's branch-point columns need;
`KEEP_WORK = False` keeps nothing and disables the acceptance test.

**Check `failures.csv` does not exist before you read any table.**

### Tracing a number back to its source

Every row carries `nucleus_uid`, `file`, `position`, `field_id`, `label` and `source_path`, plus
`condition` when files sit in per-treatment subfolders. `nucleus_uid` is the join key across
everything — and it is **the same format `nucleus3d` uses**, so:

```python
from combine import join_nuclei_and_antennas
both = join_nuclei_and_antennas("results/nuclei", "results/antennas")
```

That is a `merge(on="nucleus_uid")` with the checks that make it safe to trust — see
[Running both, and the join](#running-both-and-the-join). The bare merge works too, but returns
an empty frame rather than an error if the ids ever drift apart again.

The same id is the GraphML filename stem:

```python
import networkx as nx
G = nx.read_graphml("results/graphs/SetC_Control_004_crop_p00_nuc003.graphml")
G.graph["total_antenna_length_um"], G.graph["units"], G.graph["coordinate_order"]
```

**`nucleus_uid` maps forward only.** Build it from `(file stem, position, label)`; never parse
one back out of a path. A stem can itself contain the separator — in this project's own data a
field named `..._sphere_postfix__crop` has a doubled underscore in its stem — and inverting the
encoding once silently dropped **24 of 80** control nuclei, a fifth of the arm, with no error.
Every table carries `file`, `position` and `label` as their own columns so you never need to.

### Columns

**Per nucleus** (`antenna_nuclei.csv`). Sizes and denominators: `nucleus_volume_um3`,
`nucleus_imaged_volume_um3`, `nucleus_surface_um2`, `nucleus_surface_um2_true`,
`nucleus_surface_cap_fraction`, `nucleus_equivalent_diameter_um`, `axial_edge_fraction`.
Antennas: `n_antennas`, `n_segments`, `total_antenna_length_um`, `n_branch_points`,
`branching_ratio`, `branch_points_per_um`, `n_crossings_resolved`, `n_gap_joins`.
Ends, raw **and** cut-excluded in matched pairs: `n_roots`, `n_tips`, `n_cut_ends`,
`n_ends_total`, `n_components`, `n_components_uncut`, `length_in_uncut_components_um`.
Detection diagnostics: `noise_sd`, `null_level`, `separation`, `fg_fraction`, `threshold_low`,
`threshold_high`, `n_puncta_rejected`, `detected`.

**The three densities**, each naming its denominator:

| column | denominator | when to use it |
|---|---|---|
| `length_density_um_per_um3` | `nucleus_volume_um3` | the default |
| `length_density_um_per_um2_surface` | **`nucleus_surface_um2_true`** | when antennas are envelope-anchored |
| `antennas_per_um3` | `nucleus_volume_um3` | counts rather than length |

**Never report a total.** A nucleus here is 250–900 µm³ of imaged volume, so
`total_antenna_length_um` is mostly a statement about how big the nucleus was.

**Per edge** (`antenna_edges.csv`): `u`, `v`, `kind_u`, `kind_v`, `length`, `chord_length`,
`tortuosity`, `mean_intensity`, `intensity_cv`, `mean_width`, `mean_width_deconvolved`,
`filament_id`, `segment_id`, `orientation_{src,dst}_{z,y,x}`, `exit_angle_{src,dst}`,
`polyline`, plus the nucleus's ids and four denominators so a per-edge density needs no join.

**Node `kind`** is one of `root` (within `root_tolerance_um` of the envelope), `tip` (a free
end inside the nucleus), `branch`, `crossing` (a pass-through anchor where a crossing was
dissolved), and `cut` — an end against a **manufactured** face, which is neither a root nor a
tip. `cut` takes precedence over both, deliberately: counting a cut end as a tip inflates a
number that is then read as biology.

---

## Which tools were selected, and the numbers that chose them

Every choice was made by measurement on a phantom at this optics, not by reputation.

### Enhancement — Frangi, single scale, `gamma` in noise units

| enhancer | tracer | detected | dist. err (µm) | clDice | β₀ err | β₁ err |
|---|---|---|---|---|---|---|
| **frangi** | **NMS** | **6/6** | 0.096 | 0.839 | **0.17** | **2.00** |
| frangi | skeleton | **6/6** | 0.113 | 0.794 | 2.33 | 5.33 |
| sato | skeleton | 4/6 | 0.105 | 0.844 | 0.33 | 4.83 |
| meijering | NMS | 4/6 | 0.118 | 0.778 | 1.17 | 6.00 |

**Frangi is the only enhancer that detects in all six scenes.** Sato and Meijering fail outright
on the two faint ones (3.8σ peak SNR), and half of real nuclei sit below 4.8σ. Its advantage is
the **blobness term**: on a synthetic tube-plus-point it suppresses the point relative to the
tube by **47×**, against Sato's 5×.

**Single scale, matched to the PSF.** Filaments here measure 1.35–1.51× the PSF, so a scale
sweep would recover the PSF rather than the object. Width is still measured per segment and
flagged `width_is_psf_limited`.

**`gamma` in units of the measured noise sd, never the library default.** `skimage`'s
`frangi(gamma=None)` sets it per image to half the maximum Hessian norm, so the response is
incomparable between images: measured, a pure-noise image then scores **higher** than real data
(data/noise p99.9 = **0.73**). Fixed in noise units the ratio is **19.7**, and stable from 5 to
20 sd.

### Tracing — across-ridge NMS (Steger-style, 3D), then `skeletonize` + `skan`

NMS costs 2.4× the runtime and 0.017 µm of distance error — a tenth of the resolution — and buys
a **14× reduction in β₀ error** (0.17 against 2.33). β₀ error is the number of spurious breaks,
and a broken antenna is the failure mode that corrupts a graph: it becomes two antennas, with
two false tips and half the length each. On a one-voxel structure that error is invisible to
clDice and F1 — a broken filament scores 1.000 on both.

The reason it wins is measurable: after enhancement a filament's detected cross-section is
**4.35:1 in z**, against a PSF anisotropy of 4:1. Skeletonising that ribbon finds a medial axis that
wanders inside it; NMS takes the response maximum, which does not.

### The work grid is Nyquist-matched per axis, not isotropic in microns

`skeletonize` (Lee) has no spacing parameter — it decides topology on the index grid. What it
needs is not a grid isotropic in **microns** but structure isotropic in **voxels**. On a
0.092 µm isotropic grid the blurred filament is an 8 × 2 voxel ribbon, and skeletonising it gave
**7.86 µm** of wandering, branched centreline for a **6.55 µm** filament. Nyquist-matched per
axis, it is round, and the skeleton is a tube.

---

## Reading the QC figure

Five panels per field:

1. **Actin max projection** — the raw data.
2. **Nuclei (cyan) + traced centrelines (yellow)** — **the panel that matters.** Is what it
   traced a filament? Dense yellow scribbles inside a nucleus are not antennas; they are the
   rims of bright puncta. Nucleus label ids match the `label` column.
3. **Response against threshold** — the in-mask response p99.9 over the pure-noise null, per
   nucleus, with `low_k` and `high_k` marked. Nuclei left of the red line cannot seed a trace.
   **A high value does not mean filaments** — a control full of bright puncta scores higher
   than real signal.
4. **XZ side view** — how thick the slab is and how cut the nuclei are. `axial_edge_fraction`
   near 0 means masks close in z; near 1 means they are cut at their widest.
5. **Density against volume** — outliers are puncta-driven nuclei.

`qc/traces/<nucleus_uid>.png` gives the same for one nucleus, at full resolution. **Open a few.**
A length density cannot tell you the tracer is following the rim of a blob; this can.

`field_summary.csv` is the same information numerically, so you can sort by `separation_med` or
`n_puncta_rejected` to find the fields worth opening.

---

## Tuning

> **Read this first: no threshold is ever re-tuned per condition.** Every threshold here is a
> multiple of the response of **pure noise at the same parameters**, which is what makes one
> frozen value mean the same thing on a bright nucleus and a faint one. Tune on a control or a
> phantom, freeze the value, and use it for every condition. Tuning per condition is fitting the
> result, and it is undetectable in the output.

`UPPER_CASE` names are in `scripts/run_antennas.py`'s SETTINGS block. `lower_case` names are fields of
the `params.py` dataclasses; set them by replacing the block:

```python
from dataclasses import replace
PARAMS = PARAMS.with_(trace=replace(PARAMS.trace, max_angle_deg=25.0))
```

| Symptom | Parameter | Direction |
|---|---|---|
| dim filaments missed | `HIGH_K` | lower (2.0–3.0) |
| noise speckle traced | `HIGH_K`, `detect.min_object_voxels` | raise |
| one filament broken into pieces | `trace.max_gap_um` | raise — but see below |
| unrelated filaments joined together | `trace.max_angle_deg` | lower |
| puncta rims traced as loops | `trace.max_loop_perimeter_um`, `trace.max_loop_area_um2` | raise |
| short spurs everywhere | `MIN_BRANCH_UM` | raise, and report a curve |
| a diffuse background, not a flat one | `detect.background_sigma_um` | lower toward the object |
| the noise estimate tracks the signal | `detect.noise_estimator` | `first_difference_x` |
| skeleton wanders or self-branches | check the work grid is Nyquist-matched per axis | — |
| run is too slow | `MEASURE_WIDTH = False` | — |

**`max_gap_um` does not do what it looks like it does.** Measured: widening it does not reunite
broken filaments. The number of joins landing on the same object is flat at every gap width,
while the fraction landing on a *different* object rises from 83% to 97%. The gaps in a real
trace are mostly regions where nothing was detected at all, and nothing can be joined across a
region with no evidence in it.

Test on a few fields first — set `FIELDS` to two or three stems, or `LIMIT_NUCLEI = 5` — before
committing to a full run.

---

## Things that will bite you

**The reporter-negative control reports antennas, and on the data this was built from no
per-nucleus number is quotable as biology.** This is the headline caveat and it is not
hypothetical. On 13 fields of HA-K-actin whole-mount spheres with a **reporter-negative** arm —
same embryos, same session, same imaging, no actin construct, so the correct length density is
**zero** — **66 of 66 control nuclei reported antennas.** At `min_branch_um = 0.5`, with the
`first_difference_x` noise estimator this package ships: median density ratio **1.69**,
rank-biserial **+0.196**, p **0.013**. At 2.0 µm it is still 65 of 66. The traced objects are
geometrically indistinguishable between the arms (median edge length 0.757 vs 0.739 µm, median
tortuosity 1.035 vs 1.022).

Two things did **not** move that verdict, and both are worth knowing because both were expected
to. Repairing the noise estimator (`in_mask_mad` -> `first_difference_x`) improved every
discrimination statistic — ratio 1.56 -> 1.69, rank-biserial +0.178 -> +0.196, p 0.024 -> 0.013 —
and left `detect_neg` at 1.000. Re-segmenting the nuclei did not separate the conditions either.
**A detector that is wrong about what an object is does not become right when its denominator or
its masks improve.**

Three explanations were ruled out by measurement, not by argument: **not saturation** (5 × 10⁻⁶
of voxels clipped in the control, none in the reporter arm); **not a threshold merely set too
low** (the control's `separation` is **280** against the reporter arm's **14.5** — its puncta
clear a noise-referenced threshold by a wide margin, and raising the threshold removes the real
signal *first*); **not a mis-parameterised puncta filter** (it rejects *closed* loops, and these
rims are broken or fused to a neighbour, so they have free ends and survive — it caught 22 rims
out of well over a hundred).

The cause is the one thing single-scale ridge detection cannot help with: **a blob much wider
than the analysis scale presents its shoulder as a locally cylindrical surface, which *is* a
ridge at that scale.** Frangi's blobness term rejects a blob at the scale it is measured at; it
does not reject the flank of a blob four times wider. An independent pipeline over the same data,
failing by an entirely different route (an intensity threshold rather than a shape-selective 3D
Hessian), reached the same verdict on the same control.

**So: run `run_acceptance.py` on your own probe-absent control before you quote anything, and
open `qc/traces/` for the control nuclei.** If you have no such control, this package cannot
tell you whether its output is real, and neither can you. That is worth an acquisition.

**And read `acceptance_per_field.csv` beside the curve.** The curve pools nuclei across fields
and treats them as independent, which they are not. Median density across three fields of the
*same* condition spanned **0.0366 to 0.5581** on this assay — a 15× range — so a pooled effect
size is overconfident, and one computed from a single field per arm can come out with the wrong
**sign**. `detect_neg` is the robust part of the verdict; the effect size is not.

**The channel name is not the stain.** Check with `describe_file` on every new dataset. One
dataset in this project stores actin at index **2** and another at index **0** — same
instrument, five days apart. Resolving by index analyses the DNA channel as actin and never
errors. Channels here are resolved by **OME name**, and a missing one raises.

**A number in voxels is a bug.** Every length, area, volume and coordinate here is in microns,
and `Grid` is the only converter. Coordinates are `(z, y, x)`. Two specific instances: counting
skeleton voxels times one spacing underestimates an oblique filament by up to **1.73×**; and
summing raw voxel steps goes the other way and **overestimates by 1.33×**, because a digitised
line zigzags between the 26 available step directions. The reported `length` is measured over
smoothed polylines, which recovers a known phantom filament to 1.03×. `skeleton_length_um` is
the raw step sum and is a QC number only.

**The PSF is per experiment and per optical regime, and everything scale-like derives from it.**
The shipped defaults are widefield 100× at a 0.2302 µm lateral FWHM and **5.0 samples per FWHM**.
They do not transfer. At iSIM sampling (0.1265 µm FWHM, 1.94 samples per FWHM) the PSF-matched
lateral Hessian scale is **0.826 native pixels**, where a Gaussian derivative is not a derivative
— the kernel is narrower than the sample spacing — and `Optics.scales()` raises `SamplingError`
rather than running it. **Both axes are checked**: a filter that aliases in z aliases as badly as
one that aliases in x.

**The optical regime comes from the voxel size, never the filename.** One folder in this
project's data holds 100× and 40× stacks together. A file whose voxel size does not match the
declared `Optics` is refused; a file whose *name* lacks the declared marker is refused too. The
voxel size decides and the marker only catches a file misnamed at acquisition.

**One `.nd2` is many fields, and the position loop encloses the z loop.** A P axis read as a
plane axis turns 34 fields into one stack 34× too deep — a slab 34× too thick, unrelated nuclei
stacked on each other, and a plausible number for every one. `multi_series=True` is required
before any P axis is indexed, and the axis order is checked.

**Every observable is a density or a ratio and states its denominator.** And the denominator has
a trap in it: the marching-cubes surface of a nucleus that is cut in z includes a **manufactured
flat lid** where the mesh was closed. Measured at **42.5%** of the reported surface on a cut
nucleus, and **0.00 µm²** on one that is not. `nucleus_surface_um2_true` removes it and is what
`length_density_um_per_um2_surface` divides by; `nucleus_surface_um2` is emitted only so the two
can be compared.

**Masks get cut in z by the segmentation, not only by the slab, and it does not look wrong.**
`axial_edge_fraction` is the model-free covariate: 0 when the mask tapers away to nothing, 1 when
it is cut at its widest. On data that is **not** axially truncated by its acquisition — 20–34 µm
slabs against 7–10 µm nuclei — 2D-per-plane segmentation gave a median of **0.8764**, against
**0.0548** for a 3D flow field over the same fields. That is pure segmentation defect, and every
density computed against those masks was divided by a volume roughly 23% too small. Check panel
4 and the `axial_edge_fraction` column on any masks you bring — **including `nucleus3d`'s**,
which segments in 3D and scores 0.087 on a phantom but scored **0.90–1.00 on five of six** traced
nuclei of real 100× sphere data.

**The probe is a reagent, and the channel is not "actin".** HA-K-actin reports **F-actin**, so
the channel *is* the antenna. A genetically encoded actin chromobody reports actin **regardless
of assembly state**, so there an antenna is a sparse subset of a bright diffuse monomer pool.
**No observable transfers between the two**, and a detector validated on one cannot be validated
on the other even in principle. Whether the probe could have perturbed what it reports is decided
by one fact — its order relative to fixation. Applied post-fixation with a cross-linker present,
a filament-stabilising ligand is inert as used; applied live, it is not.

**A condition is a folder name, and nothing parses it.** In this project's own data `LatB_009`
and `LatA_004` are the **same compound** — latrunculin B, in embryos B and A — and there is no
latrunculin A in the dataset at all. Treating the trailing letter as a drug manufactures a drug
comparison out of one condition imaged twice, and treating nuclei as independent replicates when
the true n is *two embryos* overstates every p-value. Group on `condition` yourself, with the
design in front of you.

---

## Relationship to other pipelines

- **`nucleus3d`** — the sibling in this repository: 3D nucleus segmentation, per-nucleus
  quantification and table analysis. It writes a label volume per field into `<output>/labels`;
  point `LABELS_DIR` at that folder, or just run `scripts/run_combined.py`, which wires the
  handoff and the join. **The two do not accept the same data**: this package declares an
  optical regime and refuses a file whose voxel size disagrees, so `nucleus3d` will process
  stacks this one correctly declines.
- **AntEnnA** — the research repository this package was extracted from. It keeps the full
  nine-stage form with per-stage manifests, three experiment descriptors and the append-only
  decision record behind every number quoted above.

### Running both, and the join

```bash
python scripts/run_combined.py
```

Runs `nucleus3d`'s segmentation, hands its label volumes here, and joins the two tables:

```
results/
    nuclei/                        a full nucleus3d output tree
        nuclei_measurements.csv
        labels/<stem>_p00.tif      <- the handoff
        qc/
    antennas/                      a full antenna3d output tree
        antenna_nuclei.csv
        antenna_edges.csv
        graphs/<nucleus_uid>.graphml
        qc/traces/
    nuclei_and_antennas.csv        one row per nucleus both pipelines measured
```

The two trees are separate because both pipelines write a `field_summary.csv`, a
`run_parameters.json` and a `failures.csv`. Pointed at one folder they would overwrite each
other. Set `STOP_AFTER_SEGMENTATION = True` on a new dataset: this half is the long one, and
there is no point spending it on masks nobody has looked at.

`join_nuclei_and_antennas` **raises** rather than returning an empty frame when no id matches —
that is what a drifted id looks like, and it is otherwise indistinguishable from a dataset in
which nothing was detected. It also refuses a duplicated id, and disagreement on `label` or
`position` under a shared id. `file`, `condition` and `source_path` are suffixed `_nucleus` /
`_antenna` rather than collapsed: `file` is a basename with its extension on one side and a stem
on the other.

**Fewer joined rows than nucleus rows is correct, not a bug.** The two gate differently on
purpose: this package takes 100–4000 µm³ and drops nuclei touching an xy edge, because a nucleus
cut laterally has no usable denominator; `nucleus3d` takes anything over 15 µm³ and keeps border
nuclei, because a slab volume is not a nuclear volume. The join prints the reconciliation.
Measured on two real fields: 42 nucleus rows → 37 pass the volume gate → 22 survive the border
filter, and `n_antenna_only` is 0.

There is no requirement to use `run_combined.py` — `run_segmentation.py` and `run_antennas.py`
do the two halves independently; point `LABELS_DIR` at the former's `labels/` folder and the
join works the same.

**Both packages export a `SegParams` and a `describe_file`, and they are different objects.**
Different fields, different return shapes. In any script that touches both, import them
qualified — `import nucleus3d as n3`, `import antenna3d as a3` — as `run_combined.py` does. A
`from ... import *` binds one name to the other package's object and nothing errors.

### Pulling a newer antenna3d from AntEnnA

This package arrived here by `git subtree`, so its five extraction commits are part of this
repository's history rather than a squashed import. The prefix was flattened afterwards — the
package is at `antenna3d/`, its tests at `tests/antenna3d/`, its templates in `scripts/` — so a
later update is a re-split and a merge, not a plain `git subtree pull`:

```bash
# in the AntEnnA checkout
git subtree split --prefix=share/antenna3d -b antenna3d-export

# here
git fetch <antenna-remote> antenna3d-export
git subtree add --prefix=vendor/antenna3d FETCH_HEAD    # lands whole, then move the parts
```

Worth weighing against the alternative each time: upstream changes to a package that has been
adapted to this toolbox are rarely a clean apply, and reading the diff and porting it by hand is
often less work than reconciling one.
