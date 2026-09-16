# Nuclei-And-Inside

3D nucleus segmentation and per-nucleus quantification of Nikon `.nd2`
z-stacks, plus the analysis of the resulting tables. A Python port of the
lab MATLAB pipeline [VisitorGene_PseudoTime](https://github.com/lhilbert/VisitorGene_PseudoTime).

The work is **two steps that meet at one file**:

```
 .nd2 z-stacks ──┐
                 │  step 1   scripts/run_segmentation.py      hours
                 └─────────> nuclei_measurements.csv  ──┐
                             (+ optional per-nucleus     │
                              substacks, see below)      │
                                                         │  step 2
                                                         │  scripts/run_analysis.py
                                                         └─────────> figures   seconds
```

Step 1 needs the image data and a lot of RAM; step 2 needs only the table
(except the image mosaics, which re-read crops). Splitting them means a
long segmentation run happens once, and the analysis can then be iterated
in seconds -- on a different machine if you like.

A second, optional pipeline sits alongside these: `antenna3d` traces nuclear
F-actin antennas and returns one graph per nucleus. It does not segment
nuclei — it reads step 1's label volumes and joins on `nucleus_uid`. It is
documented on its own in [`docs/antenna3d.md`](docs/antenna3d.md), and it is
not part of the path above; see [the second
pipeline](#the-second-pipeline-antennas).

## New to Git and GitHub?

You don't need prior experience to contribute.

- **No GitHub account yet?** [Create one](https://docs.github.com/en/get-started/start-your-journey/creating-an-account-on-github),
  then work through GitHub's
  [Hello World quickstart](https://docs.github.com/en/get-started/start-your-journey/hello-world)
  (~15 min, entirely in the browser).
- **Learning Git itself:** the Carpentries'
  [Version Control with Git](https://swcarpentry.github.io/git-novice/)
  lesson is written for researchers and assumes no background.
- **Prefer a graphical client?** We recommend
  [GitHub Desktop](https://docs.github.com/en/desktop) or the Source Control
  panel in VS Code.
- **Why we work this way:** see the version control chapter of
  [The Turing Way](https://book.the-turing-way.org/reproducible-research/vcs).

Questions are welcome — open an issue rather than getting stuck.

---

## Install

To manage Python dependencies for a single project (think: additional functions
needed only for this project), it is highly recommended that you use a virtual
environment manager.

Below are the commands for either `uv` or Conda. Both work — you only need one.
We recommend `uv`.

### Option A — uv (recommended)

Install `uv`: https://docs.astral.sh/uv/getting-started/installation/

`uv` creates and manages a project-local virtual environment in `.venv`:

```bash
git clone https://github.com/<org>/nucleus3d.git
cd nucleus3d
uv sync
```

`uv sync` installs the package in editable mode and records the resolved dependency versions in `uv.lock`, so everyone working on the project gets an identical environment.

You do not need to activate anything — prefix commands with `uv run`:

```bash
uv run python scripts/run_segmentation.py
```

If you prefer an activated shell (e.g. for an interactive interpreter):

```bash
source .venv/bin/activate      # macOS / Linux
.venv\Scripts\activate         # Windows
```

### Option B — Conda (fallback)

Use this if you need packages that are awkward to install via pip on your
platform, or if your group already standardises on Conda.

Install Miniforge: https://github.com/conda-forge/miniforge#install

We recommend Miniforge over Anaconda or Miniconda: it provides the same `conda`
command, comes preconfigured for the free conda-forge channel, and avoids the
Anaconda Terms of Service, which can require a paid licence for institutional
use.

```bash
git clone https://github.com/<org>/nucleus3d.git
cd nucleus3d
conda create -n nucleus3d python=3.11
conda activate nucleus3d
pip install -e .
```

Note that this path does not use `uv.lock`, so dependency versions are resolved
fresh at install time and may differ from those used during development.

### Dependencies

`nd2`, `numpy`, `scipy`, `scikit-image`, `pandas`, `tifffile`, `matplotlib`,
`networkx`, `skan`

`networkx` and `skan` belong to `antenna3d` alone. It also has two optional
extras, neither of which the normal path needs — see
[`docs/antenna3d.md`](docs/antenna3d.md#install).

### Check it worked

```bash
uv run pytest
```

---

## Example data and code validation

Validated on two datasets: vt-iSIM fixed zebrafish embryos (JF646-Hoechst) and cultured cells (DAPI, drug conditions).

You can download one of these data sets from the following, publicly shared Zenodo repository. The analysis should run fine on these data.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5242952.svg)](https://doi.org/10.5281/zenodo.5242952)



---

## How the code is organised

```
nucleus3d/
    core/          step 1: images -> tables
        io.py          .nd2 access
        segment.py     3D segmentation
        quantify.py    per-nucleus measurements
        export.py      per-nucleus substacks (optional output)
        validate.py    per-field QC figures
        parallel.py    memory-aware multiprocessing
        pipeline.py    the driver: run(), process_field()
    analysis/      step 2: tables -> figures
        features.py    all-feature PCA
        plots.py       point figures: scatter, slab diagnostic, PCA summary
        mosaic.py      image figures: gallery, mosaic
        style.py       shared palette and font sizes
antenna3d/         the second pipeline -- its own package; see docs/antenna3d.md
combine.py         the join between the two pipelines' tables
scripts/
    run_segmentation.py   step 1 template -- copy per experiment, edit, run
    run_analysis.py       step 2 template
    run_antennas.py       antenna3d template
    run_combined.py       both pipelines, then the join
    run_acceptance.py     the antenna probe-absent control test
```

`antenna3d` is a sibling package, not a third `nucleus3d` subpackage: it
meets `nucleus3d` at the label volume and at `nucleus_uid`, and `combine.py`
holds that join outside both, so neither imports the other.

**The import direction is the contract:** `analysis` imports from `core`,
and `core` never imports from `analysis`. That is what lets step 1 run on a
headless node with no interest in plotting, and step 2 run wherever the
table is. If you ever need a core module to import a plotting function, the
thing you want probably belongs in `analysis`.

Only a handful of names are re-exported at the top level -- the two entry
points and the figure functions. Everything else is reached through the
half it belongs to, so the import says which step a function is part of:

```python
from nucleus3d import SegParams, run              # step 1
from nucleus3d.analysis import midplane_scatter   # step 2
from nucleus3d.core.quantify import midplane_metrics   # internals
```

## Step 1 — from images to tables

```bash
cp scripts/run_segmentation.py my_experiment.py    # optional: keep a record
python scripts/run_segmentation.py
```

Edit the `SETTINGS` block and nothing else. The three settings that matter
most on a new dataset:

| setting | why |
|---|---|
| `INPUT_DIR` | folder of `.nd2`, searched recursively; one subfolder per treatment gives you a `condition` column for free |
| `DNA_CHANNEL` | **check this every time.** Folder names lie -- one example dataset is called `...JF646Hoechst...` and has no Hoechst channel at all. Run `describe_file()` first |
| `N_WORKERS` | `"auto"` sizes the pool from measured memory, not from the core count. See the parallel section |

What it writes into `OUTPUT_DIR`:

| file | contents |
|---|---|
| `nuclei_measurements.csv` | one row per nucleus, ~54 columns |
| `field_summary.csv` | one row per field: nuclei found, unsegmented fraction |
| `run_parameters.json` | every setting used, plus the worker plan |
| `qc/*.png` | per-field validation figures -- look at these before trusting the table |
| `labels/*.tif` | per-field label volumes (see below) |
| `field_cache/` | per-field tables for fast re-runs (see below) |
| `nucleus_boxes/` | per-nucleus substacks, **only if you ask for them** |
| `failures.csv` | fields that errored, if any |

### Label volumes

`labels/<stem>_p<NN>.tif` is one compressed integer volume per field,
carrying the voxel size in its OME metadata. On by default
(`SAVE_LABEL_VOLUMES`), because it is cheap and it is what the second
pipeline reads. It is written even for a field that segmented to nothing, so
that an **empty** label volume means "nothing here" and a **missing** one
means something broke.

One consequence worth knowing: **it disables the measurement cache.** The
cache stores the table, not the label image, so a cached field returns
before segmentation runs and there would be no labels to write. Set
`SAVE_LABEL_VOLUMES = False` to get cached re-runs of the table alone.

### Per-nucleus substacks (optional)

Off by default, because it is the disk-hungry output: roughly 1-3 MB per
nucleus, so ~1-2 GB for a dataset of 800. Turn it on in the settings block:

```python
SAVE_NUCLEUS_BOXES = True
BOX_PAD_UM = 1.0            # margin around the bounding box
BOX_INCLUDE_MASK = True     # segmentation mask as an extra channel
```

Where they land and what they are:

```
results/
    nucleus_boxes/
        SetC_Control_004_crop_p00_nuc001.ome.tif     <file>_p<position>_nuc<label>
        SetC_Control_004_crop_p00_nuc002.ome.tif
        ...
    nucleus_boxes_index.csv                          nucleus_uid -> file name
```

Each file is a 3D OME-TIFF holding **all channels** of the cropped box,
plus the segmentation mask as one more channel when
`BOX_INCLUDE_MASK = True`. Physical voxel size travels in the OME metadata,
so ImageJ/Fiji opens them at the right scale. The provenance -- source
file, stage position, label, and the bounding box in source coordinates --
is embedded in the OME description:

```python
from nucleus3d.core.export import read_box_provenance
read_box_provenance("results/nucleus_boxes/SetC_Control_004_crop_p00_nuc001.ome.tif")
```

`nucleus_boxes_index.csv` shares the `nucleus_uid` key with
`nuclei_measurements.csv`, so a substack and its measurements join
directly:

```python
import pandas as pd
m = pd.read_csv("results/nuclei_measurements.csv")
b = pd.read_csv("results/nucleus_boxes_index.csv")
joined = m.merge(b, on="nucleus_uid")        # adds the ome_tiff column
```

On a re-run, a substack whose provenance and mask hash still match what
would be written is left in place rather than rewritten. The check is
content-addressed, not a timestamp -- see the re-running section.

### Running fields in parallel

Set `N_WORKERS` in `scripts/run_segmentation.py` (or `n_workers=` on `run`): `1` for
serial, an integer, or `"auto"`.

**The limit is memory, not cores.** One 31 × 716 × 794 field of the example
data is ~106 MB as uint16, but segmenting it peaks at **2.8–3.0 GB**
resident (measured). With ~30 GB genuinely free that is six workers at the
default `MEMORY_FRACTION = 0.75`; with 10 GB free it is one. What limits the
pool is free RAM, not installed RAM, so `n_workers=cpu_count()` is still
not a safe default.

Most of that used to be one array. `core.segment._distance_map` treats nuclei
clipped at the slab faces as interior by replicating the first and last
slices `pad = clip(round(12 µm / dz), 1, 64)` planes deep, and
`distance_transform_edt` works in float64: at `dz = 0.1 µm` the pad
saturates at 64 planes each side, so a 31-plane stack was transformed as
159 planes, and peak hit **7.1–7.5 GB**.

`SegParams.blockwise_distance` (default `True`) computes that transform one
connected component at a time on small lateral crops instead. This is not
an approximation — for a voxel inside a maximal connected component the
nearest background voxel lies inside that component's own bounding box, so
a crop cannot hide it, and the code additionally checks at runtime that no
distance exceeds its voxel's distance to the crop edge, widening the crop
or falling back to the whole-volume transform if it does. Verified on two
fields: label images and measurement tables bit-identical, runtime
unchanged (38.6 vs 40.3 s, 36.4 vs 39.6 s), peak 3.03 vs 7.07 GB and 2.82
vs 7.50 GB. Set it to `False` to reproduce the old path, e.g. to check the
equality on your own data.

`"auto"` therefore sizes the pool in two steps, and prints both:

1. `estimate_peak_bytes` predicts peak RAM from the file geometry alone
   (metadata only, no pixels read).
2. The first field is then processed in an isolated subprocess and its
   **actual** peak RSS measured; the pool is sized from the measurement.
   The calibration field is a real field — its results are kept, nothing is
   processed twice.

```
memory: 39 GB installed, 9 GB available, budget 7 GB
predicted peak 3.1 GB/field -> 2 worker(s) (memory allows 2, cores allow 13)
calibrating on field 1 ...
measured peak 3.0 GB/field (predicted 3.1) in 39s -> 2 worker(s), limited by memory
```

Note `available`, not `installed`: the run above had 9 GB free on a 39 GB
machine, so it got two workers. Closing other applications buys more
workers than any setting here.

`MEMORY_FRACTION` (default 0.75) is the share of *currently available* RAM
the pool may claim, so the count depends on what else is running — close
your other work before a big run, or lower the fraction to stay usable.
`MAX_WORKERS` caps it outright. The chosen plan is recorded in
`run_parameters.json` under `parallel_plan`.

Two practical notes. Workers are spawned, so the entry point must stay
guarded by `if __name__ == "__main__":` — both templates in `scripts/` already are, but
a script of your own that calls `run(..., n_workers="auto")` at module level
will fail. And if the platform refuses to create a process pool (some
sandboxes and hardened containers deny the POSIX semaphore probe
`ProcessPoolExecutor` makes at construction), the run falls back to serial
with a printed warning rather than failing.

### Re-running: what gets reused

Two independent caches, because the two expensive things are not the same
thing. Segmentation is **~38 s per field**; exporting that field's nucleus
boxes is **~0.4 s**. Reusing boxes alone therefore saves almost nothing —
what makes a repeated or resumed run cheap is not re-segmenting.

**The field cache** (`use_cache=True`, on by default when `save_qc`,
`save_boxes` and `save_labels` are all off) writes each field's measurement
table under
`<outdir>/field_cache/` and reads it back instead of re-segmenting. Measured
on two fields: 85 s cold, 0.53 s warm, and the returned table is
*bit-identical* to the freshly computed one.

All three of those settings need the label image, which is not cached --
only the table is. `SAVE_LABEL_VOLUMES` defaults to **on**, so a plain run
does not use the field cache; turn it off when you want a cached re-run of
the table alone.

The cache key covers everything that can change the numbers: the source
file's path, size and mtime, the stage position, the DNA channel,
`min_blob_um3`, a fingerprint of every `SegParams` field, and a SHA-1 of
`core/segment.py`, `core/quantify.py` and `core/io.py`. Editing a metric therefore
invalidates every entry automatically — no manual cache clearing, and no
risk of a stale table outliving the code that made it. It is bypassed when
`save_qc` or `save_boxes` is set, since both need the label image and only
the table is cached.

Getting bit-identity out of a CSV took two non-default settings, both in
`pipeline._write_cache` / `_read_cache`: `float_format="%.17g"` on write
(the default repr is shorter than float64 needs) and
`float_precision="round_trip"` on read (pandas' default C parser is fast
but not correctly rounded, and silently introduced differences of ~1e-14).
A sidecar `.schema.json` restores the column dtypes, which CSV cannot
carry — without it a column of whole-valued floats returns as `int64`.

**Box reuse** (`reuse_boxes=True`) leaves an exported OME-TIFF alone when
its embedded provenance shows it would be written identically: same source
file, crop, padding, channels, voxel size, `SegParams` fingerprint, and the
same mask content by SHA-1. The mask hash is what makes this safe —
parameters alone are not enough, since a box whose pixels no longer match
the current segmentation is worse than no box at all. Each row of the
returned index carries `reused` and `reuse_reason`, and the reason names
the first field that differed:

```
second pass : 8 boxes, reused=8, reasons=['match']
changed params -> reused=0, reasons=['seg_params_sha1 differs']
changed pad    -> reused=0, reasons=['bbox_yx_in_source differs']
```

### Reading the QC figure

Five panels per field:

1. **DNA max projection** — the raw data
2. **Segmentation** — outlines with label IDs matching the `label` column
3. **Unsegmented DNA** — DNA signal assigned to no nucleus, in cyan. This
   is the missed-nucleus check; a whole cyan nucleus means a real miss
4. **XZ side view** — how thick the slab is and how clipped the nuclei are
5. **Size vs solidity** — outliers are merges, fragments or debris

`field_summary.csv` gives the same information numerically, so you can
sort by `unsegmented_fraction` to find the fields worth opening.

---

### Tuning

| Symptom | Parameter | Direction |
|---|---|---|
| dim nuclei missed | `thresh_factor` | lower (0.4–0.6) |
| one nucleus split into pieces | `seed_depth_um` | raise |
| two nuclei merged into one | `seed_depth_um` | lower |
| boundaries inside the chromatin edge | `refine_factor` | lower |
| boundaries bleeding into background | `refine_factor` | raise |
| debris counted as nuclei | `min_area_um2`, `min_solidity` | raise |
| condensed/mitotic figures rejected | `min_solidity` | lower (~0.3) |

Test on a few fields first — set `POSITIONS = [0, 5, 10]` — before committing to a full run.

---

## Step 2 — from tables to figures

```bash
python scripts/run_analysis.py
```

Point `RESULTS_DIR` at what step 1 wrote. Everything here reads the table
only, so it runs in seconds and can be re-run as often as you like -- with
one exception: the mosaics re-read mid-plane crops from the `.nd2` files,
so set `MAKE_MOSAICS = False` if the image data is not reachable.

Every figure function drops two groups of nuclei by default, and both
matter more than they sound -- nuclei whose widest plane is a slab face
(`mid_at_z_border`) and nuclei touching an image edge in xy
(`touches_xy_border`). On the example data that leaves 241 of 797. The
`feature_pca` fit applies the same two filters, so the PCA and the scatter
describe the same population.

### Check the slab first: `zclip_diagnostics`

```python
from nucleus3d import zclip_diagnostics
zclip_diagnostics(table, "results/zclip_diagnostics.png")
```

Three panels: where each nucleus's widest plane sits in the stack, the
percentage of each condition clipped at a slab face, and whether being
clipped biases the measurements. It excludes xy-border nuclei by default,
so it measures the z effect alone rather than a mixture of the two
truncations. On the full example dataset (196 fields, 475 xy-interior
nuclei of 797):

| | |
|---|---|
| widest plane at a slab face | 234 of 475 interior nuclei (49%) |
| clipping rate by condition | Flavopiridol 44%, Triptolide 46%, Embryo sphere 47%, Control 57% |
| clipped ÷ unclipped median | mid-plane area **0.92** (p = 7e-04), CV 1.05 (p = 0.12), solidity 1.00 (p = 0.19), radial index 1.00 (p = 0.12) |

**Only mid-plane area is biased, and modestly.** Clipped nuclei measure 8%
smaller, which is what a grazing cut through the nuclear cap should do. CV,
solidity and the radial index show no detectable bias — all three are
ratios normalised by the mask itself, so a partial section rescales
numerator and denominator together.

Note how much the xy filter matters to this estimate: the same area
comparison run over all 797 nuclei gives 0.82 at p = 4e-07, because
laterally cut nuclei are both smaller *and* somewhat more likely to be
z-clipped. Measuring one truncation without excluding the other overstates
it by more than a factor of two.

Two warnings from getting this wrong first time. Both conclusions above are
*reversals* of what the first 51 fields (166 nuclei) suggested: that subset
put the CV ratio at 0.92 with p = 0.021, and made the clipping rate look
like it varied 2.3-fold between conditions (26% to 60%) rather than the
1.2-fold (45% to 55%) seen across all 196 fields. A per-condition rate
estimated from two or three files is dominated by which part of each
specimen happened to be in focus, and an apparently significant bias at
n = 166 was not reproducible at n = 797. Run `zclip_diagnostics` on the
whole dataset, not a pilot subset, before deciding what to filter.

---

### The condition figure

```python
from nucleus3d.analysis import midplane_scatter
midplane_scatter(table, "results/midplane_scatter.png")
```

One panel per condition, every nucleus a point, all other nuclei in grey
behind for context, and a cross at each condition median. CV is on the
x-axis of every row; the rows differ in what is plotted against it —
solidity, the radial index, and mid-plane area. Pass `metrics=` for other
pairs. `scripts/run_analysis.py` writes it automatically.

### Looking at the nuclei behind the points: `nucleus_gallery`

```python
from nucleus3d.analysis import nucleus_gallery
nucleus_gallery(table, "results/nucleus_gallery.png",
                x="mid_dna_cv_corr", y="mid_solidity", grid=(3, 4))
```

Tiles the scatter plane into a grid, picks the nucleus nearest each cell
centre, and shows its mid-plane xy section beside a scatter marking which
nuclei were picked — laid out in the same geometry as the plane, so the
gallery reads as the scatter itself. Crops come from the raw `.nd2` on
demand (`source_path`, `position` and `mid_z` from the table), each field
opened once however many nuclei it contributes, so the image files must be
reachable. Empty cells are labelled as such: a shape/texture combination
that does not occur is as informative as a tile.

Every tile is the same physical size (`window_um`, default 16 µm) and is
zero-padded if the window runs off the field edge — rare now that
xy-border nuclei are excluded by default. Nuclei are picked without
replacement, since the cell tolerance otherwise lets one nucleus win two
adjacent cells and appear as two examples. Contrast is stretched per tile,
which is right for judging chromatin texture and wrong for comparing
brightness between tiles. Neighbouring nuclei often appear at a tile's
edge; only the centred one is the subject.

On the example data the gallery makes one thing plain immediately: the
high-CV region (CV ≳ 0.4) is **morphologically heterogeneous**. It holds
both mitotic figures — condensed, individually resolved chromosomes — and
interphase nuclei with bright chromatin puncta, and solidity separates the
two only partially: mitotic figures appear at solidity 0.89 *and* 0.97,
punctate interphase nuclei at 0.88 *and* 0.98. So a high CV on its own does
not identify a cell state, and neither does the CV–solidity pair. That is
worth knowing before reading a shift in median CV as a chromatin-texture
effect: it may be a change in the mitotic fraction instead. Look at the
tiles.

(An earlier reading of a 12-nucleus gallery from 51 fields had solidity
cleanly separating the two morphologies. Drawing the gallery from all 196
fields contradicted it. The lesson is the same one as in "What the example
data taught us": estimate this kind of thing from the whole dataset, never
from a pilot subset.)

### The continuum, not examples: `nucleus_mosaic`

```python
from nucleus3d.analysis import nucleus_mosaic
# the plane itself, densely tiled
nucleus_mosaic(table, "results/nucleus_mosaic.png", grid=(6, 9))
# or one-dimensional, ranked by chromatin contrast
nucleus_mosaic(table, "results/nucleus_pseudotime.png", grid=(8, 12),
               sort_by="mid_dna_cv_corr")
```

Where `nucleus_gallery` gives a dozen captioned examples with a scatter
beside them, this composites many mid-plane sections into a single image
whose axes are the metrics — the layout to use when the plane appears to
order nuclei along something monotonic. On the example data the ramp is
legible: uniform, low-contrast interphase chromatin at low CV, then
increasingly punctate nuclei, then condensed mitotic figures past CV ≈ 0.5.
Whether that is cell-cycle progression is a hypothesis the mosaic makes
visible, not one it tests — the DNA channel alone cannot date a nucleus,
and a mitotic classifier or a cell-cycle marker would be needed to say so.

Three implementation choices are worth knowing, because each one is a
trade:

1. `binning="quantile"` (default) puts cell edges at equal-count quantiles, so the mosaic is dense where the nuclei are — a linear 6 × 9 grid fills 30 of 54 cells on this data, a quantile grid 52 of 54. The price is that axis *spacing* is rank, not value; the tick labels carry the real metric values at the cell edges. `binning="linear"` keeps the honest plane and the empty cells.
2. The tiles are pasted into one array and drawn with a single `imshow`, so cell boundaries and grid lines coincide by construction. Note the consequence: whenever cell edges are unevenly spaced (quantile or `sort_by` mode) the axes must be in cell units, since stretching the canvas linearly over a value range would put the grid lines off the tiles they mark.
3. `window_um=None` (default) sizes every tile from the largest nucleus drawn (1.2 × its mid-plane equivalent diameter, 22 µm here). A fixed 16 µm window clips the biggest 1.2% of nuclei at the *tile* edge, which looks exactly like the xy-border filter having failed when it is only framing.

That last point answers a question that will come up: nuclei in a tile may
appear to run off the edge. Check `touches_xy_border` before suspecting the
filter — a nucleus 11 µm across in a 16 µm tile fills two thirds of the
frame and its neighbours intrude at the margins, which is a framing effect,
not a truncated mask.

### All features at once: `feature_pca`

```python
from nucleus3d.analysis import feature_pca, pca_summary, nucleus_mosaic
res = feature_pca(table)                       # z-scored, 6 components
pca_summary(res, "results/pca_summary.png")
nucleus_mosaic(res["scores"], "results/pca_mosaic.png", x="PC1", y="PC2")
```

`feature_pca` is an SVD of the z-scored feature matrix — identical to
`sklearn.decomposition.PCA` up to component sign, which is pinned so that
each component's largest loading is positive (otherwise the same data can
produce a mirrored figure on another machine). No new dependency; see
`pca.py` for why.

**"All features" needs two qualifications, and both change the answer.**
Position columns and per-field backgrounds are excluded: stage coordinates
would invent a component that maps the field layout, and a `*_bg` value is
one number per field shared by every nucleus in it, so it aliases onto
condition. Then near-duplicates are pruned at |r| ≥ 0.98, keeping one
representative per group by a declared priority. This is not cosmetic:
`n_voxels` and `volume_um3` are the same measurement (r = 1.0000), as are
`max_area_um2` and `max_area_um2`, and each channel contributes median,
mean, p90, std and background-corrected twins above r = 0.98. Unpruned, a
quantity written down five times gets five votes and PC1 becomes a
bookkeeping artefact. On the example data 41 numeric columns reduce to
**22 features**; `res["dropped"]` names every removal and its reason.

On the 241-nucleus analysis population the first three components carry
41%, 19% and 12% of the variance (72% cumulative):

| | strongest positive loadings | strongest negative | reading |
|---|---|---|---|
| PC1 | `DAPI_std`, `dna_cv`, `DAPI_p90`, `dna_p90_p10` | `mCherry_std`, `mCherry_mean_corr`, `Cy5_enrichment`, `solidity` | chromatin contrast against Pol II channel signal and smooth outline |
| PC2 | all three channels' intensity, `volume_um3` | `solidity`, `dna_radial_norm` | overall brightness and size |
| PC3 | `dna_radial_um`, `volume_um3`, `DAPI_integrated_corr` | — | size and radial spread |

PC1 recovers the axis the hand-picked CV/solidity plane was tracking
(r = +0.75 with `mid_dna_cv_corr`, −0.52 with `mid_solidity`), and it is
independent of mid-plane area (r = +0.04).

**The caveat that matters here.** 48.6% of the variance *along PC1* is the
culture-versus-embryo split — a different specimen imaged in a different
session, so intensity-derived features carry acquisition differences as
well as biology. PC1 is therefore half specimen and half chromatin state on
this dataset, and a PC1 difference between a culture condition and the
embryos is not interpretable. Within the culture set alone the chromatin
reading holds up (r = +0.80 with `mid_dna_cv_corr`, n = 172). Fit the PCA
per specimen type, or drop the intensity features, before comparing across
them. PC2 and PC3 are nearly free of the split (1.5% and 4.9%).

## The measurement table, column by column

Reference for `nuclei_measurements.csv`: one row per nucleus, ~54 columns.
The families are geometry (from `core/segment.py`), per-channel intensity,
chromatin texture, radial distribution, and the mid-plane read-out. Two
columns are flags you should almost always filter on —
`mid_at_z_border` and `touches_xy_border`; see "What the example data
taught us" for why.

### Tracing a number back to its source

Every measurement row carries `nucleus_uid`, `file`, `position`, `field`, `source_path` and `label`, plus `condition` when files sit in per-treatment subfolders. The `nucleus_uid` is the join key across everything:

```python
import pandas as pd
m = pd.read_csv("results/nuclei_measurements.csv")
b = pd.read_csv("results/nucleus_boxes_index.csv")
both = m.merge(b, on="nucleus_uid")       # table row -> OME-TIFF on disk
```

The same ID is the filename stem of the exported crop, and it is embedded in the OME-TIFF header along with the bounding box in source coordinates:

```python
from nucleus3d import read_box_provenance
read_box_provenance("results/nucleus_boxes/SetC_Control_004_crop_p00_nuc003.ome.tif")
```

### Intensity columns

For each channel `<C>`: `<C>_median`, `<C>_mean`, `<C>_p90`, `<C>_std`
(raw), `<C>_bg` (field background), `<C>_median_corr`, `<C>_mean_corr`,
`<C>_integrated_corr`, `<C>_enrichment`.

**Use `<C>_integrated_corr` when comparing nuclei of different sizes** —
it is the background-subtracted sum over the nucleus in intensity·µm³, so it scales with total molecule number rather than concentration.

Chromatin texture on the DNA channel: `dna_cv`, `dna_p90_p50`,
`dna_p90_p10`, `dna_top10_fraction`. Mitotic figures score high on all four (`dna_cv` ≈ 0.4 versus ≈ 0.2 for interphase).

### Radial distribution columns

`dna_radial_um` is the intensity-weighted mean radius: for every voxel in
the nucleus, its distance from the mask centre weighted by the
background-subtracted DNA intensity there, summed and normalised by total
intensity. `dna_radial_norm` divides that by the *unweighted* mean radius
of the same mask, which is the form to compare between nuclei — it is
dimensionless and independent of nucleus size:

| `dna_radial_norm` | meaning |
|---|---|
| ≈ 1.0 | DNA spread through the mask like the mask itself |
| > 1.0 | signal pushed toward the nuclear rim |
| < 1.0 | signal concentrated centrally |

Negative weights (a few dim voxels after background subtraction) are
clipped to zero; without that they would pull the intensity centroid the
wrong way.

### Mid-plane columns — the 2D read-out

`mid_*` columns describe the single z-plane where each nucleus is widest,
rather than an average over the stack: `mid_z`, `max_area_um2`,
`mid_solidity`, `mid_dna_cv`, `mid_dna_cv_corr`, `mid_dna_radial_um`,
`mid_dna_radial_norm`, `mid_at_z_border`.

Why one plane instead of the whole stack: in a thin slab the top and
bottom planes are grazing cuts through the nuclear cap — small, ragged,
and dominated by partial-volume effects. A solidity or CV averaged over z
mixes real chromatin structure with how much of the nucleus the slab
happened to catch. The widest plane is the nearest reproducible
approximation to an equatorial section.

**`mid_at_z_border` is the column to filter on.** It is `True` when the
widest plane is the first or last plane of the stack, which means the
nucleus is cut off in z and its true widest section was never imaged. In
the bundled example data (3.1 µm slabs, nuclei ~10 µm across) this is a
large fraction of all nuclei, and `midplane_scatter` drops them by
default. Solidity and CV for those nuclei describe a grazing cut, not a
mid-nuclear section.

`mid_dna_cv` and `mid_dna_cv_corr` differ only in whether the camera
offset is subtracted from the mean. The offset inflates the mean but not
the standard deviation, so the raw CV is systematically damped and depends
on detector settings; the corrected one estimates chromatin contrast and
is what to compare between experiments.

## What the example data taught us

These are findings about this dataset and this imaging geometry, not
properties of the code -- but each one cost a day to establish, and each
one changes how a number should be read.

### Nuclei cut at the image edge

`touches_xy_border` (from `core/segment.py`, true when a mask reaches x = 0,
y = 0, or the far edge) flags nuclei truncated laterally. All three figure
functions now drop them by default (`exclude_xy_border=True`).

This is the larger of the two truncations, and on the example data it is
not subtle — border nuclei against interior ones, median ratios over all
797 nuclei:

| metric | ratio | p |
|---|---|---|
| mid-plane area | **0.61** | 5e-27 |
| nuclear volume | **0.58** | 3e-24 |
| CV | 0.95 | 0.001 |
| solidity | 0.99 | 1e-08 |
| radial index | 1.01 | 3e-07 |

A laterally cut nucleus is a fragment, so every shape and texture read-out
on it describes the fragment, not the nucleus — and unlike the z case, even
the mask-normalised metrics shift measurably. 322 of 797 nuclei (40%) are
affected, which is why they appeared in nearly every gallery tile before
this filter existed.

Requiring both filters — interior in xy *and* not z-clipped — leaves **241
of 797 nuclei** (Control 59, Flavopiridol 35, Triptolide 78, Embryo sphere
69). That is the population the figures use, and it is the honest cost of
asking for a complete equatorial section in a 3.1 µm slab of a 46 × 52 µm
field.

### Which features survive a change of imaging session

Refitting the PCA on the culture set alone (172 nuclei, 25 features — the
redundancy pruning is data-dependent, so three intensity columns that were
duplicates across specimens survive here) splits the components cleanly by
what they are made of, and the result is a caution about DNA content:

| component | built from | variance by imaging set | by condition |
|---|---|---|---|
| PC1 (30%) | DNA-channel contrast: `DAPI_std`, `dna_cv`, `dna_p90_p10` | **0.0%** | 23.7% |
| PC2 (29%) | `Cy5_*` intensity, `mCherry_*`, `volume_um3` | 1.5% | 20.9% |
| PC3 (11%) | `DAPI_median`, `DAPI_integrated_corr`, solidity | **13.6%** | 1.2% |

The direct test is the control nuclei, which are the same biology imaged in
two sessions (`SetC` n = 34, `SetD` n = 25):

| | SetC | SetD | p |
|---|---|---|---|
| `mid_dna_cv_corr` | 0.274 | 0.283 | 0.9 |
| `DAPI_integrated_corr` | 4.3e5 | 9.0e5 | 2e-07 |

**Intensity ratios travel between sessions; absolute intensities do not.**
CV is a within-nucleus ratio and is indistinguishable across the two
sessions, while integrated DAPI — the total DNA amount — differs
two-fold on identical biology. That kills the obvious cell-cycle approach:
DNA content is exactly the quantity that would separate G1 from G2, and it
is the one this data cannot compare across sessions. Staging the cycle
would need intensities normalised per session (e.g. each field or set
rescaled by its own control median) before `DAPI_integrated_corr` means
anything, or a cycle marker in one of the free channels.

#### Why DNA content cannot stage the cell cycle here (tested, negative)

The obvious next move is to normalise intensities per session and read
`DAPI_integrated_corr` as DNA content, expecting G1 and G2 as two clouds a
factor two apart. It does not work, for three reasons measured on this
data — worth recording so nobody spends a week on it:

1. **Every nucleus is truncated in z.** The imaged mask is 2.94 ± 0.11 µm thick, pinned at the 3.1 µm slab, while an ellipsoid of the same equatorial area and volume needs 4.45 µm (median; >3.1 µm for 100% of nuclei, and 96% of masks reach within 0.2 µm of both slab faces). `volume_um3` is therefore mid-plane area × slab thickness — which is why it correlates at r = 0.9946 with `max_area_um2` — and the integrated signal is DNA in a slice, not DNA in a nucleus. `mid_at_z_border` was never a containment test; it only says the widest plane is interior.
2. **Field-to-field brightness swamps the signal.** Among control nuclei, field explains R² = 0.887 of the variance in log2 DNA amount against 0.602 expected from the group count alone (permutation p < 0.001; restricted to fields with ≥ 2 controls, 0.816 vs 0.441, p = 0.0005). Per-field medians span 5.9× after per-session normalisation and 3.8× after dividing by each field's own background, against 2.0× for a G1→G2 doubling. The mechanism is visible: nuclear brightness follows its field's background (Spearman ρ = +0.57, p = 3e-06), and background is empty medium carrying no DNA.
3. **The apparent two clouds are fields.** A pooled KDE does show two modes 2.1× apart, which looks exactly like a doubling — but the low mode has a prominence of 0.012 of peak density, it is absent in `SetC` and marginal in `SetD` (prominence 0.09, n = 25), and centring each nucleus on its own field's mean removes it entirely, collapsing the spread from 1.56 to 0.60 in log2. The split also runs on concentration (1.97×) rather than area (1.44×), which is the signature of exposure, not of replication.

What would be needed: an intensity standard imaged with the sample, or a
cell-cycle marker in one of the two free channels. Absent that, stage the
cycle by morphology — mitotic figures are unambiguous in the DNA channel —
and treat interphase as one class.

Two consequences for analysis of this dataset:

1. Prefer the ratio metrics (`*_cv*`, solidity, `*_radial_norm`, enrichment) for anything compared across files, and treat `*_mean_corr` / `*_integrated_corr` as within-session quantities.
2. `Flavopiridol` was imaged only in `SetC`, so condition and session are partly confounded. Compare it against the `SetC` controls, not the pooled ones. Doing that, its chromatin-contrast effect holds: CV 0.274 → 0.405 (p = 4e-09), PC1 −2.05 → +2.18 (p = 2e-09), mid-plane solidity 0.983 → 0.971 (p = 0.02).

## The second pipeline: antennas

`antenna3d` traces nuclear F-actin antennas and writes one graph per nucleus
(`graphs/<nucleus_uid>.graphml`) alongside a per-nucleus and a per-edge
table. It reads step 1's `labels/` folder rather than segmenting anything
itself, so its rows carry the same `nucleus_uid` and the two tables join.

```bash
python scripts/run_antennas.py     # antennas alone; point LABELS_DIR at step 1's output
python scripts/run_combined.py     # both pipelines in order, then the join
```

[`docs/antenna3d.md`](docs/antenna3d.md) is the reference for all of it: the
parameters, the measurements behind every default, the combined run and the
join, and the caveats. **Read the caveats before quoting a number from it** —
on the data this package was built from, the probe-absent control test fails,
and it declines some data `nucleus3d` accepts.

---

## Things that will bite you

**The channel name is not the stain.** Check with `describe_file` on every new dataset. One dataset in this project is named `...JF646Hoechst...` and contains no Hoechst channel — the DNA stain is on `Cy5`. Another has no DNA channel at all, and the corresponding MATLAB script sets `NucSegChannel = S5P_SegChannel`, i.e. it segments nuclei from the Pol II signal because there is nothing else. A wrong channel produces a full, confident, meaningless table.

**Volumes are slab volumes, not nuclear volumes.** These stacks are thin optical slabs through taller nuclei, truncated top and bottom. `volume_um3`is the volume *within the imaged slab*. This is why `min_volume_um3` defaults to 15 rather than the MATLAB `Nuc_min_vol = 40`, and why `max_area_um2` — the largest cross-section — is the more portable size measure. Do not compare volumes across datasets with different slab thicknesses.

**`ignore_z_border` matters more than it looks.** With nuclei cut off in z, a plain distance transform treats the cut faces as background, so the distance map becomes a z-dominated plateau and the watershed shreds single nuclei. The option replicates the end slices so the faces read as interior. Leave it on unless your stacks fully contain every nucleus.

**Otsu is bimodal; your field may not be.** A field with both bright mitotic figures and dim interphase nuclei is effectively trimodal and the threshold lands too high. On one validated control field this missed two dim nuclei out of ten; `thresh_factor=0.4` recovered both (enrichment 2.2× and 2.5× over background, solidity 0.82 and 0.97) without inflating the
rest. Panel 3 of the QC figure is what catches this.

**Mitotic cells segment correctly, but check anyway.** Condensed chromosome masses score solidity ≈ 0.84–0.85 and pass the default 0.7 filter — the worry that they would be rejected as non-convex turned out to be unfounded at this resolution. Verified on cultured-cell data containing obvious
mitotic figures.

**A "clean" missed-nucleus check can be wrong.** During development the first version of this check thresholded at the 99th percentile, found only bright objects, and reported zero missed nuclei while two dim ones sat unsegmented. The QC figure is what exposed it. `unsegmented_fraction` now uses Otsu and dilates the masks by 0.5 µm so boundary halo does not drown the signal.

**`nucleus_uid` maps forward only.** Build it from
`(file stem, position, label)`; never parse one back out of a path. A stem
can itself contain the separator — in this project's own data a field named
`..._sphere_postfix__crop` has a doubled underscore — and inverting the
encoding once silently dropped **24 of 80** control nuclei, a fifth of the
arm, with no error. Every table carries `file`, `position` and `label` as
their own columns so you never need to.

**Segmentation at native 100× sampling is impractical.** A 2280×2588×123
field is 726 M voxels, and the coarse DoG arm runs at σ = 217 px laterally:
measured, **~24 min per field** for the DoG alone, with the process reaching
18.8 GB and swapping on a 23 GB machine. `plan_workers` does not help — it
would allocate one worker for a field that size. Segmenting a 5× laterally
binned stack (0.23 µm, the grid `antenna3d` uses itself, and which its
`bin_factor` path already accepts) takes **~0.3 min** and 0.12 GB, and a
nucleus is still 35 px across. The toolbox does not yet offer that as a
setting; see the open issue.
