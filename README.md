# Nuclei-And-Inside
Repository to support co-development of microscopy image analysis of nuclei as well as objects and structures within them.

3D nucleus segmentation and per-nucleus quantification for Nikon `.nd2` z-stacks.

Python port of the lab MATLAB pipeline (`PseudoTimeCourse_CroppedImages.m` + `otsuLimit.m`), with a watershed split for touching nuclei and a boundary-refinement step that the MATLAB version does not have.
Original pipeline: [https://github.com/lhilbert/VisitorGene_PseudoTime](https://github.com/lhilbert/VisitorGene_PseudoTime)

---

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

```bash
conda create -n nucleus3d python=3.11
conda activate nucleus3d
pip install -e .
```

Dependencies: `nd2`, `numpy`, `scipy`, `scikit-image`, `pandas`, `tifffile`, `matplotlib`.

---

## Example data and code validation

Validated on two datasets: vt-iSIM fixed zebrafish embryos (JF646-Hoechst) and cultured cells (DAPI, drug conditions).

You can download one of these data sets from the following, publicly shared Zenodo repository. The analysis should run fine on these data.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5242952.svg)](https://doi.org/10.5281/zenodo.5242952)



---

## Use

**1. Check your channels first.** This takes two seconds and prevents the
most expensive mistake in the pipeline:

```bash
python -c "from nucleus3d import describe_file; print(describe_file('yourfile.nd2'))"
```

**2. Copy `run_analysis.py`**, edit the `SETTINGS` block (input folder, output folder, DNA channel, whether to export nucleus boxes), and run it:

```bash
python run_analysis.py
```

**3. Look at the QC figures** in `<output>/qc/` before you use the table.

---

## Layout

```
nucleus3d/
    io.py         read .nd2 -> Stack (pixels + voxel size + channel names)
    segment.py    segmentation: DoG -> Otsu -> watershed -> refine -> filter
    quantify.py   per-nucleus intensities, all channels, background-corrected
    export.py     one 3D OME-TIFF per nucleus (all channels + mask)
    validate.py   per-field QC figure and QC numbers
    pipeline.py   batch driver: walks files/fields, writes the output tree
run_analysis.py   TEMPLATE -- the only file you normally edit
```

The split is deliberate: `run_analysis.py` holds parameters and no logic, the modules hold logic and no parameters. To script something custom, import the package and skip the template:

```python
from nucleus3d import load_field, segment_nuclei, quantify_nuclei, SegParams

stack = load_field("field.nd2", position=3)
labels, props = segment_nuclei(stack.channel("DAPI"), stack.voxel_um,
                               SegParams(thresh_factor=0.5))
table = quantify_nuclei(stack, labels, "DAPI")
```

---

## Output

```
<output>/
    nuclei_measurements.csv   one row per nucleus
    field_summary.csv         one row per field: QC numbers
    run_parameters.json       every parameter used
    qc/                       one validation figure per field
    nucleus_boxes/            one OME-TIFF per nucleus    (if enabled)
    nucleus_boxes_index.csv   nucleus -> file map          (if enabled)
    failures.csv              fields that errored          (if any)
```

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

---

## Reading the QC figure

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

## Tuning

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

## Things that will bite you

**The channel name is not the stain.** Check with `describe_file` on every new dataset. One dataset in this project is named `...JF646Hoechst...` and contains no Hoechst channel — the DNA stain is on `Cy5`. Another has no DNA channel at all, and the corresponding MATLAB script sets `NucSegChannel = S5P_SegChannel`, i.e. it segments nuclei from the Pol II signal because there is nothing else. A wrong channel produces a full, confident, meaningless table.

**Volumes are slab volumes, not nuclear volumes.** These stacks are thin optical slabs through taller nuclei, truncated top and bottom. `volume_um3`is the volume *within the imaged slab*. This is why `min_volume_um3` defaults to 15 rather than the MATLAB `Nuc_min_vol = 40`, and why `max_area_um2` — the largest cross-section — is the more portable size measure. Do not compare volumes across datasets with different slab thicknesses.

**`ignore_z_border` matters more than it looks.** With nuclei cut off in z, a plain distance transform treats the cut faces as background, so the distance map becomes a z-dominated plateau and the watershed shreds single nuclei. The option replicates the end slices so the faces read as interior. Leave it on unless your stacks fully contain every nucleus.

**Otsu is bimodal; your field may not be.** A field with both bright mitotic figures and dim interphase nuclei is effectively trimodal and the threshold lands too high. On one validated control field this missed two dim nuclei out of ten; `thresh_factor=0.4` recovered both (enrichment 2.2× and 2.5× over background, solidity 0.82 and 0.97) without inflating the
rest. Panel 3 of the QC figure is what catches this.

**Mitotic cells segment correctly, but check anyway.** Condensed chromosome masses score solidity ≈ 0.84–0.85 and pass the default 0.7 filter — the worry that they would be rejected as non-convex turned out to be unfounded at this resolution. Verified on cultured-cell data containing obvious
mitotic figures.

**A "clean" missed-nucleus check can be wrong.** During development the first version of this check thresholded at the 99th percentile, found only bright objects, and reported zero missed nuclei while two dim ones sat unsegmented. The QC figure is what exposed it. `unsegmented_fraction` now uses Otsu and dilates the masks by 0.5 µm so boundary halo does not drown the signal.
