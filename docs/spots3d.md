# spots3d — Oligopaint spots inside nuclei

A port of the Oligopaint (OP) block of the Hilbert-lab MATLAB script
[`Segmentation_FeatureExtraction.m`](https://github.com/lhilbert/VisitorGene_PseudoTime/blob/10fc17207e20162d354ba53ea3fcaeabdf45ca1e/Segmentation_FeatureExtraction.m).
Line numbers below refer to that file at commit `10fc172` (the file was last
changed in `669688e`).

`spots3d` is a sibling of `antenna3d`. It does **not** segment nuclei. It reads
nucleus3d's `labels/<stem>_pNN.tif`, and its rows carry the same `nucleus_uid`.

```bash
python scripts/run_segmentation.py    # nuclei -> results/labels/
python scripts/run_spots.py           # spots  -> results_spots/
```

## What is ported, in MATLAB order

| Step | MATLAB | spots3d |
|---|---|---|
| Nuclear mask eroded in-plane by `round(0.5 µm / pixel)` (11 px at 0.046 µm) | l.455–459 | `detect.erode_mask_union` |
| Crop every channel to the bounding box of the eroded nucleus | l.656–701 | `pipeline.detect_field` |
| Background correction: 2D DoG, `imgaussfilt(σs) − imgaussfilt(σl)`, per z-plane, on the crop; `I − G(σl)` when σs < 0.3 px | l.725–732 | `filters.dog_matlab` |
| Threshold `DoG .* mask > mean + 5·SD`, mean/SD in float64 over the nucleus mask (SD with N−1) | l.755–759 | `detect.threshold_in_nucleus` |
| 18-connected components | l.795 | `detect.components_gated` |
| Volume gate `n_voxels > 0.05 µm³ / voxel volume` (strict; 118.15 → ≥ 119 voxels here) | l.796–800 | `detect.components_gated` |
| DBSCAN, eps 0.2 µm, minPts 1, on unweighted component centroids in µm; voxel union per cluster. Runs **after** the volume gate | l.1223–1243 | `detect.dbscan_merge` |
| Volume, solidity, centroid | l.1278–1286 | `features.spot_features` |
| Elongation: major/minor axis (MATLAB ellipse moments, +1/12) on plane `ceil(nz/2)` of the spot's bbox; NaN if that plane is empty | l.1307–1359 | `features.elongation_matlab` |
| Intensity: median(spot) / median(nucleus) in every quantified channel, after that channel's blur (0.07 µm for OP) | l.1380–1406 | `features.spot_intensities` |
| Nucleus mean intensity per channel (raw, over the eroded mask) | l.632–638 | `nuc_mean_<channel>` |
| ">4 spots in a nucleus" | l.1993–2050 | flag only: `qc_more_than_max_spots`, `qc_any_nucleus_more_than_max` |

**There is no shape filter.** MATLAB measures solidity and elongation but gates OP spots on volume only, and so does this port.

Component and spot numbering follows MATLAB's `bwconncomp` (column-major `[Y X Z]` order) and `dbscan`/`unique` (first member). This is meant to make spot *k* here the same object as spot *k* in a MATLAB run on the same mask. That has not been checked against MATLAB yet (see "Validating against MATLAB").

## Not ported (yet)

- Pol II S5P/S2P cluster detection.
- OP–cluster pairing.
- The random-position reference (l.1436–1786).
- Cytoplasmic shell intensity, condition pooling, statistics and figures (l.2298 onwards).
- MATLAB's own nuclear segmentation: nucleus3d's labels are used instead.

## Defaults (all in `spots3d/params.py`, all configurable)

| Parameter | Default | MATLAB |
|---|---|---|
| `sigma_small_um` | 0.07 | `OP_segBlurSigma_small`, l.100 |
| `sigma_large_um` | 0.3 | `OP_segBlurSigma_large`, l.101 |
| `n_sd` | 5 | `OP_seg_numStdDev`, l.102 |
| `dbscan_eps_um` | 0.2 | `OP_DBSCAN_epsilon`, l.104 |
| `min_vol_um3` | 0.05 | `OP_minVol`, l.110 |
| `quant_blur_um` | 0.07 | `quantBlurSigma`, l.66 |
| `erosion_um` | 0.5 | `nuc_segErosion`, l.70 |
| `qc_max_spots` | 4 | `fractionAboveFour`, l.1993 |
| `image_dtype` | float32 | `imageClass = 'single'`, l.46 |

The TripleOP channel layout (verified on all 15 files) is `[0 Blank, 1 Cy5, 2 Blank, 3 Mcherry, 4 GFP, 5 Blank, 6 Dapi, 7 Blank]`:

| Probe | Dye | Channel index | ND2 name |
|---|---|---|---|
| OP1 | Cy3 | 3 | Mcherry |
| OP2 | A488 | 4 | GFP |
| OP3 | Cy5 | 1 | Cy5 |

Channels are read **by index**, and the name is only checked.

## Grid

Detection runs on the native .nd2 grid (0.0459985 × 0.0459985 × 0.2 µm). A label volume is accepted only if it has OME PhysicalSize equal to the image's and the same shape. nucleus3d's `labels/` meet this condition.

Binned or uncalibrated masks, such as a Cellpose run on 2×-downsampled DAPI (0.091997 µm), are refused rather than rescaled.

## Known deviations from MATLAB

- **Erosion structuring element:** an exact Euclidean disk. MATLAB's `strel('disk', r)` is a polygonal (N=4) approximation, so the eroded rim differs by up to a pixel or two along diagonals.
- **Erosion is of the union of the *kept* nucleus3d labels.** MATLAB erodes its own binary mask before filtering nuclei.
  - Faces between touching nuclei are not eroded in either.
  - A nucleus that erosion splits stays one `nucleus_uid`; MATLAB would re-label the pieces.
- **Nuclear masks themselves:** watershed and boundary refinement in nucleus3d against DoG/Otsu in MATLAB. This is the largest source of difference on real data, and it is intentional.
- **Solidity:** skimage's voxel-count convex hull against MATLAB's `regionprops3` ConvexVolume. They differ by a few % for small objects. Solidity is reported only, never used for gating.
- **Comparison precision:** the comparison of the float32 DoG with a float64 threshold is done in float64. MATLAB mixes single/double. This affects only exact ties.
- **DBSCAN:** uses `distance <= eps` (rangesearch semantics). Not verifiable without MATLAB; it only matters at exact ties.
- **QC flag:** MATLAB counts only spots in nuclei that also contain a large Pol II cluster ("paired"). Here every detected spot counts.
- **Empty nuclei:** a nucleus with 0 spots gets `median_spot_volume_um3 = NaN`. MATLAB leaves its preallocated 0.

## Validating against MATLAB (not done; needs a licence)

**Numerical MATLAB–Python parity has not been validated.** The port follows the MATLAB code line by line, and the unit tests check each step against what the MATLAB code does, but no output has yet been compared with a MATLAB run on the same data.

1. MATLAB R2019b+ with the Image Processing and Statistics & Machine Learning toolboxes.
2. A `.mat` per field in the extraction format (`imgStack`, `imgSize`, `pixelSize`, `zStepSize`, `condInd`, `condName`). It can be written from Python with `scipy.io.savemat`.
3. A harness that:
   - loads that `.mat` and **this repo's label TIFF** (`tiffreadVolume`);
   - builds `comps` from the labels;
   - runs l.686–759, 795–800 and 1220–1406 verbatim;
   - saves the per-spot results.

   Using the same masks isolates the spot step from the nuclear-segmentation difference.
4. Expected agreement:
   - counts and voxel counts per spot: identical;
   - centroids: compare `matlab_cent_{x,y,z}_um` with `OP_centCell`, ≤ 1e-6 µm;
   - intensity ratios: relative ≤ 1e-5;
   - elongation: ≤ 1e-6;
   - solidity: documented tolerance only.

## A bug in nucleus3d this exposed (not fixed here)

`nucleus3d.core.quantify.quantify_nuclei` keys channels by **name**. The TripleOP files have four channels called `Blank` (indices 0, 2, 5, 7), with two effects:

- `stack.channel("Blank")` always returns index 0, so `Blank_bg` is index 0's background.
- The per-channel loop overwrites the `Blank_*` columns, so they hold index 7's values.

`results/nuclei_measurements.csv` shows the mismatch: `Blank_bg = 274` against `Blank_median = 50`, giving `Blank_enrichment = 0.18`. Columns of uniquely named channels (Cy5, Mcherry, GFP, Dapi) are unaffected. spots3d never looks channels up by name.
