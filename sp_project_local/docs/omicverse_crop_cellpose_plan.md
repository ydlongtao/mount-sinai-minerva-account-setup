# OmicVerse Crop, Rotation, and Cellpose Pilot

This is the next pilot stage after the validated Visium HD v2 workflow.

## Method Basis

- Crop and rotate follows:
  https://omicverse.readthedocs.io/zh-cn/latest/Tutorials-space/t_crop_rotate.html
- Cellpose and bin-to-cell follows:
  https://omicverse.readthedocs.io/zh-cn/latest/Tutorials-space/t_cellpose.html
- The workflow uses 2 um bins, H&E as the primary segmentation source, GEX
  segmentation as a secondary rescue source, and bin2cell aggregation.

## Staged Scripts

1. 11_pilot_crop_rotate_visium_hd.py tests coordinate crop and rotation on the
   validated 8 um mapped object without altering v2.
2. 12_pilot_cellpose_he.py loads the R4 2 um matrix and tissue image, crops a
   1000 um pilot region, and runs H&E Cellpose plus label expansion.
3. 13_pilot_cellpose_gex_bin2cell.py adds GEX labels, reconciles labels,
   aggregates bins to cells, and optionally exports SpaceRanger v4 format.
4. 14_cellpose_preflight.py checks the Cellpose dependency, image path, and
   2 um input path before any segmentation job.

## Current Pilot Defaults

- Sample: SC000895-R4
- Crop: x 2000-3000 um, y 2000-3000 um
- Image scale: 0.3239732935 um/pixel
- H&E buffer: 150 bins
- H&E flow threshold: 0.4
- GEX sigma: 5
- GPU host: lg02e02 H100

The crop coordinates are placeholders for the first visual test and must be
reviewed against the H&E overlay before using another region.

## Acceptance Checks

- Crop and rotation return valid spatial coordinates and expected bin counts.
- H&E output contains labels_he and labels_he_expanded.
- GEX output contains labels_gex and labels_joint.
- Cell-level output contains geometry or spatial cell coordinates.
- Exported SpaceRanger directory can be read back with
  ov.io.read_visium_hd(data_type='cellseg').
- H&E overlay and cell boundaries are visually registered before biological
  interpretation.

Do not submit these LSF jobs until the crop window and Cellpose dependencies
have been reviewed.

The current omicverse-gpu environment does not yet contain the cellpose Python
package. The preflight job is expected to fail until that dependency is
installed and verified.
