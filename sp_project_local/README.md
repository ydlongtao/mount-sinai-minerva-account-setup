# Minerva Visium HD Prostate Cancer Pipeline

This package is a local, reviewable script bundle for the 10x Visium HD human prostate cancer dataset on Minerva.

It is designed to be uploaded to:

```bash
/sc/arion/work/huangl21/sp_project/
```

Production jobs are submitted through LSF with project account `acc_DiseaseGeneCell`. Raw data are treated as read-only under:

```bash
/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome
```

Scratch and temporary files are placed under:

```bash
/sc/arion/scratch/huangl21/sp_project/
```

## Contents

- `scripts/`: Python scripts for setup, preflight, manifest, QC, pilot analysis, per-sample processing, marker scoring, integration, and reporting.
- `lsf/`: LSF job scripts with fixed Minerva paths and `#BSUB -P acc_DiseaseGeneCell`.
- `config/`: path configuration, prostate cancer marker genes, and a sample manifest template.
- `notebooks/`: a read-only monitoring notebook for LSF status, logs, tables, figures, and reports.
- `docs/`: runbook and data structure notes.
- `submit_pipeline.sh`: remote submission helper. Run only after uploading to Minerva.

## First Remote Commands After Upload

```bash
cd /sc/arion/work/huangl21/sp_project
chmod +x submit_pipeline.sh
./submit_pipeline.sh preflight
```

Then inspect:

```bash
results/preflight/environment_report.md
logs/
```

Continue only after preflight passes.

## Current Pilot Status

The validated pilot sample is `SC000895-R4` from the 10x Visium HD human
prostate cancer dataset.

The current branch includes:

- v2 pilot QC and modeling at `square_016um` with mapped high-resolution review
  at `square_008um`.
- Crop/rotate testing for the Visium HD analysis window.
- H&E Cellpose segmentation on the cropped `square_002um` object.
- bin2cell aggregation to cell-level objects.
- A patched Space Ranger-style GeoJSON export that supports both `Polygon` and
  `MultiPolygon` geometries via Shapely. This avoids the OmicVerse simple WKT
  parser limitation that only exported simple polygons.
- A raw-count export step that re-aggregates counts from the original
  `square_002um/filtered_feature_bc_matrix.h5` instead of exporting normalized
  float values.
- Official Space Ranger segmented-cell validation and an eight-sample review
  decision record in `docs/official_segmented_review_decisions.md`.
- A CPU array for the next per-sample cell-level pass:
  `submit_official_cell_analysis.sh`.

Local result downloads and generated reports are intentionally ignored by git
under `sp_project_local/local_results/`.

## Cellpose / bin2cell Pilot Order

Run these after the v2 pilot inputs are present on Minerva:

```bash
cd /sc/arion/work/huangl21/sp_project
bsub < lsf/04g_cellpose_preflight.lsf
bsub < lsf/04d_pilot_crop_rotate.lsf
bsub < lsf/04e_pilot_cellpose_he.lsf
bsub < lsf/04f_pilot_cellpose_gex_bin2cell.lsf
bsub < lsf/04h_export_cellpose_raw_counts.lsf
bsub < lsf/04i_gex_cellpose_parameter_sweep.lsf
bsub -w "done(<SWEEP_JOB_ID>)" < lsf/04j_cellpose_qc_marker_report.lsf
```

The successful cellpose/bin2cell pilot produced `879` cell-level objects and a
complete GeoJSON with `879` features and `0` skipped geometries.

## Scope

BAM, molecule info, and cloupe files are not processed in this pilot branch.
The `square_002um` data are used only for the crop-level Cellpose/bin2cell
pilot; broader all-sample production should be submitted only after reviewing
the pilot QC/marker report.

## Official Segmented-Cell Next Stage

The official segmented-cell review confirmed that R4, R2, and R6 have lower
complexity consistent with tissue and sequencing context. They are retained,
and no special low-RNA immune-cell threshold is applied. The next stage keeps
all official cells, normalizes each sample, selects 2,000 HVGs, computes PCA,
neighbors, Scanpy/igraph Leiden clusters, UMAP, prostate marker scores, and
spatial cluster previews. Zero-count objects remain in the saved H5AD but are
excluded from the graph calculation only.

After uploading this bundle to Minerva:

```bash
cd /sc/arion/work/huangl21/sp_project
chmod +x submit_official_cell_analysis.sh
./submit_official_cell_analysis.sh
```

Review the per-sample cluster and marker outputs before starting cross-sample
integration, spatial domains, and sample-level differential analysis.
