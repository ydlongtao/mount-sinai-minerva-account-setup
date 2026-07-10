# Agent.md

Operational notes for future agents working on `sp_project_local`.

## Invariants

- Keep raw data read-only under
  `/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome`.
- Keep production work under `/sc/arion/work/huangl21/sp_project`.
- Keep temporary files under `/sc/arion/scratch/huangl21/sp_project`.
- Every LSF script must include `#BSUB -P acc_DiseaseGeneCell`.
- Do not commit `local_results/`, `results/`, `logs/`, `work/`, H5AD/H5, BAM,
  SVS, TIFF, or credentials.

## Validated Pilot

The current successful pilot is `SC000895-R4`.

Key outputs produced on Minerva:

- `results/pilot_v2/SC000895-R4/cellpose_gex/SC000895-R4_cell_level_cellpose.h5ad`
- `results/pilot_v2/SC000895-R4/cellpose_gex/SC000895-R4_cell_level_cellpose_raw_counts.h5ad`
- `results/pilot_v2/SC000895-R4/cellpose_gex/cellpose_spaceranger_output/graphclust_annotated_cell_segmentations.geojson`
- `results/pilot_v2/SC000895-R4/cellpose_gex/cellpose_spaceranger_raw_counts_output/filtered_feature_cell_matrix.h5`

The OmicVerse Space Ranger export originally skipped complex geometries because
its internal WKT parser only handled simple `POLYGON` strings. The local script
`scripts/13_pilot_cellpose_gex_bin2cell.py` overwrites that GeoJSON with a
Shapely-based exporter that supports `Polygon` and `MultiPolygon`; keep this
patch unless upstream OmicVerse fixes the parser.

## Current Follow-up

- `scripts/15_export_cellpose_raw_counts.py` aggregates raw integer 2um counts
  to the Cellpose/bin2cell cells and should be used for formal count matrix
  export.
- `scripts/16_gex_cellpose_parameter_sweep.py` tests whether GEX Cellpose
  parameters can produce non-zero `labels_gex`; the earlier baseline produced
  zero GEX labels, so H&E-expanded labels drove the successful `labels_joint`.
- `scripts/17_cellpose_qc_marker_report.py` creates the QC, marker score, and
  spatial overlay HTML report for review.
