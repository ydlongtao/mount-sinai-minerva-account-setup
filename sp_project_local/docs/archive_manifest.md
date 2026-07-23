# Project Archive Manifest

Archive date: 2026-07-10

## Included

- scripts/: v1 and v2 Python analysis scripts, H5AD inspection, and local
  report generation.
- lsf/: preflight, setup, pilot, v2 QC, v2 GPU model, report, array, and
  downstream job wrappers.
- config/: path template, sample manifest template, pilot configuration, and
  prostate cancer marker panels.
- notebooks/: read-only LSF and result monitoring notebook.
- docs/: runbook, data structure notes, review checklist, and this archive
  boundary.
- README.md, AGENTS.md, and Agent.md: usage and maintenance instructions.

## Excluded

- Raw Visium HD matrices and barcode_mappings.parquet.
- H5AD, H5, BAM, BAI, cloupe, SVS, TIFF, and other large data files.
- local_results/, LSF logs, scratch files, checkpoints, and generated work
  directories.
- Passwords, MFA material, SSH private keys, tokens, and browser session data.

## Validated Pilot Outputs

Local results are kept outside GitHub at:

    /Users/huangfulongtao/Documents/Minerva服务器账户配置/sp_project_local/local_results/SC000895-R4_v2_20260710/

Important local files:

- SC000895-R4_v2_analysis_report.html
- model/SC000895-R4_016um_v2_full.h5ad
- model/SC000895-R4_008um_v2_mapped.h5ad
- model/clustering_resolution_summary.csv
- model/marker_availability.csv
- h5ad_structure.json

## Reproduction Boundary

The v2 pilot was validated on Minerva with QC job 256092651, GPU model job
256092679, and report job 256092680. These IDs are provenance only; future
runs must create new job IDs and must not reuse old outputs without checking
the input and parameter state.
