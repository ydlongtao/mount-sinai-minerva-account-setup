# AGENTS.md

## Project Scope

This directory is a reproducible code archive for the Minerva 10x Visium HD
human prostate cancer pilot. It contains analysis code and execution
documentation, not raw data or production outputs.

## Required Invariants

- Raw data are read-only under
  /sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome.
- Work files use /sc/arion/work/huangl21/sp_project/.
- Temporary files use /sc/arion/scratch/huangl21/sp_project/.
- Every LSF script uses project account acc_DiseaseGeneCell.
- Do not silently overwrite v1 or v2 outputs.
- Do not commit H5AD, BAM, SVS, raw matrices, logs, or local results.

## V2 Decisions

- Use 16 um for global modeling.
- Map labels to 8 um for high-resolution review.
- Preserve full genes and the counts layer in the final 16 um H5AD.
- Check graph connectivity before UMAP and Leiden.
- Record GPU/CPU fallback behavior in run_state.json.

## Validation

    python3 -m py_compile scripts/*.py
    for f in lsf/*.lsf; do bash -n "$f"; done
    rg -n '#BSUB -P acc_DiseaseGeneCell' lsf

Use docs/runbook.md for Minerva submission and monitoring. Use
docs/archive_manifest.md before publishing changes.
