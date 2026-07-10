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

## Scope

Version 1 uses `square_008um` as the primary binned resolution and segmented cell output as a parallel read/process check. `square_016um` is recorded in the manifest for quick comparison. `square_002um`, BAM, molecule info, and cloupe files are intentionally not processed in this first pass.
