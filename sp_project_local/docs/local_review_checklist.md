# Local Review Checklist

- [ ] Package is local only; nothing has been uploaded or submitted.
- [ ] `config/paths.env` points to `/sc/arion/work/huangl21/sp_project`.
- [ ] `config/paths.env` points scratch to `/sc/arion/scratch/huangl21/sp_project`.
- [ ] `config/paths.env` points raw data to `/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome`.
- [ ] Every LSF file contains `#BSUB -P acc_DiseaseGeneCell`.
- [ ] LSF stdout/stderr point to `/sc/arion/work/huangl21/sp_project/logs/%J.%I.stdout` and `.stderr`.
- [ ] Preflight is the first required compute-node job.
- [ ] Jupyter notebook is monitoring-only.
- [ ] `.gitignore` excludes results, work files, logs, and large data extensions.
