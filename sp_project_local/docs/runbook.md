# Runbook

## Local Review

Before upload, check this local package:

```bash
cd /Users/huangfulongtao/Documents/Minerva服务器账户配置
find sp_project_local -maxdepth 3 -type f | sort
python3 -m py_compile sp_project_local/scripts/*.py
grep -R "/Users/huangfulongtao" sp_project_local/scripts sp_project_local/lsf sp_project_local/config
grep -R "#BSUB -P acc_DiseaseGeneCell" sp_project_local/lsf
```

Expected:

- Python scripts compile.
- No local Mac path appears in remote scripts, LSF scripts, or config.
- Every LSF script contains `#BSUB -P acc_DiseaseGeneCell`.
- The package contains no raw data or generated analysis outputs.

## Upload

Upload the contents of `sp_project_local/` to:

```bash
/sc/arion/work/huangl21/sp_project/
```

Do not upload raw data into the project directory. The pipeline reads raw data from:

```bash
/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome
```

## Submit Order

Run from Minerva:

```bash
cd /sc/arion/work/huangl21/sp_project
chmod +x submit_pipeline.sh
./submit_pipeline.sh preflight
```

Inspect `results/preflight/environment_report.md` and corresponding logs. Continue only if the status is `pass`.

```bash
./submit_pipeline.sh manifest
./submit_pipeline.sh pilot
```

Inspect the pilot outputs for `SC000895-R4`:

```bash
find results/pilot -maxdepth 3 -type f | sort
```

If the pilot is acceptable:

```bash
./submit_pipeline.sh samples
```

After all array jobs complete:

```bash
./submit_pipeline.sh downstream
```

## Monitoring

Use Terminal or the opened Jupyter session for monitoring only:

```bash
bjobs -u $USER
bjobs -l <JOBID>
bhist <JOBID>
bpeek <JOBID>
tail -n 100 logs/<JOBID>.<INDEX>.stderr
```

The notebook `notebooks/monitor_lsf_and_results.ipynb` previews logs, preflight reports, QC tables, figures, and final reports.

## Expected Outputs

- `results/preflight/environment_report.md`
- `config/sample_manifest.csv`
- `results/sample_qc_summary.csv`
- `work/<sample>/<sample>_008um.h5ad`
- `work/<sample>/<sample>_segmented.h5ad`
- `results/samples/<sample>/*.png`
- `results/prostate_marker_scores.csv`
- `results/integrated_visium_hd.h5ad`
- `results/integrated_umap.png`
- `results/summary_report.md`
- `results/summary_report.html`
