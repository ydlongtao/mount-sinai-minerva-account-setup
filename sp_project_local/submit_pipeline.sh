#!/bin/bash
set -euo pipefail

PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"

mkdir -p logs results work /sc/arion/scratch/huangl21/sp_project

usage() {
  cat <<'EOF'
Usage: ./submit_pipeline.sh <stage>

Stages:
  setup       Submit project directory setup only.
  preflight   Submit setup and compute-node environment preflight.
  manifest    Submit sample manifest and QC summary jobs.
  pilot       Submit SC000895-R4 pilot read/process job.
  samples     Submit 10-sample LSF array job.
  downstream  Submit marker annotation, integration, and report jobs.

Recommended order:
  ./submit_pipeline.sh preflight
  inspect results/preflight/environment_report.md and logs/
  ./submit_pipeline.sh manifest
  ./submit_pipeline.sh pilot
  inspect pilot outputs
  ./submit_pipeline.sh samples
  ./submit_pipeline.sh downstream
EOF
}

stage=${1:-}
case "$stage" in
  setup)
    bsub < lsf/01_setup_project.lsf
    ;;
  preflight)
    setup_job=$(bsub < lsf/01_setup_project.lsf | awk -F'[<>]' '{print $2}')
    bsub -w "done(${setup_job})" < lsf/00_preflight_env.lsf
    ;;
  manifest)
    manifest_job=$(bsub < lsf/02_build_sample_manifest.lsf | awk -F'[<>]' '{print $2}')
    bsub -w "done(${manifest_job})" < lsf/03_collect_qc_metrics.lsf
    ;;
  pilot)
    bsub < lsf/04_pilot_read_one_sample.lsf
    ;;
  samples)
    bsub < lsf/05_process_each_sample_array.lsf
    ;;
  downstream)
    marker_job=$(bsub < lsf/06_marker_annotation_prostate_cancer.lsf | awk -F'[<>]' '{print $2}')
    integrate_job=$(bsub < lsf/07_integrate_samples.lsf | awk -F'[<>]' '{print $2}')
    bsub -w "done(${marker_job}) && done(${integrate_job})" < lsf/08_make_summary_report.lsf
    ;;
  -h|--help|help|"")
    usage
    ;;
  *)
    echo "Unknown stage: $stage" >&2
    usage >&2
    exit 2
    ;;
esac
