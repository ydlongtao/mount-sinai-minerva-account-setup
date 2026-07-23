#!/bin/bash
set -euo pipefail

PROJECT_HOME=${PROJECT_HOME:-/sc/arion/work/huangl21/sp_project}
cd "$PROJECT_HOME"
mkdir -p logs results/batch work/batch

submit() {
  local label=$1
  shift
  local output
  output=$(bsub "$@")
  local jobid
  jobid=$(printf '%s\n' "$output" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p')
  if [[ -z "$jobid" ]]; then
    echo "Failed to parse job id for $label" >&2
    echo "$output" >&2
    exit 1
  fi
  echo "$label $jobid"
}

qc_job=$(submit "batch_qc" < lsf/05_batch_qc.lsf | awk '{print $2}')
model16_job=$(submit "model_016um" -w "done($qc_job)" < lsf/06_process_016um_array.lsf | awk '{print $2}')
map8_job=$(submit "map_008um" -w "done($model16_job)" < lsf/07_process_008um_array.lsf | awk '{print $2}')
cellpose_job=$(submit "cellpose_he_gpu" -w "done($qc_job)" < lsf/08_cellpose_he_array_gpu.lsf | awk '{print $2}')
bin2cell_job=$(submit "bin2cell_raw_counts" -w "done($cellpose_job)" < lsf/09_bin2cell_raw_counts_array.lsf | awk '{print $2}')
cell_qc_job=$(submit "cell_qc_marker" -w "done($bin2cell_job)" < lsf/10_cell_qc_marker_array.lsf | awk '{print $2}')
integrate16_job=$(submit "integrate_016um" -w "done($model16_job)" < lsf/11_integrate_016um.lsf | awk '{print $2}')
integrate_cell_job=$(submit "integrate_cell_level" -w "done($cell_qc_job)" < lsf/12_integrate_cell_level.lsf | awk '{print $2}')
domains_job=$(submit "spatial_domains" -w "done($cell_qc_job)" < lsf/13_spatial_domains_array.lsf | awk '{print $2}')
report_job=$(submit "batch_report" -w "done($map8_job) && done($integrate16_job) && done($integrate_cell_job) && done($domains_job)" < lsf/14_batch_summary_report.lsf | awk '{print $2}')

cat <<EOF

Submitted Visium HD batch v1:
  batch_qc:              $qc_job
  model_016um array:     $model16_job
  map_008um array:       $map8_job
  cellpose_he GPU array: $cellpose_job
  bin2cell raw array:    $bin2cell_job
  cell QC array:         $cell_qc_job
  integrate 016um:       $integrate16_job
  integrate cell level:  $integrate_cell_job
  spatial domains array: $domains_job
  final report:          $report_job

Monitor:
  bjobs $qc_job $model16_job $map8_job $cellpose_job $bin2cell_job $cell_qc_job $integrate16_job $integrate_cell_job $domains_job $report_job
  bpeek <JOBID>

Final report:
  $PROJECT_HOME/results/batch/reports/spatial_batch_summary.html
EOF
