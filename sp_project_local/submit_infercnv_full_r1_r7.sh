#!/usr/bin/env bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
job=$(bsub < "$PROJECT_HOME/lsf/57_infercnv_full_r1_r7_array.lsf" | awk -F'[<>]' '{print $2}')
printf 'infercnvpy full array submitted: %s\n' "$job"
printf 'Samples: R1 and R7; monitor with: bjobs %s; bpeek %s\n' "$job" "$job"
