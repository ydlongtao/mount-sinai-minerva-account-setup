#!/usr/bin/env bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
job=$(bsub < "$PROJECT_HOME/lsf/58_infercnv_sensitivity_r1_r7_array.lsf" | awk -F'[<>]' '{print $2}')
printf 'infercnvpy sensitivity array submitted: %s\n' "$job"
printf 'Monitor with: bjobs %s; bpeek %s[1]\n' "$job" "$job"
