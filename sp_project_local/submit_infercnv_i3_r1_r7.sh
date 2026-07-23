#!/usr/bin/env bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
job=$(bsub < "$PROJECT_HOME/lsf/59_infercnv_i3_r1_r7_array.lsf" | awk -F'[<>]' '{print $2}')
printf 'R inferCNV i3 array submitted: %s\n' "$job"
printf 'Monitor with: bjobs %s; bpeek %s[1]\n' "$job" "$job"
