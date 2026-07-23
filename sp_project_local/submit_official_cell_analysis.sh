#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"
job=$(bsub < lsf/45_official_cell_level_analysis_array.lsf | awk -F'[<>]' '{print $2}')
printf 'Official segmented cell-level analysis submitted: %s\n' "$job"
printf 'Monitor: bjobs %s; bpeek %s\n' "$job" "$job"
