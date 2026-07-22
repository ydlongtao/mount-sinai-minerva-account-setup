#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"
job=$(bsub < lsf/47_repair_coarse_leiden_array.lsf | awk -F'[<>]' '{print $2}')
printf 'Repair coarse Leiden submitted: %s\n' "$job"
printf 'Monitor: bjobs %s\n' "$job"
