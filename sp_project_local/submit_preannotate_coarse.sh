#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"
job=$(bsub < lsf/48_preannotate_coarse_clusters_array.lsf | awk -F'[<>]' '{print $2}')
printf 'Preannotation array submitted: %s\n' "$job"
printf 'Monitor: bjobs %s\n' "$job"
