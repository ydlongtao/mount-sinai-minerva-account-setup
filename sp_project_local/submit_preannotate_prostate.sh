#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"
job=$(bsub < lsf/49_preannotate_prostate_clusters_array.lsf | awk -F'[<>]' '{print $2}')
printf 'Prostate-aware preannotation submitted: %s\n' "$job"
printf 'Monitor: bjobs %s\n' "$job"
