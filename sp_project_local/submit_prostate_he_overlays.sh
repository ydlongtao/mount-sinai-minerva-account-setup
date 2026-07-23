#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"
job=$(bsub < lsf/50_make_prostate_annotation_he_overlays_array.lsf | awk -F'[<>]' '{print $2}')
printf 'H&E annotation overlay job: %s\n' "$job"
printf 'Monitor: bjobs %s; bpeek %s\n' "$job" "$job"
