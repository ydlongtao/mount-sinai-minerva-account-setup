#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"
job=$(bsub < lsf/53_niche_marker_pathway_r1_r7.lsf | awk -F'[<>]' '{print $2}')
printf 'Niche marker/pathway job: %s\n' "$job"
printf 'Monitor: bjobs %s; bpeek %s\n' "$job" "$job"
