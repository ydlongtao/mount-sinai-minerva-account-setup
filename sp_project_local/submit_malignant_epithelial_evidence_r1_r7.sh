#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"
job=$(bsub < lsf/54_malignant_epithelial_evidence_r1_r7.lsf | awk -F'[<>]' '{print $2}')
printf 'Malignant epithelial evidence job: %s\n' "$job"
printf 'Monitor: bjobs %s; bpeek %s\n' "$job" "$job"
