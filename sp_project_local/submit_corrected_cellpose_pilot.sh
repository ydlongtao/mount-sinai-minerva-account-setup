#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
PREP=$(bsub < "$PROJECT_HOME/lsf/27_prepare_corrected_workspace.lsf" | sed -n 's/Job <\([0-9][0-9]*\)>.*/\1/p')
for sample in SC000895-R4 SC000895-R5 SC000895-R9; do
  bsub -w "done($PREP)" \
    -J "sp_corrected_cellpose_${sample}" \
    -P acc_DiseaseGeneCell -q gpuexpress -n 1 -gpu "num=1" \
    -R "span[hosts=1] rusage[mem=64000]" -W 04:00 -L /bin/bash \
    -o "$PROJECT_HOME/logs/%J.0.stdout" -e "$PROJECT_HOME/logs/%J.0.stderr" \
    "export SP_PILOT_SAMPLE=$sample; bash $PROJECT_HOME/lsf/28_corrected_cellpose_pilot_gpu.lsf"
done
printf 'Preparation job: %s\nPilot samples: SC000895-R4 SC000895-R5 SC000895-R9\n' "$PREP"
