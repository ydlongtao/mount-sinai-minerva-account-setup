#!/bin/bash
set -euo pipefail
PROJECT_HOME=${PROJECT_HOME:-/sc/arion/work/huangl21/sp_project}
export SP_PROJECT_HOME="$PROJECT_HOME"
export SP_PROJECT_SCRATCH=/sc/arion/scratch/huangl21/sp_project
export SP_PROJECT_MANIFEST="$PROJECT_HOME/config/sample_manifest_v2_no_legacy.csv"
export SP_PROJECT_BATCH_RESULTS="$PROJECT_HOME/results/batch_v2_no_legacy"
export SP_PROJECT_BATCH_WORK="$PROJECT_HOME/work/batch_v2_no_legacy"
export SP_PROJECT_CONFIG="$PROJECT_HOME/config/batch_analysis_v2_no_legacy.yaml"
export PYTHONPATH="$PROJECT_HOME/scripts:${PYTHONPATH:-}"
