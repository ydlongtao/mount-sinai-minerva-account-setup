#!/bin/bash
set -euo pipefail
PROJECT_HOME=${PROJECT_HOME:-/sc/arion/work/huangl21/sp_project}
export SP_PROJECT_HOME="$PROJECT_HOME"
export SP_PROJECT_SCRATCH=/sc/arion/scratch/huangl21/sp_project
export SP_PROJECT_MANIFEST="$PROJECT_HOME/config/sample_manifest_v2_corrected.csv"
export SP_OFFICIAL_RESULTS="$PROJECT_HOME/results/segmented_official_v1"
export SP_OFFICIAL_WORK="$PROJECT_HOME/work/segmented_official_v1"
export PYTHONPATH="$PROJECT_HOME/scripts:${PYTHONPATH:-}"
