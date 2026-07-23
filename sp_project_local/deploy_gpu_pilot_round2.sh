#!/bin/bash
set -euo pipefail

LOCAL_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SSH_TARGET=${MINERVA_TARGET:-huangl21@minerva.hpc.mssm.edu}

FILES=(
  scripts/04_pilot_read_one_sample_gpu.py
  scripts/sp_utils.py
  lsf/04_pilot_read_one_sample_gpu.lsf
  lsf/04_pilot_read_one_sample.lsf
)

for file in "${FILES[@]}"; do
  if [[ ! -f "$LOCAL_ROOT/$file" ]]; then
    echo "Missing local file: $LOCAL_ROOT/$file" >&2
    exit 1
  fi
done

echo "Uploading round-2 pilot scripts to $SSH_TARGET"
echo "Complete the Minerva password and MFA prompts when requested."

tar -C "$LOCAL_ROOT" -czf - "${FILES[@]}" | ssh "$SSH_TARGET" '
set -euo pipefail

PROJECT=/sc/arion/work/huangl21/sp_project
SCRATCH=/sc/arion/scratch/huangl21/sp_project
SAMPLE=SC000895-R4
STAGE="$SCRATCH/deploy_gpu_pilot_round2_$$"

cleanup_stage() {
    rm -rf "$STAGE"
}
trap cleanup_stage EXIT

ACTIVE=$(bjobs -u "$USER" 2>/dev/null | awk "NR > 1 && \$3 ~ /sp_pilot/ {print}")
if [[ -n "$ACTIVE" ]]; then
    echo "A pilot job is still active; refusing cleanup:" >&2
    echo "$ACTIVE" >&2
    exit 1
fi

mkdir -p "$STAGE"
tar -xzf - -C "$STAGE"

STAMP=$(date +%Y%m%d_%H%M%S)
BACKUP="$PROJECT/work/deploy_backups/gpu_pilot_round2_$STAMP"
mkdir -p "$BACKUP" "$PROJECT/scripts" "$PROJECT/lsf" "$PROJECT/logs"

for file in \
    scripts/04_pilot_read_one_sample_gpu.py \
    scripts/sp_utils.py \
    lsf/04_pilot_read_one_sample_gpu.lsf \
    lsf/04_pilot_read_one_sample.lsf
do
    if [[ -f "$PROJECT/$file" ]]; then
        cp -p "$PROJECT/$file" "$BACKUP/$(basename "$file")"
    fi
    cp -p "$STAGE/$file" "$PROJECT/$file"
done

rm -rf \
    "$PROJECT/results/pilot/$SAMPLE" \
    "$PROJECT/results/pilot_gpu/$SAMPLE" \
    "$PROJECT/work/$SAMPLE" \
    "$PROJECT/work_gpu/$SAMPLE" \
    "$SCRATCH/256024573" \
    "$SCRATCH/256035794_gpu_pilot"

mkdir -p \
    "$PROJECT/results/pilot_gpu/$SAMPLE" \
    "$PROJECT/work_gpu/$SAMPLE" \
    "$SCRATCH"

chmod 755 \
    "$PROJECT/scripts/04_pilot_read_one_sample_gpu.py" \
    "$PROJECT/lsf/04_pilot_read_one_sample_gpu.lsf" \
    "$PROJECT/lsf/04_pilot_read_one_sample.lsf"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /sc/arion/work/huangl21/conda_envs/omicverse-gpu

python -m py_compile \
    "$PROJECT/scripts/04_pilot_read_one_sample_gpu.py" \
    "$PROJECT/scripts/sp_utils.py"

bash -n \
    "$PROJECT/lsf/04_pilot_read_one_sample_gpu.lsf" \
    "$PROJECT/lsf/04_pilot_read_one_sample.lsf"

python - <<"PY"
import rapids_singlecell as rsc

checks = [
    rsc.pp.normalize_total,
    rsc.pp.log1p,
    rsc.pp.highly_variable_genes,
    rsc.pp.scale,
    rsc.pp.pca,
    rsc.pp.neighbors,
    rsc.tl.umap,
    rsc.tl.leiden,
]
print("RAPIDS callable checks:", all(callable(fn) for fn in checks), len(checks))
PY

echo "Old scripts backed up at: $BACKUP"
echo "Old pilot output folders cleaned."

cd "$PROJECT"
SUBMIT_OUTPUT=$(bsub < lsf/04_pilot_read_one_sample_gpu.lsf)
echo "$SUBMIT_OUTPUT"

JOB_ID=$(printf "%s\n" "$SUBMIT_OUTPUT" | sed -n "s/.*Job <\([0-9][0-9]*\)>.*/\1/p")
if [[ -z "$JOB_ID" ]]; then
    echo "Unable to parse the submitted job ID." >&2
    exit 1
fi

echo "GPU pilot round 2 job: $JOB_ID"
bjobs "$JOB_ID"
echo "Monitor: bpeek $JOB_ID"
echo "GPU usage: tail -f $PROJECT/logs/$JOB_ID.gpu_dmon.log"
'
