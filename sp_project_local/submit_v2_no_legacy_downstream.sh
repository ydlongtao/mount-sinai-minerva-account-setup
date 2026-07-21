#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"

extract_job_id() { sed -n 's/Job <\([0-9][0-9]*\)>.*/\1/p'; }

J016=$(bsub < lsf/22_integrate_016um_v2_no_legacy.lsf | extract_job_id)
JCELL=$(bsub < lsf/23_integrate_cell_v2_no_legacy.lsf | extract_job_id)
JDOM=$(bsub -w "done(${J016}) && done(${JCELL})" < lsf/24_spatial_domains_v2_no_legacy.lsf | extract_job_id)
JREPORT=$(bsub -w "done(${JDOM})" < lsf/25_batch_report_v2_no_legacy.lsf | extract_job_id)

printf 'v2 downstream jobs:\n  integrate_016um: %s\n  integrate_cell: %s\n  spatial_domains: %s\n  report: %s\n' "$J016" "$JCELL" "$JDOM" "$JREPORT"
