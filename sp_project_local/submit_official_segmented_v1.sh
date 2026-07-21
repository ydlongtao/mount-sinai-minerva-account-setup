#!/bin/bash
set -euo pipefail
PROJECT_HOME=/sc/arion/work/huangl21/sp_project
cd "$PROJECT_HOME"
validate=$(bsub < lsf/42_validate_official_segmented_inputs.lsf | awk -F'[<>]' '{print $2}')
cells=$(bsub -w "done($validate)" < lsf/43_process_official_segmented_array.lsf | awk -F'[<>]' '{print $2}')
report=$(bsub -w "done($cells)" < lsf/44_make_official_segmented_report.lsf | awk -F'[<>]' '{print $2}')
printf 'Official segmented v1 submitted:\n  validation: %s\n  sample array: %s\n  report: %s\n' "$validate" "$cells" "$report"
