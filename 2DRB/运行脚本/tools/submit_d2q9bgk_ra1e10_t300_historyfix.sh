#!/usr/bin/env bash
set -euo pipefail

new_root=/data2/XLLi/LBMCDE/THRERMAL-D2Q9BGK-TEST/FLOW_BGK_CHI_EQUAL_SCAN/Ra1e10_restart100_latest_t300_grid2049_historyfix
cases=(
  chinu0.0_chikappa0.0
  chinu0.2_chikappa0.2
  chinu0.5_chikappa0.5
  chinu0.7_chikappa0.7
  chinu0.9_chikappa0.9
  chinu0.99_chikappa0.99
)

previous_job=
for case_name in "${cases[@]}"; do
  pbs="${new_root}/${case_name}/run_resume.pbs"
  if [[ -z "${previous_job}" ]]; then
    previous_job=$(qsub "${pbs}")
  else
    previous_job=$(qsub -W "depend=afterany:${previous_job}" "${pbs}")
  fi
  printf '%s %s\n' "${case_name}" "${previous_job}"
done
