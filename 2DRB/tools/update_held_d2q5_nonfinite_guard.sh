#!/usr/bin/env bash
# Refresh only dependency-held D2Q5 Luo-TRT cases from the current guarded template.
set -euo pipefail

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
TEMPLATE=${ROOT}/analysis/current_d2q5_luo_trt_stats750_1000.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90

refresh_held_case() {
    local job_id=$1
    local case_dir=$2
    local grid=$3
    local rayleigh=$4
    local chi_nu=$5
    local state source_file temporary_file

    state=$(qstat -f "${job_id}" | awk '/job_state =/{print $3; exit}')
    if [[ "${state}" != H ]]; then
        printf 'SKIP %s: scheduler state is %s, expected H\n' "${job_id}" "${state:-unknown}"
        return 0
    fi

    source_file=${case_dir}/source/${SOURCE_NAME}
    temporary_file=${source_file}.guarded.$$
    cp "${TEMPLATE}" "${temporary_file}"
    sed -i \
        -e "s@integer(kind=4), parameter :: nx=2048, ny=2048@integer(kind=4), parameter :: nx=${grid}, ny=${grid}@" \
        -e "s@real(kind=8), parameter :: Rayleigh=1.0d10@real(kind=8), parameter :: Rayleigh=${rayleigh}@" \
        -e "s@real(kind=8), parameter :: chi_nu=0.5d0@real(kind=8), parameter :: chi_nu=${chi_nu}@" \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^!#define FLOW_ODD_EFFECTIVE_MAGIC/#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        "${temporary_file}"

    grep -q "nx=${grid}, ny=${grid}" "${temporary_file}"
    grep -q "Rayleigh=${rayleigh}" "${temporary_file}"
    grep -q "chi_nu=${chi_nu}" "${temporary_file}"
    grep -q '^#define FLOW_ODD_EFFECTIVE_MAGIC' "${temporary_file}"
    grep -q 'subroutine check_nonfinite_state' "${temporary_file}"
    grep -q 'error stop 86' "${temporary_file}"

    mv "${temporary_file}" "${source_file}"
    sha256sum "${TEMPLATE}" > "${case_dir}/source/source_template.sha256"
    sha256sum "${source_file}" > "${case_dir}/source/source.sha256"
    printf 'nonfinite_guard=macro_fields_every_1_tff_and_statistics_before_write\n' >> \
        "${case_dir}/case_manifest.txt"
    printf 'UPDATED %s %s\n' "${job_id}" "${case_dir}"
}

[[ -f "${TEMPLATE}" ]]
grep -q 'subroutine check_nonfinite_state' "${TEMPLATE}"

refresh_held_case 6035.master "${ROOT}/EFFECTIVE/Ra1e10/chinu0.0" 2049 1.0d10 0.0d0
refresh_held_case 6110.master "${ROOT}/EFFECTIVE/Ra1e8/chinu0.8" 513 1.0d8 0.8d0
refresh_held_case 6111.master "${ROOT}/EFFECTIVE/Ra1e8/chinu0.9" 513 1.0d8 0.9d0
