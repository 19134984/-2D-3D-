#!/usr/bin/env bash
set -euo pipefail

# 补齐两处重启时间容差后，只重提 Sq=0.05 和 0.08；Sq=0.10 保持原任务 6269 不动。
ROOT=/data2/XLLi/LBMCDE/FLOW-TEST/SQ_SCAN/Ra1e8/chinu0.4
STAGED_SOURCE="${ROOT}/_staging/first_solver_restart_tolerance_both_fixed.F90"
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
SEED_005="${ROOT}/sq0.0500000_v2/results"
SEED_008="${ROOT}/sq0.0800000/results"
CASE_005="${ROOT}/sq0.0500000_t1000_v2"
CASE_008="${ROOT}/sq0.0800000_t1000_v3"

if [[ ! -f "${STAGED_SOURCE}" ]]; then
    printf 'Missing staged source: %s\n' "${STAGED_SOURCE}" >&2
    exit 2
fi
if [[ -e "${CASE_005}" ]] || [[ -e "${CASE_008}" ]]; then
    printf 'Refusing to overwrite an existing retry directory.\n' >&2
    exit 3
fi

configure_source() {
    local source_file=$1
    local sq_value=$2
    cp "${STAGED_SOURCE}" "${source_file}"
    sed -i 's/\r$//' "${source_file}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e "s/real(kind=8), parameter :: Sq=1.0d0 ! 固定奇模态松弛率对照：tau_q=1/real(kind=8), parameter :: Sq=${sq_value} ! 固定奇模态松弛率长时间测试/" \
        -e 's/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=1/' \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=513, ny=513/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d8/' \
        -e 's/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=0.4d0/' \
        "${source_file}"

    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${source_file}"
    grep -Fq "parameter :: Sq=${sq_value}" "${source_file}"
    grep -Fq 'parameter :: loadInitField=1' "${source_file}"
    [[ "$(grep -Fc '0.5d0/timeUnit+64.0d0*epsilon(1.0d0)' "${source_file}")" -eq 2 ]]
}

make_retry() {
    local case_dir=$1
    local seed_dir=$2
    local old_case_dir=$3
    local old_job_name=$4
    local new_job_name=$5
    local sq_value=$6
    local keep_diagnostics=$7

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    cp -a "${seed_dir}/." "${case_dir}/results/"
    cp "${seed_dir}/run.status" "${case_dir}/results/seed_300tff_run.status"
    cp "${seed_dir}/timing.txt" "${case_dir}/results/seed_300tff_timing.txt"
    cp "${seed_dir}/runtime_source.sha256" "${case_dir}/results/seed_300tff_runtime_source.sha256"
    if [[ "${keep_diagnostics}" -eq 0 ]]; then
        rm -f "${case_dir}/results/PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q5.dat"
        rm -f "${case_dir}/results/FlowOddMomentHistory_2DOpenaccLBMCDE_D2Q5.dat"
    fi

    configure_source "${case_dir}/source/${SOURCE_NAME}" "${sq_value}"
    sha256sum "${STAGED_SOURCE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${case_dir}/source/${SOURCE_NAME}" > "${case_dir}/source/case_source.sha256"
    {
        printf 'node=node07\nRa=1.0d8\nPr=0.7d0\nnx=513\nny=513\nchi_nu=0.4d0\nchi_b=0.0d0\n'
        printf 'flow_policy=FLOW_ODD_FIXED_SQ\nSq=%s\n' "${sq_value}"
        printf 'unsteadyRunDuration=1000.0d0\nloadInitField=1\nrestart_seed=%s\n' "${seed_dir}"
        printf 'restart_tolerance_fix=NuRe_and_dissipation_profile\n'
    } > "${case_dir}/case_settings.txt"

    cp "${old_case_dir}/pbs/run.pbs" "${case_dir}/pbs/run.pbs"
    sed -i \
        -e "s|${old_case_dir}|${case_dir}|g" \
        -e "s/^#PBS -N ${old_job_name}$/#PBS -N ${new_job_name}/" \
        "${case_dir}/pbs/run.pbs"
    chmod 750 "${case_dir}/pbs/run.pbs"
}

make_retry "${CASE_005}" "${SEED_005}" "${ROOT}/sq0.0500000_t1000" SQ8C05T1000 SQ8C05T1V2 0.05d0 1
make_retry "${CASE_008}" "${SEED_008}" "${ROOT}/sq0.0800000_t1000_v2" SQ8C08T1V2 SQ8C08T1V3 0.08d0 0

sq005_job=$(qsub -W depend=afterany:6269.master "${CASE_005}/pbs/run.pbs")
sq008_job=$(qsub -W depend=afterany:${sq005_job} "${CASE_008}/pbs/run.pbs")

{
    printf 'sq0.0500000_t1000_v2=%s\n' "${sq005_job}"
    printf 'sq0.0800000_t1000_v3=%s\n' "${sq008_job}"
} | tee "${ROOT}/restart_retry_submission_manifest.txt"
