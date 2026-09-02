#!/usr/bin/env bash
set -euo pipefail

# Ra=1e9 固定 Sq 控制算例：从 100 t_ff 续算到 1400 t_ff。
# 目标统计总窗为 700--1400 t_ff，最终结果取 1050--1400 t_ff。

FLOW_ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
ROOT=${FLOW_ROOT}/SQ_SCAN/Ra1e9
STAGE=${FLOW_ROOT}/_staging/20260831_ra1e9_selected_sq_t1400/latest_local.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90

TEMPLATE_CASE=${ROOT}/chinu0.4/sq1.0000000_restart1000_latest_t1400
PBS_TEMPLATE=${TEMPLATE_CASE}/pbs/run.pbs
TEMPLATE_SEED=${ROOT}/chinu0.4/sq1.0000000_restart300_latest_t1000/results

C05Q1_SEED=${ROOT}/chinu0.5/sq1.0000000_t100_latest/results
C09Q05_SEED=${ROOT}/chinu0.9/sq0.5000000_t100_latest/results
C09Q1_SEED=${ROOT}/chinu0.9/sq1.0000000_t100_latest/results
C05Q1_CASE=${ROOT}/chinu0.5/sq1.0000000_restart100_latest_t1400
C09Q05_CASE=${ROOT}/chinu0.9/sq0.5000000_restart100_latest_t1400
C09Q1_CASE=${ROOT}/chinu0.9/sq1.0000000_restart100_latest_t1400

test -f "${STAGE}"
test -f "${PBS_TEMPLATE}"
for target in "${C05Q1_CASE}" "${C09Q05_CASE}" "${C09Q1_CASE}"; do
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_source() {
    local src=$1
    local chi_value=$2
    local sq_value=$3

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=1025, ny=1025/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d9/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi_value}/" \
        -e 's/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=1/' \
        -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=1400.0d0/' \
        "${src}"
    sed -i "/^#ifdef FLOW_ODD_FIXED_SQ$/,/^#endif$/ s|^[[:space:]]*real(kind=8), parameter :: Sq=.*$|        real(kind=8), parameter :: Sq=${sq_value} ! 固定奇矩松弛率控制实验|" "${src}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq 'parameter :: nx=1025, ny=1025' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d9' "${src}"
    grep -Fq "parameter :: chi_nu=${chi_value}" "${src}"
    grep -Fq "parameter :: Sq=${sq_value}" "${src}"
    grep -Fq 'parameter :: loadInitField=1' "${src}"
    grep -Fq 'parameter :: unsteadyRunDuration=1400.0d0' "${src}"
}

make_case() {
    local case_dir=$1
    local seed_dir=$2
    local chi_value=$3
    local sq_value=$4
    local job_name=$5
    local node_name=$6
    local src=${case_dir}/source/${SOURCE_NAME}
    local result=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${result}"
    configure_source "${src}" "${chi_value}" "${sq_value}"
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${src}" > "${case_dir}/source/case_source.sha256"

    cp "${PBS_TEMPLATE}" "${pbs}"
    sed -i \
        -e "s|#PBS -N S9C4Q1T1400|#PBS -N ${job_name}|" \
        -e "s|#PBS -l nodes=node07:ppn=1|#PBS -l nodes=${node_name}:ppn=1|" \
        -e "s|${TEMPLATE_CASE}/source/${SOURCE_NAME}|${src}|g" \
        -e "s|${TEMPLATE_CASE}/results|${result}|g" \
        -e "s|${TEMPLATE_SEED}|${seed_dir}|g" \
        -e 's/EXPECTED_SEED_TF=1000/EXPECTED_SEED_TF=100/' \
        -e 's/EXPECTED_SEED_COUNT=1001/EXPECTED_SEED_COUNT=101/' \
        "${pbs}"
    chmod 700 "${pbs}"

    grep -Fq "#PBS -N ${job_name}" "${pbs}"
    grep -Fq "#PBS -l nodes=${node_name}:ppn=1" "${pbs}"
    grep -Fq "SRC=${src}" "${pbs}"
    grep -Fq "SEED_DIR=${seed_dir}" "${pbs}"
    grep -Fq 'EXPECTED_SEED_TF=100' "${pbs}"
    grep -Fq 'EXPECTED_SEED_COUNT=101' "${pbs}"
    grep -Fq "RESULT_DIR=${result}" "${pbs}"

    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job_name}
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=1.0d9
Pr=0.7d0
nx=1025
ny=1025
chi_nu=${chi_value}
flow_policy=FLOW_ODD_FIXED_SQ
Sq=${sq_value}
loadInitField=1
target_time_tf=1400.0d0
restart_seed=${seed_dir}
restart_seed_time_tf=100
restart_seed_sample_count=101
statistics_window_tf=700--1400
final_statistics_window_tf=1050--1400
nonfinite_population_interval_tff=1.0d0
flow_odd_interval_tff=0.1d0
node=${node_name}
EOF
}

make_case "${C05Q1_CASE}" "${C05Q1_SEED}" 0.5d0 1.0d0 S9C05Q1T1400 node05
make_case "${C09Q05_CASE}" "${C09Q05_SEED}" 0.9d0 0.5d0 S9C09Q05T1400 node05
make_case "${C09Q1_CASE}" "${C09Q1_SEED}" 0.9d0 1.0d0 S9C09Q1T1400 node07

jid_c05q1=$(qsub "${C05Q1_CASE}/pbs/run.pbs")
jid_c09q05=$(qsub -W depend=afterany:${jid_c05q1} "${C09Q05_CASE}/pbs/run.pbs")
jid_c09q1=$(qsub -W depend=afterany:6546.master "${C09Q1_CASE}/pbs/run.pbs")

printf 'node05 case=%s job=%s depends_on=none\n' "${C05Q1_CASE}" "${jid_c05q1}"
printf 'node05 case=%s job=%s depends_on=%s\n' "${C09Q05_CASE}" "${jid_c09q05}" "${jid_c05q1}"
printf 'node07 case=%s job=%s depends_on=6546.master\n' "${C09Q1_CASE}" "${jid_c09q1}"
