#!/usr/bin/env bash
set -euo pipefail

# Ra=1e10 原始魔法参数极限 chi_nu 短筛选。
# 两个算例均使用 2049x2049 网格，从完全初始状态运行到 100 t_ff。

FLOW_ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
ROOT=${FLOW_ROOT}/ORIGINAL/Ra1e10
STAGE=${FLOW_ROOT}/_staging/20260830_ra1e10_original_chi0995_0999_t100/latest_local.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
PBS_TEMPLATE=${FLOW_ROOT}/FIXED_SQ1/Ra1e10/chinu0.99_t100_rerun_from0_latest/pbs/run.pbs

CASE_0995=${ROOT}/chinu0.995_t100_latest
CASE_0999=${ROOT}/chinu0.999_t100_latest

test -f "${STAGE}"
test -f "${PBS_TEMPLATE}"
for target in "${CASE_0995}" "${CASE_0999}"; do
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_source() {
    local src=$1
    local chi_value=$2

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e 's/^!#define FLOW_ODD_ORIGINAL_MAGIC$/#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^#define FLOW_ODD_FIXED_SQ$/!#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=2049, ny=2049/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi_value}/" \
        -e 's/integer(kind=4), parameter :: loadInitField=1/integer(kind=4), parameter :: loadInitField=0/' \
        -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=100.0d0/' \
        "${src}"

    grep -Fxq '#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq 'parameter :: nx=2049, ny=2049' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d10' "${src}"
    grep -Fq "parameter :: chi_nu=${chi_value}" "${src}"
    grep -Fq 'parameter :: loadInitField=0' "${src}"
    grep -Fq 'parameter :: unsteadyRunDuration=100.0d0' "${src}"
    grep -Fq 'Sq=1.0d0/(0.5d0+flowMagicParameter/(tauf-0.5d0))' "${src}"
}

make_case() {
    local case_dir=$1
    local chi_value=$2
    local job_name=$3
    local node_name=$4
    local old_case=/data2/XLLi/LBMCDE/FLOW-TEST/FIXED_SQ1/Ra1e10/chinu0.99_t100_rerun_from0_latest
    local src=${case_dir}/source/${SOURCE_NAME}
    local result=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${result}"
    configure_source "${src}" "${chi_value}"
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${src}" > "${case_dir}/source/case_source.sha256"

    cp "${PBS_TEMPLATE}" "${pbs}"
    sed -i \
        -e "s|#PBS -N S10C099R100|#PBS -N ${job_name}|" \
        -e "s|#PBS -l nodes=node05:ppn=1|#PBS -l nodes=${node_name}:ppn=1|" \
        -e "s|${old_case}/source/${SOURCE_NAME}|${src}|g" \
        -e "s|${old_case}/results|${result}|g" \
        "${pbs}"
    chmod 700 "${pbs}"

    grep -Fq "#PBS -N ${job_name}" "${pbs}"
    grep -Fq "#PBS -l nodes=${node_name}:ppn=1" "${pbs}"
    grep -Fq "SRC=${src}" "${pbs}"
    grep -Fq "RESULT_DIR=${result}" "${pbs}"

    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job_name}
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=1.0d10
Pr=0.7d0
nx=2049
ny=2049
chi_nu=${chi_value}
flow_policy=FLOW_ODD_ORIGINAL_MAGIC
loadInitField=0
target_time_tf=100.0d0
nonfinite_population_interval_tff=1.0d0
flow_odd_interval_tff=0.1d0
node=${node_name}
purpose=near_unity_chinu_original_magic_stability_screen
EOF
}

make_case "${CASE_0995}" 0.995d0 O10C0995T100 node07
make_case "${CASE_0999}" 0.999d0 O10C0999T100 node05

jid_0995=$(qsub -W depend=afterany:6515.master "${CASE_0995}/pbs/run.pbs")
jid_0999=$(qsub -W depend=afterany:6511.master "${CASE_0999}/pbs/run.pbs")
printf 'node07 case=%s job=%s depends_on=6515.master\n' "${CASE_0995}" "${jid_0995}"
printf 'node05 case=%s job=%s depends_on=6511.master\n' "${CASE_0999}" "${jid_0999}"
