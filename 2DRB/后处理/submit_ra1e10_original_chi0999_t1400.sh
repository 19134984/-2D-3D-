#!/usr/bin/env bash
set -euo pipefail

# Ra=1e10、原始魔法参数、chi_nu=0.999：从 100 续算到 1400 t_ff。

FLOW_ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
STAGE=${FLOW_ROOT}/_staging/20260831_ra1e10_original_chi0999_t1400/latest_local.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
SEED=${FLOW_ROOT}/ORIGINAL/Ra1e10/chinu0.999_t100_latest/results
CASE=${FLOW_ROOT}/ORIGINAL/Ra1e10/chinu0.999_restart100_latest_t1400
TEMPLATE_CASE=${FLOW_ROOT}/FIXED_SQ1/Ra1e10/chinu0.8_restart100_latest_t300
PBS_TEMPLATE=${TEMPLATE_CASE}/pbs/run.pbs
TEMPLATE_SEED=${FLOW_ROOT}/FIXED_SQ1/Ra1e10/chinu0.8_t100_latest/results

test -f "${STAGE}"
test -f "${PBS_TEMPLATE}"
if [[ -e "${CASE}" ]]; then
    printf 'Refusing to overwrite %s\n' "${CASE}" >&2
    exit 3
fi

mkdir -p "${CASE}/source" "${CASE}/pbs" "${CASE}/results"
SRC=${CASE}/source/${SOURCE_NAME}
RESULT=${CASE}/results
PBS=${CASE}/pbs/run.pbs

cp "${STAGE}" "${SRC}"
sed -i 's/\r$//' "${SRC}"
sed -i \
    -e 's/^!#define FLOW_ODD_ORIGINAL_MAGIC$/#define FLOW_ODD_ORIGINAL_MAGIC/' \
    -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
    -e 's/^#define FLOW_ODD_FIXED_SQ$/!#define FLOW_ODD_FIXED_SQ/' \
    -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
    -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=2049, ny=2049/' \
    -e 's/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=0.999d0/' \
    -e 's/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=1/' \
    -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=1400.0d0/' \
    "${SRC}"

grep -Fxq '#define FLOW_ODD_ORIGINAL_MAGIC' "${SRC}"
grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${SRC}"
grep -Fxq '!#define FLOW_ODD_FIXED_SQ' "${SRC}"
grep -Fxq '!#define FLOW_BGK' "${SRC}"
grep -Fq 'parameter :: nx=2049, ny=2049' "${SRC}"
grep -Fq 'parameter :: Rayleigh=1.0d10' "${SRC}"
grep -Fq 'parameter :: chi_nu=0.999d0' "${SRC}"
grep -Fq 'parameter :: loadInitField=1' "${SRC}"
grep -Fq 'parameter :: unsteadyRunDuration=1400.0d0' "${SRC}"
grep -Fq 'Sq=1.0d0/(0.5d0+flowMagicParameter/(tauf-0.5d0))' "${SRC}"

sha256sum "${STAGE}" > "${CASE}/source/local_master.sha256"
sha256sum "${SRC}" > "${CASE}/source/case_source.sha256"

cp "${PBS_TEMPLATE}" "${PBS}"
sed -i \
    -e 's|#PBS -N S10C08T300|#PBS -N O10C0999T1400|' \
    -e "s|${TEMPLATE_CASE}/source/${SOURCE_NAME}|${SRC}|g" \
    -e "s|${TEMPLATE_CASE}/results|${RESULT}|g" \
    -e "s|${TEMPLATE_SEED}|${SEED}|g" \
    "${PBS}"
chmod 700 "${PBS}"

grep -Fq '#PBS -N O10C0999T1400' "${PBS}"
grep -Fq '#PBS -l nodes=node05:ppn=1' "${PBS}"
grep -Fq "SRC=${SRC}" "${PBS}"
grep -Fq "SEED_DIR=${SEED}" "${PBS}"
grep -Fq 'EXPECTED_SEED_TF=100' "${PBS}"
grep -Fq 'EXPECTED_SEED_COUNT=101' "${PBS}"
grep -Fq "RESULT_DIR=${RESULT}" "${PBS}"

cat > "${CASE}/case_settings.txt" <<EOF
case=O10C0999T1400
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=1.0d10
Pr=0.7d0
nx=2049
ny=2049
chi_nu=0.999d0
tau0=3.469283095d0
flow_policy=FLOW_ODD_ORIGINAL_MAGIC
Sq=1.775736689d0
loadInitField=1
target_time_tf=1400.0d0
restart_seed=${SEED}
restart_seed_time_tf=100
restart_seed_sample_count=101
statistics_window_tf=700--1400
final_statistics_window_tf=1050--1400
nonfinite_population_interval_tff=1.0d0
flow_odd_interval_tff=0.1d0
node=node05
purpose=near_unity_chinu_delayed_instability_and_long_statistics_test
EOF

jid=$(qsub -W depend=afterany:6586.master "${PBS}")
printf 'node05 case=%s job=%s depends_on=6586.master\n' "${CASE}" "${jid}"
