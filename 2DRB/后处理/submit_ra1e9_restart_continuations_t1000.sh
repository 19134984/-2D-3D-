#!/usr/bin/env bash
set -euo pipefail

# Ra=1e9 重启续算：
# 1) 原始魔法参数 chi_nu=0.1、0.99 从 300 t_ff 续算到 1000 t_ff；
# 2) 原始魔法参数 chi_nu=0.95 从发散前的 100 t_ff 快照重启到 300 t_ff；
# 3) chi_nu=0.4、固定 Sq=0.8、0.9 从 100 t_ff 续算到 1000 t_ff，Sq=1 从 300 t_ff 续算到 1000 t_ff。
# 所有任务使用提交时的本地最新源文件；任一任务非数退出后，afterany 链继续下一个算例。

FLOW_ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
STAGE=${FLOW_ROOT}/_staging/20260828_ra1e9_restart_continuations/latest_local.F90
ORIGINAL_ROOT=${FLOW_ROOT}/ORIGINAL/Ra1e9
SQ_ROOT=${FLOW_ROOT}/SQ_SCAN/Ra1e9/chinu0.4
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90

NODE07_DEPEND_AFTER=6360.master
NODE05_DEPEND_AFTER=6366.master

O01_SEED=${ORIGINAL_ROOT}/chinu0.1_t300_latest/results
O095_SEED=${ORIGINAL_ROOT}/chinu0.95_t300_latest/results
O099_SEED=${ORIGINAL_ROOT}/chinu0.99_t300_latest/results
SQ08_SEED=${SQ_ROOT}/sq0.8000000_t100/results
SQ09_SEED=${SQ_ROOT}/sq0.9000000_t100/results
SQ1_SEED=${SQ_ROOT}/sq1.0000000_t300/results

O01_CASE=${ORIGINAL_ROOT}/chinu0.1_restart300_latest_t1000
O095_CASE=${ORIGINAL_ROOT}/chinu0.95_restart100_recheck_t300_latest
O099_CASE=${ORIGINAL_ROOT}/chinu0.99_restart300_latest_t1000
SQ08_CASE=${SQ_ROOT}/sq0.8000000_restart100_latest_t1000
SQ09_CASE=${SQ_ROOT}/sq0.9000000_restart100_latest_t1000
SQ1_CASE=${SQ_ROOT}/sq1.0000000_restart300_latest_t1000

test -f "${STAGE}"
for target in "${O01_CASE}" "${O095_CASE}" "${O099_CASE}" \
              "${SQ08_CASE}" "${SQ09_CASE}" "${SQ1_CASE}"; do
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_common() {
    local src=$1
    local chi=$2
    local target_tf=$3

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=1025, ny=1025/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d9/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi}/" \
        -e 's/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=1/' \
        -e "s/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=${target_tf}/" \
        "${src}"

    grep -Fq 'parameter :: nx=1025, ny=1025' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d9' "${src}"
    grep -Fq "parameter :: chi_nu=${chi}" "${src}"
    grep -Fq 'parameter :: loadInitField=1' "${src}"
    grep -Fq "parameter :: unsteadyRunDuration=${target_tf}" "${src}"
    grep -Fq 'parameter :: nonFiniteCheckIntervalTf=1.0d0' "${src}"
    grep -Fq 'parameter :: flowOddMomentDiagnosticIntervalTf=0.1d0' "${src}"
}

configure_original() {
    local src=$1
    local chi=$2
    local target_tf=$3

    configure_common "${src}" "${chi}" "${target_tf}"
    sed -i \
        -e 's/^!#define FLOW_ODD_ORIGINAL_MAGIC$/#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^#define FLOW_ODD_FIXED_SQ$/!#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        "${src}"
    grep -Fxq '#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq 'Sq=1.0d0/(0.5d0+flowMagicParameter/(tauf-0.5d0))' "${src}"
}

configure_fixed_sq() {
    local src=$1
    local sq=$2
    local target_tf=$3

    configure_common "${src}" 0.4d0 "${target_tf}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        "${src}"
    sed -i "/^#ifdef FLOW_ODD_FIXED_SQ$/,/^#endif$/ s|^[[:space:]]*real(kind=8), parameter :: Sq=.*$|        real(kind=8), parameter :: Sq=${sq} ! 固定奇矩松弛率诊断|" "${src}"
    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq "parameter :: Sq=${sq}" "${src}"
}

write_pbs() {
    local case_dir=$1
    local seed_dir=$2
    local expected_seed_tf=$3
    local expected_seed_count=$4
    local accepted_seed_status=$5
    local job_name=$6
    local node_name=$7
    local src=${case_dir}/source/${SOURCE_NAME}
    local result=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs

    cat > "${pbs}" <<'PBSEOF'
#!/usr/bin/env bash
#PBS -N __JOB__
#PBS -q batch
#PBS -l nodes=__NODE__:ppn=1
#PBS -l walltime=96:00:00
#PBS -o __RESULT__/pbs.stdout
#PBS -e __RESULT__/pbs.stderr

set -u
SRC=__SOURCE__
SEED_DIR=__SEED__
EXPECTED_SEED_TF=__SEED_TF__
EXPECTED_SEED_COUNT=__SEED_COUNT__
ACCEPTED_SEED_STATUS=__SEED_STATUS__
RESULT_DIR=__RESULT__
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=${NVHPC_ROOT}/compilers/bin/nvfortran
RELOAD_PREFIX=reloadFile2DOpenaccLBMCDE_D2Q5

mkdir -p "${RESULT_DIR}"
cd "${RESULT_DIR}" || exit 90

if [[ ! -f "${SEED_DIR}/run.status" ]] || \
   [[ "$(<"${SEED_DIR}/run.status")" != "${ACCEPTED_SEED_STATUS}" ]]; then
    printf 'Unexpected restart seed run.status: %s\n' "${SEED_DIR}" > restart_preflight.error
    printf '121\n' > run.status
    exit 121
fi

meta=${SEED_DIR}/${RELOAD_PREFIX}-latest.meta
if [[ ! -f "${meta}" ]]; then
    printf 'Missing restart metadata: %s\n' "${meta}" > restart_preflight.error
    printf '122\n' > run.status
    exit 122
fi
if ! awk -v target="${EXPECTED_SEED_TF}" '$1=="time_tf" {t=$2} END {exit !(t>target-0.01 && t<target+0.01)}' "${meta}"; then
    printf 'Unexpected restart seed time: %s\n' "${meta}" > restart_preflight.error
    printf '122\n' > run.status
    exit 122
fi
if ! grep -Fq "cumulativeStatisticSampleCount ${EXPECTED_SEED_COUNT}" "${meta}"; then
    printf 'Unexpected restart seed sample count: %s\n' "${meta}" > restart_preflight.error
    printf '123\n' > run.status
    exit 123
fi

required_histories=(
    Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat
    Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat
    DissipationHistory_2DOpenaccLBMCDE.dat
    TemperatureProfileHistory_2DOpenaccLBMCDE.bin
    PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q5.dat
    FlowOddMomentHistory_2DOpenaccLBMCDE_D2Q5.dat
)
for required in "${required_histories[@]}"; do
    if [[ ! -f "${SEED_DIR}/${required}" ]]; then
        printf 'Missing restart history: %s\n' "${SEED_DIR}/${required}" > restart_preflight.error
        printf '124\n' > run.status
        exit 124
    fi
    cp -a "${SEED_DIR}/${required}" "${RESULT_DIR}/"
done
cp -a "${meta}" "${RESULT_DIR}/${RELOAD_PREFIX}-latest.meta"

reload_id=$(awk '$1=="reloadFileName" {print $2}' "${meta}")
restart_bin=${SEED_DIR}/${RELOAD_PREFIX}-${reload_id}.bin
if [[ -z "${reload_id}" ]] || [[ ! -f "${restart_bin}" ]]; then
    printf 'Missing referenced restart binary: %s\n' "${restart_bin}" > restart_preflight.error
    printf '124\n' > run.status
    exit 124
fi
cp -a "${restart_bin}" "${RESULT_DIR}/"
sha256sum "${meta}" "${restart_bin}" > input_restart.sha256
cp "${SEED_DIR}/run.status" input_seed_run.status
if [[ -f "${SEED_DIR}/runtime_source.sha256" ]]; then
    cp "${SEED_DIR}/runtime_source.sha256" input_seed_runtime_source.sha256
fi

export PATH="${NVHPC_ROOT}/compilers/bin:${PATH}"
export LD_LIBRARY_PATH="${NVHPC_ROOT}/compilers/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1
export ACC_DEVICE_TYPE=nvidia

printf 'PBS_JOBID=%s\n' "${PBS_JOBID:-unknown}" > job_identity.txt
printf 'HOSTNAME=%s\n' "$(hostname)" >> job_identity.txt
date '+CASE_START %F %T %z' > timing.txt
sha256sum "${SRC}" > runtime_source.sha256
"${NVFORTRAN}" --version > compiler_version.txt 2>&1
nvidia-smi > nvidia_smi_before.log 2>&1 || true

BUILD_FLAGS=(-cpp -O3 -acc -gpu=cc60 -Minfo=accel)
printf '%q ' "${BUILD_FLAGS[@]}" > build_flags.txt
printf '\n' >> build_flags.txt
"${NVFORTRAN}" "${BUILD_FLAGS[@]}" "${SRC}" -o solver.exe > compile.log 2>&1
compile_status=$?
printf '%s\n' "${compile_status}" > compile.status
if [[ "${compile_status}" -ne 0 ]]; then
    printf '125\n' > run.status
    date '+CASE_END %F %T %z' >> timing.txt
    exit "${compile_status}"
fi

set +e
./solver.exe > solver.stdout 2>&1
run_status=$?
set -e
if [[ "${run_status}" -eq 0 ]] && grep -Eqi \
    '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)' \
    solver.stdout Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat \
    Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat DissipationHistory_2DOpenaccLBMCDE.dat \
    PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q5.dat \
    FlowOddMomentHistory_2DOpenaccLBMCDE_D2Q5.dat 2>/dev/null; then
    printf 'Non-finite value found in text output/history.\n' >> solver.stdout
    run_status=86
fi
printf '%s\n' "${run_status}" > run.status
nvidia-smi > nvidia_smi_after.log 2>&1 || true
date '+CASE_END %F %T %z' >> timing.txt
exit "${run_status}"
PBSEOF

    sed -i \
        -e "s|__JOB__|${job_name}|g" \
        -e "s|__NODE__|${node_name}|g" \
        -e "s|__SOURCE__|${src}|g" \
        -e "s|__SEED__|${seed_dir}|g" \
        -e "s|__SEED_TF__|${expected_seed_tf}|g" \
        -e "s|__SEED_COUNT__|${expected_seed_count}|g" \
        -e "s|__SEED_STATUS__|${accepted_seed_status}|g" \
        -e "s|__RESULT__|${result}|g" \
        "${pbs}"
    chmod 700 "${pbs}"
}

create_case() {
    local case_dir=$1
    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
}

write_provenance() {
    local case_dir=$1
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${case_dir}/source/${SOURCE_NAME}" > "${case_dir}/source/case_source.sha256"
}

make_original_case() {
    local case_dir=$1
    local chi=$2
    local target_tf=$3
    local seed_dir=$4
    local seed_tf=$5
    local seed_count=$6
    local seed_status=$7
    local job_name=$8
    local node_name=$9

    create_case "${case_dir}"
    configure_original "${case_dir}/source/${SOURCE_NAME}" "${chi}" "${target_tf}"
    write_provenance "${case_dir}"
    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job_name}
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=1.0d9
Pr=0.7d0
nx=1025
ny=1025
chi_nu=${chi}
flow_policy=FLOW_ODD_ORIGINAL_MAGIC
flow_magic_parameter=3/16
loadInitField=1
target_time_tf=${target_tf}
restart_seed=${seed_dir}
restart_seed_time_tf=${seed_tf}
accepted_seed_run_status=${seed_status}
nonfinite_population_interval_tff=1.0d0
flow_odd_interval_tff=0.1d0
node=${node_name}
purpose=restart_stability_and_convergence_check
EOF
    write_pbs "${case_dir}" "${seed_dir}" "${seed_tf}" "${seed_count}" \
        "${seed_status}" "${job_name}" "${node_name}"
}

make_fixed_case() {
    local case_dir=$1
    local sq=$2
    local target_tf=$3
    local seed_dir=$4
    local seed_tf=$5
    local seed_count=$6
    local job_name=$7
    local node_name=$8

    create_case "${case_dir}"
    configure_fixed_sq "${case_dir}/source/${SOURCE_NAME}" "${sq}" "${target_tf}"
    write_provenance "${case_dir}"
    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job_name}
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=1.0d9
Pr=0.7d0
nx=1025
ny=1025
chi_nu=0.4d0
flow_policy=FLOW_ODD_FIXED_SQ
Sq=${sq}
loadInitField=1
target_time_tf=${target_tf}
restart_seed=${seed_dir}
restart_seed_time_tf=${seed_tf}
accepted_seed_run_status=0
nonfinite_population_interval_tff=1.0d0
flow_odd_interval_tff=0.1d0
node=${node_name}
purpose=restart_stability_and_convergence_check
EOF
    write_pbs "${case_dir}" "${seed_dir}" "${seed_tf}" "${seed_count}" \
        0 "${job_name}" "${node_name}"
}

make_original_case "${O01_CASE}" 0.1d0 1000.0d0 "${O01_SEED}" 300 301 0 O9C01T1000 node07
make_original_case "${O095_CASE}" 0.95d0 300.0d0 "${O095_SEED}" 100 101 86 O9C095R300 node07
make_original_case "${O099_CASE}" 0.99d0 1000.0d0 "${O099_SEED}" 300 301 0 O9C099T1000 node07

make_fixed_case "${SQ08_CASE}" 0.8d0 1000.0d0 "${SQ08_SEED}" 100 101 S9C4Q8T1000 node05
make_fixed_case "${SQ09_CASE}" 0.9d0 1000.0d0 "${SQ09_SEED}" 100 101 S9C4Q9T1000 node05
make_fixed_case "${SQ1_CASE}" 1.0d0 1000.0d0 "${SQ1_SEED}" 300 301 S9C4Q1T1000 node05

previous=${NODE07_DEPEND_AFTER}
for case_dir in "${O01_CASE}" "${O095_CASE}" "${O099_CASE}"; do
    jid=$(qsub -W depend=afterany:${previous} "${case_dir}/pbs/run.pbs")
    printf 'node07 case=%s job=%s depends_on=%s\n' "${case_dir}" "${jid}" "${previous}"
    previous=${jid}
done

previous=${NODE05_DEPEND_AFTER}
for case_dir in "${SQ08_CASE}" "${SQ09_CASE}" "${SQ1_CASE}"; do
    jid=$(qsub -W depend=afterany:${previous} "${case_dir}/pbs/run.pbs")
    printf 'node05 case=%s job=%s depends_on=%s\n' "${case_dir}" "${jid}" "${previous}"
    previous=${jid}
done
