#!/usr/bin/env bash
set -euo pipefail

# 第一套 D2Q9-TRT 流场 + D2Q5-Luo-TRT 温度场后续稳定性测试：
# 1. 固定 Sq=1、chi_nu=0.8/0.9/0.95，从 100 续算到 300 t_ff；
# 2. 固定 Sq=1、chi_nu=0.99，从完全初始状态重新计算到 100 t_ff；
# 3. Ra=1e9、chi_nu=0.5/0.9、固定 Sq=0.5/1，从完全初始状态计算到 100 t_ff。
# 每个任务检测到非数后退出；同节点后续任务使用 afterany，仍会继续执行。

FLOW_ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
STAGE=${FLOW_ROOT}/_staging/20260830_ra1e9_ra1e10_sq_followup/latest_local.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90

SQ10_ROOT=${FLOW_ROOT}/FIXED_SQ1/Ra1e10
SQ08_SEED=${SQ10_ROOT}/chinu0.8_t100_latest/results
SQ09_SEED=${SQ10_ROOT}/chinu0.9_t100_latest/results
SQ095_SEED=${SQ10_ROOT}/chinu0.95_t100_latest/results
SQ08_CASE=${SQ10_ROOT}/chinu0.8_restart100_latest_t300
SQ09_CASE=${SQ10_ROOT}/chinu0.9_restart100_latest_t300
SQ095_CASE=${SQ10_ROOT}/chinu0.95_restart100_latest_t300
SQ099_RERUN_CASE=${SQ10_ROOT}/chinu0.99_t100_rerun_from0_latest
SQ9_ROOT=${FLOW_ROOT}/SQ_SCAN/Ra1e9
SQ9_C05_Q05_CASE=${SQ9_ROOT}/chinu0.5/sq0.5000000_t100_latest
SQ9_C05_Q1_CASE=${SQ9_ROOT}/chinu0.5/sq1.0000000_t100_latest
SQ9_C09_Q05_CASE=${SQ9_ROOT}/chinu0.9/sq0.5000000_t100_latest
SQ9_C09_Q1_CASE=${SQ9_ROOT}/chinu0.9/sq1.0000000_t100_latest

test -f "${STAGE}"
for target in "${SQ08_CASE}" "${SQ09_CASE}" "${SQ095_CASE}" \
    "${SQ099_RERUN_CASE}" "${SQ9_C05_Q05_CASE}" "${SQ9_C05_Q1_CASE}" \
    "${SQ9_C09_Q05_CASE}" "${SQ9_C09_Q1_CASE}"; do
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_common() {
    local src=$1
    local nx_value=$2
    local ra_value=$3
    local chi_value=$4
    local target_tf=$5
    local load_value=$6

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e "s/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=${nx_value}, ny=${nx_value}/" \
        -e "s/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=${ra_value}/" \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi_value}/" \
        -e "s/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=${load_value}/" \
        -e "s/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=${target_tf}/" \
        "${src}"

    grep -Fq "parameter :: nx=${nx_value}, ny=${nx_value}" "${src}"
    grep -Fq "parameter :: Rayleigh=${ra_value}" "${src}"
    grep -Fq "parameter :: chi_nu=${chi_value}" "${src}"
    grep -Fq "parameter :: loadInitField=${load_value}" "${src}"
    grep -Fq "parameter :: unsteadyRunDuration=${target_tf}" "${src}"
    grep -Fq 'parameter :: nonFiniteCheckIntervalTf=1.0d0' "${src}"
    grep -Fq 'parameter :: flowOddMomentDiagnosticIntervalTf=0.1d0' "${src}"
}

configure_fixed_sq() {
    local src=$1
    local nx_value=$2
    local ra_value=$3
    local chi_value=$4
    local target_tf=$5
    local load_value=$6
    local sq_value=$7

    configure_common "${src}" "${nx_value}" "${ra_value}" "${chi_value}" "${target_tf}" "${load_value}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        "${src}"
    sed -i "/^#ifdef FLOW_ODD_FIXED_SQ$/,/^#endif$/ s|^[[:space:]]*real(kind=8), parameter :: Sq=.*$|        real(kind=8), parameter :: Sq=${sq_value} ! 固定奇矩松弛率对照|" "${src}"
    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq "parameter :: Sq=${sq_value}" "${src}"
}

write_runtime_body() {
    local pbs=$1
    cat >> "${pbs}" <<'PBSEOF'

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
}

write_restart_pbs() {
    local case_dir=$1
    local seed_dir=$2
    local seed_tf=$3
    local seed_count=$4
    local job_name=$5
    local node_name=$6
    local src=${case_dir}/source/${SOURCE_NAME}
    local result=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs

    cat > "${pbs}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=${node_name}:ppn=1
#PBS -l walltime=96:00:00
#PBS -o ${result}/pbs.stdout
#PBS -e ${result}/pbs.stderr

set -u
SRC=${src}
SEED_DIR=${seed_dir}
EXPECTED_SEED_TF=${seed_tf}
EXPECTED_SEED_COUNT=${seed_count}
RESULT_DIR=${result}
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=\${NVHPC_ROOT}/compilers/bin/nvfortran
RELOAD_PREFIX=reloadFile2DOpenaccLBMCDE_D2Q5

mkdir -p "\${RESULT_DIR}"
cd "\${RESULT_DIR}" || exit 90

if [[ ! -f "\${SEED_DIR}/run.status" ]] || [[ "\$(<"\${SEED_DIR}/run.status")" != 0 ]]; then
    printf 'Invalid restart seed run.status: %s\n' "\${SEED_DIR}" > restart_preflight.error
    printf '121\n' > run.status
    exit 121
fi
meta=\${SEED_DIR}/\${RELOAD_PREFIX}-latest.meta
if [[ ! -f "\${meta}" ]]; then
    printf 'Missing restart metadata: %s\n' "\${meta}" > restart_preflight.error
    printf '122\n' > run.status
    exit 122
fi
if ! awk -v target="\${EXPECTED_SEED_TF}" '\$1=="time_tf" {t=\$2} END {exit !(t>target-0.01 && t<target+0.01)}' "\${meta}"; then
    printf 'Unexpected restart seed time: %s\n' "\${meta}" > restart_preflight.error
    printf '122\n' > run.status
    exit 122
fi
if ! grep -Fq "cumulativeStatisticSampleCount \${EXPECTED_SEED_COUNT}" "\${meta}"; then
    printf 'Unexpected restart seed sample count: %s\n' "\${meta}" > restart_preflight.error
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
for required in "\${required_histories[@]}"; do
    if [[ ! -f "\${SEED_DIR}/\${required}" ]]; then
        printf 'Missing restart history: %s\n' "\${SEED_DIR}/\${required}" > restart_preflight.error
        printf '124\n' > run.status
        exit 124
    fi
    cp -a "\${SEED_DIR}/\${required}" "\${RESULT_DIR}/"
done
cp -a "\${meta}" "\${RESULT_DIR}/\${RELOAD_PREFIX}-latest.meta"

reload_id=\$(awk '\$1=="reloadFileName" {print \$2}' "\${meta}")
restart_bin=\${SEED_DIR}/\${RELOAD_PREFIX}-\${reload_id}.bin
if [[ -z "\${reload_id}" ]] || [[ ! -f "\${restart_bin}" ]]; then
    printf 'Missing referenced restart binary: %s\n' "\${restart_bin}" > restart_preflight.error
    printf '124\n' > run.status
    exit 124
fi
cp -a "\${restart_bin}" "\${RESULT_DIR}/"
sha256sum "\${meta}" "\${restart_bin}" > input_restart.sha256
cp "\${SEED_DIR}/run.status" input_seed_run.status
if [[ -f "\${SEED_DIR}/runtime_source.sha256" ]]; then
    cp "\${SEED_DIR}/runtime_source.sha256" input_seed_runtime_source.sha256
fi
EOF
    write_runtime_body "${pbs}"
    chmod 700 "${pbs}"
}

write_scratch_pbs() {
    local case_dir=$1
    local job_name=$2
    local node_name=$3
    local src=${case_dir}/source/${SOURCE_NAME}
    local result=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs

    cat > "${pbs}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=${node_name}:ppn=1
#PBS -l walltime=96:00:00
#PBS -o ${result}/pbs.stdout
#PBS -e ${result}/pbs.stderr

set -u
SRC=${src}
RESULT_DIR=${result}
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=\${NVHPC_ROOT}/compilers/bin/nvfortran

mkdir -p "\${RESULT_DIR}"
cd "\${RESULT_DIR}" || exit 90
EOF
    write_runtime_body "${pbs}"
    chmod 700 "${pbs}"
}

create_case_dirs() {
    local case_dir=$1
    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
}

record_source() {
    local case_dir=$1
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${case_dir}/source/${SOURCE_NAME}" > "${case_dir}/source/case_source.sha256"
}

make_sq10_restart() {
    local chi_label=$1
    local chi_value=$2
    local seed_dir=$3
    local case_dir=$4
    local job_name=$5

    create_case_dirs "${case_dir}"
    configure_fixed_sq "${case_dir}/source/${SOURCE_NAME}" 2049 1.0d10 "${chi_value}" 300.0d0 1 1.0d0
    record_source "${case_dir}"
    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job_name}
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=1.0d10
Pr=0.7d0
nx=2049
ny=2049
chi_nu=${chi_value}
flow_policy=FLOW_ODD_FIXED_SQ
Sq=1.0d0
loadInitField=1
target_time_tf=300.0d0
restart_seed=${seed_dir}
restart_seed_time_tf=100
restart_seed_sample_count=101
nonfinite_population_interval_tff=1.0d0
flow_odd_interval_tff=0.1d0
node=node05
purpose=delayed_instability_screen
EOF
    write_restart_pbs "${case_dir}" "${seed_dir}" 100 101 "${job_name}" node05
}

make_fixed_sq_scratch() {
    local case_dir=$1
    local nx_value=$2
    local ra_value=$3
    local chi_value=$4
    local sq_value=$5
    local job_name=$6
    local node_name=$7
    local purpose=$8

    create_case_dirs "${case_dir}"
    configure_fixed_sq "${case_dir}/source/${SOURCE_NAME}" "${nx_value}" "${ra_value}" \
        "${chi_value}" 100.0d0 0 "${sq_value}"
    record_source "${case_dir}"
    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job_name}
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=${ra_value}
Pr=0.7d0
nx=${nx_value}
ny=${nx_value}
chi_nu=${chi_value}
flow_policy=FLOW_ODD_FIXED_SQ
Sq=${sq_value}
loadInitField=0
target_time_tf=100.0d0
nonfinite_population_interval_tff=1.0d0
flow_odd_interval_tff=0.1d0
node=${node_name}
purpose=${purpose}
EOF
    write_scratch_pbs "${case_dir}" "${job_name}" "${node_name}"
}

make_sq10_restart 0.8 0.8d0 "${SQ08_SEED}" "${SQ08_CASE}" S10C08T300
make_sq10_restart 0.9 0.9d0 "${SQ09_SEED}" "${SQ09_CASE}" S10C09T300
make_sq10_restart 0.95 0.95d0 "${SQ095_SEED}" "${SQ095_CASE}" S10C095T300
make_fixed_sq_scratch "${SQ099_RERUN_CASE}" 2049 1.0d10 0.99d0 1.0d0 S10C099R100 node05 reproduce_chinu099_early_failure_from_identical_initial_state
make_fixed_sq_scratch "${SQ9_C05_Q05_CASE}" 1025 1.0d9 0.5d0 0.5d0 S9C05Q05 node07 fixed_sq_stability_control
make_fixed_sq_scratch "${SQ9_C05_Q1_CASE}" 1025 1.0d9 0.5d0 1.0d0 S9C05Q1 node07 fixed_sq_stability_control
make_fixed_sq_scratch "${SQ9_C09_Q05_CASE}" 1025 1.0d9 0.9d0 0.5d0 S9C09Q05 node07 fixed_sq_stability_control
make_fixed_sq_scratch "${SQ9_C09_Q1_CASE}" 1025 1.0d9 0.9d0 1.0d0 S9C09Q1 node07 fixed_sq_stability_control

previous=''
for case_dir in "${SQ08_CASE}" "${SQ09_CASE}" "${SQ095_CASE}" "${SQ099_RERUN_CASE}"; do
    if [[ -z "${previous}" ]]; then
        jid=$(qsub "${case_dir}/pbs/run.pbs")
        dependency=none
    else
        jid=$(qsub -W depend=afterany:${previous} "${case_dir}/pbs/run.pbs")
        dependency=${previous}
    fi
    printf 'node05 case=%s job=%s depends_on=%s\n' "${case_dir}" "${jid}" "${dependency}"
    previous=${jid}
done

previous=''
for case_dir in "${SQ9_C05_Q05_CASE}" "${SQ9_C05_Q1_CASE}" \
    "${SQ9_C09_Q05_CASE}" "${SQ9_C09_Q1_CASE}"; do
    if [[ -z "${previous}" ]]; then
        jid=$(qsub "${case_dir}/pbs/run.pbs")
        dependency=none
    else
        jid=$(qsub -W depend=afterany:${previous} "${case_dir}/pbs/run.pbs")
        dependency=${previous}
    fi
    printf 'node07 case=%s job=%s depends_on=%s\n' "${case_dir}" "${jid}" "${dependency}"
    previous=${jid}
done
