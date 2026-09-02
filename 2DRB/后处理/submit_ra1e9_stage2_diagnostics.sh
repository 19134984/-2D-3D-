#!/usr/bin/env bash
set -euo pipefail

# Ra=1e9 第二阶段诊断：
# 1) 原始魔法参数 chi_nu=0.4，以 0.1 t_ff 监测 f_i/g_i，运行上限 60 t_ff；
# 2) chi_nu=0.4 固定 Sq=0.7/0.8/0.9，各从静止初场运行 100 t_ff；
# 3) 原始 chi_nu=0.1/0.95/0.99 和固定 Sq=1 从 t=100 续算到 t=300。
# 所有目标目录均为新目录；任一任务非数退出后，afterany 依赖继续后续任务。

STAGE=/data2/XLLi/LBMCDE/FLOW-TEST/_staging/20260828_ra1e9_stage2/latest_local.F90
FLOW_ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
ORIGINAL_ROOT=${FLOW_ROOT}/ORIGINAL/Ra1e9
SQ_ROOT=${FLOW_ROOT}/SQ_SCAN/Ra1e9/chinu0.4
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90

DIAG_CASE=${ORIGINAL_ROOT}/chinu0.4_diag_population0p1_t60_latest
SQ07_CASE=${SQ_ROOT}/sq0.7000000_t100
SQ08_CASE=${SQ_ROOT}/sq0.8000000_t100
SQ09_CASE=${SQ_ROOT}/sq0.9000000_t100
O01_CASE=${ORIGINAL_ROOT}/chinu0.1_t300_latest
O095_CASE=${ORIGINAL_ROOT}/chinu0.95_t300_latest
O099_CASE=${ORIGINAL_ROOT}/chinu0.99_t300_latest
SQ1_CASE=${SQ_ROOT}/sq1.0000000_t300

O01_SEED=${ORIGINAL_ROOT}/chinu0.1_t100_latest/results
O095_SEED=${ORIGINAL_ROOT}/chinu0.95_t100_latest/results
O099_SEED=${ORIGINAL_ROOT}/chinu0.99_t100_latest/results
SQ1_SEED=${SQ_ROOT}/sq1.0000000_t100/results

test -f "${STAGE}"
for target in "${DIAG_CASE}" "${SQ07_CASE}" "${SQ08_CASE}" "${SQ09_CASE}" \
              "${O01_CASE}" "${O095_CASE}" "${O099_CASE}" "${SQ1_CASE}"; do
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_common() {
    local src=$1
    local chi=$2
    local load=$3
    local duration=$4
    local nonfinite_interval=$5

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=1025, ny=1025/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d9/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi}/" \
        -e "s/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=${load}/" \
        -e "s/real(kind=8), parameter :: nonFiniteCheckIntervalTf=1.0d0/real(kind=8), parameter :: nonFiniteCheckIntervalTf=${nonfinite_interval}/" \
        -e "s/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=${duration}/" \
        "${src}"

    grep -Fq 'parameter :: nx=1025, ny=1025' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d9' "${src}"
    grep -Fq "parameter :: chi_nu=${chi}" "${src}"
    grep -Fq "parameter :: loadInitField=${load}" "${src}"
    grep -Fq "parameter :: nonFiniteCheckIntervalTf=${nonfinite_interval}" "${src}"
    grep -Fq "parameter :: unsteadyRunDuration=${duration}" "${src}"
    grep -Fq 'parameter :: flowOddMomentDiagnosticIntervalTf=0.1d0' "${src}"
}

configure_original() {
    local src=$1
    local chi=$2
    local load=$3
    local duration=$4
    local nonfinite_interval=$5

    configure_common "${src}" "${chi}" "${load}" "${duration}" "${nonfinite_interval}"
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
    local chi=$2
    local sq=$3
    local load=$4
    local duration=$5
    local nonfinite_interval=$6

    configure_common "${src}" "${chi}" "${load}" "${duration}" "${nonfinite_interval}"
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
    local job_name=$3
    local src=${case_dir}/source/${SOURCE_NAME}
    local result=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs

    cat > "${pbs}" <<'PBSEOF'
#!/usr/bin/env bash
#PBS -N __JOB__
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=48:00:00
#PBS -o __RESULT__/pbs.stdout
#PBS -e __RESULT__/pbs.stderr

set -u
SRC=__SOURCE__
SEED_DIR=__SEED__
RESULT_DIR=__RESULT__
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=${NVHPC_ROOT}/compilers/bin/nvfortran

mkdir -p "${RESULT_DIR}"
cd "${RESULT_DIR}" || exit 90

if [[ -n "${SEED_DIR}" ]]; then
    if [[ ! -f "${SEED_DIR}/run.status" ]] || [[ "$(<"${SEED_DIR}/run.status")" != 0 ]]; then
        printf 'Invalid restart seed run.status: %s\n' "${SEED_DIR}" > restart_preflight.error
        printf '121\n' > run.status
        exit 121
    fi
    if ! awk '$1=="time_tf" {t=$2} END {exit !(t>99.9 && t<100.1)}' \
        "${SEED_DIR}/reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta"; then
        printf 'Restart seed is not at 100 t_ff: %s\n' "${SEED_DIR}" > restart_preflight.error
        printf '122\n' > run.status
        exit 122
    fi
    if ! grep -Fq 'cumulativeStatisticSampleCount 101' \
        "${SEED_DIR}/reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta"; then
        printf 'Restart seed sample count is not 101: %s\n' "${SEED_DIR}" > restart_preflight.error
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
        reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta
    )
    for required in "${required_histories[@]}"; do
        if [[ ! -f "${SEED_DIR}/${required}" ]]; then
            printf 'Missing restart file: %s\n' "${SEED_DIR}/${required}" > restart_preflight.error
            printf '124\n' > run.status
            exit 124
        fi
        cp -a "${SEED_DIR}/${required}" "${RESULT_DIR}/"
    done
    mapfile -t restart_bins < <(find "${SEED_DIR}" -maxdepth 1 -type f \
        -name 'reloadFile2DOpenaccLBMCDE_D2Q5-*.bin' -print)
    if [[ "${#restart_bins[@]}" -lt 1 ]]; then
        printf 'Missing restart binary: %s\n' "${SEED_DIR}" > restart_preflight.error
        printf '124\n' > run.status
        exit 124
    fi
    cp -a "${restart_bins[@]}" "${RESULT_DIR}/"
    cp "${SEED_DIR}/run.status" seed_100tff_run.status
    cp "${SEED_DIR}/runtime_source.sha256" seed_100tff_runtime_source.sha256
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
        -e "s|__SOURCE__|${src}|g" \
        -e "s|__SEED__|${seed_dir}|g" \
        -e "s|__RESULT__|${result}|g" \
        "${pbs}"
    chmod 700 "${pbs}"
}

create_case_dir() {
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
    local load=$3
    local duration=$4
    local nonfinite_interval=$5
    local seed_dir=$6
    local job_name=$7
    local purpose=$8

    create_case_dir "${case_dir}"
    configure_original "${case_dir}/source/${SOURCE_NAME}" "${chi}" "${load}" \
        "${duration}" "${nonfinite_interval}"
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
loadInitField=${load}
target_time_tf=${duration}
restart_seed=${seed_dir:-none}
nonfinite_population_interval_tff=${nonfinite_interval}
flow_odd_interval_tff=0.1d0
purpose=${purpose}
EOF
    write_pbs "${case_dir}" "${seed_dir}" "${job_name}"
}

make_fixed_case() {
    local case_dir=$1
    local sq=$2
    local load=$3
    local duration=$4
    local nonfinite_interval=$5
    local seed_dir=$6
    local job_name=$7
    local purpose=$8

    create_case_dir "${case_dir}"
    configure_fixed_sq "${case_dir}/source/${SOURCE_NAME}" 0.4d0 "${sq}" "${load}" \
        "${duration}" "${nonfinite_interval}"
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
loadInitField=${load}
target_time_tf=${duration}
restart_seed=${seed_dir:-none}
nonfinite_population_interval_tff=${nonfinite_interval}
flow_odd_interval_tff=0.1d0
purpose=${purpose}
EOF
    write_pbs "${case_dir}" "${seed_dir}" "${job_name}"
}

make_original_case "${DIAG_CASE}" 0.4d0 0 60.0d0 0.1d0 '' O9C04D01 population_first_nonfinite_diagnostic
make_fixed_case "${SQ07_CASE}" 0.7d0 0 100.0d0 1.0d0 '' S9C4Q70 fixed_sq_stability_screen
make_fixed_case "${SQ08_CASE}" 0.8d0 0 100.0d0 1.0d0 '' S9C4Q80 fixed_sq_stability_screen
make_fixed_case "${SQ09_CASE}" 0.9d0 0 100.0d0 1.0d0 '' S9C4Q90 fixed_sq_stability_screen
make_original_case "${O01_CASE}" 0.1d0 1 300.0d0 1.0d0 "${O01_SEED}" O9C01R300 delayed_instability_continuation
make_original_case "${O095_CASE}" 0.95d0 1 300.0d0 1.0d0 "${O095_SEED}" O9C095R300 delayed_instability_continuation
make_original_case "${O099_CASE}" 0.99d0 1 300.0d0 1.0d0 "${O099_SEED}" O9C099R300 delayed_instability_continuation
make_fixed_case "${SQ1_CASE}" 1.0d0 1 300.0d0 1.0d0 "${SQ1_SEED}" S9C4Q1R300 delayed_instability_continuation

cases=("${DIAG_CASE}" "${SQ07_CASE}" "${SQ08_CASE}" "${SQ09_CASE}" \
       "${O01_CASE}" "${O095_CASE}" "${O099_CASE}" "${SQ1_CASE}")
previous=''
for case_dir in "${cases[@]}"; do
    if [[ -z "${previous}" ]]; then
        jid=$(qsub "${case_dir}/pbs/run.pbs")
    else
        jid=$(qsub -W depend=afterany:${previous} "${case_dir}/pbs/run.pbs")
    fi
    printf '%s=%s\n' "${case_dir#${FLOW_ROOT}/}" "${jid}"
    previous=${jid}
done
