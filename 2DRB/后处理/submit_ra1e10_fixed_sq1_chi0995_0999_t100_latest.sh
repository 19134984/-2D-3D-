#!/usr/bin/env bash
set -euo pipefail

# Ra=1e10、固定 Sq=1 的近极限 chi_nu 短筛选。
# 两个算例均使用当前本地最新版源代码，从静止初场运行到 100 t_ff；
# 若当前算例检测到非数则退出，afterany 任务链仍继续下一个参数点。

FLOW_ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
ROOT=${FLOW_ROOT}/FIXED_SQ1/Ra1e10
STAGE=${FLOW_ROOT}/_staging/20260901_ra1e10_fixed_sq1_chi0995_0999_t100/latest_local.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90

CASE_0995=${ROOT}/chinu0.995_t100_latest
CASE_0999=${ROOT}/chinu0.999_t100_latest

test -f "${STAGE}"
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
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=2049, ny=2049/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi_value}/" \
        -e 's/integer(kind=4), parameter :: loadInitField=1/integer(kind=4), parameter :: loadInitField=0/' \
        -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=100.0d0/' \
        "${src}"
    sed -i "/^#ifdef FLOW_ODD_FIXED_SQ$/,/^#endif$/ s|^[[:space:]]*real(kind=8), parameter :: Sq=.*$|        real(kind=8), parameter :: Sq=1.0d0 ! 固定奇矩松弛率对照：tau_q=1|" "${src}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq 'parameter :: nx=2049, ny=2049' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d10' "${src}"
    grep -Fq "parameter :: chi_nu=${chi_value}" "${src}"
    grep -Fq 'parameter :: Sq=1.0d0' "${src}"
    grep -Fq 'parameter :: loadInitField=0' "${src}"
    grep -Fq 'parameter :: unsteadyRunDuration=100.0d0' "${src}"
    grep -Fq 'parameter :: nonFiniteCheckIntervalTf=1.0d0' "${src}"
    grep -Fq 'parameter :: flowOddMomentDiagnosticIntervalTf=0.1d0' "${src}"
}

write_pbs() {
    local case_dir=$1
    local job_name=$2
    local src=${case_dir}/source/${SOURCE_NAME}
    local result=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs

    cat > "${pbs}" <<'PBSEOF'
#!/usr/bin/env bash
#PBS -N __JOB__
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=96:00:00
#PBS -o __RESULT__/pbs.stdout
#PBS -e __RESULT__/pbs.stderr

set -u
SRC=__SOURCE__
RESULT_DIR=__RESULT__
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=${NVHPC_ROOT}/compilers/bin/nvfortran

mkdir -p "${RESULT_DIR}"
cd "${RESULT_DIR}" || exit 90
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
        -e "s|__RESULT__|${result}|g" \
        "${pbs}"
    chmod 700 "${pbs}"
}

make_case() {
    local case_dir=$1
    local chi_value=$2
    local job_name=$3
    local src=${case_dir}/source/${SOURCE_NAME}

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    configure_source "${src}" "${chi_value}"
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${src}" > "${case_dir}/source/case_source.sha256"
    write_pbs "${case_dir}" "${job_name}"

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
loadInitField=0
target_time_tf=100.0d0
nonfinite_population_interval_tff=1.0d0
flow_odd_interval_tff=0.1d0
node=node07
purpose=fixed_sq1_near_unity_chinu_short_stability_screen
EOF
}

make_case "${CASE_0995}" 0.995d0 S10C0995Q1
make_case "${CASE_0999}" 0.999d0 S10C0999Q1

jid_0995=$(qsub "${CASE_0995}/pbs/run.pbs")
jid_0999=$(qsub -W depend=afterany:${jid_0995} "${CASE_0999}/pbs/run.pbs")
printf 'node07 chi=0.995 Sq=1 job=%s\n' "${jid_0995}"
printf 'node07 chi=0.999 Sq=1 job=%s depends_on=%s\n' "${jid_0999}" "${jid_0995}"
