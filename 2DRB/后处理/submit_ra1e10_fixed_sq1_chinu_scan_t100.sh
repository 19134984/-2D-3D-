#!/usr/bin/env bash
set -euo pipefail

# Ra=1e10、固定 Sq=1 的 chi_nu 短扫描。
# chi_nu=0、0.1--0.9、0.95、0.99；全部从静止初场计算到 100 t_ff。
# 任一算例检测到非数后退出，afterany 任务链继续下一个参数点。

FLOW_ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
ROOT=${FLOW_ROOT}/FIXED_SQ1/Ra1e10
STAGE=${FLOW_ROOT}/_staging/20260829_ra1e10_fixed_sq1_chi_t100/latest_local.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
NODE07_DEPEND_AFTER=6413.master
NODE05_DEPEND_AFTER=6414.master

CHIS=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99)
CHI_FORTRAN=(0.0d0 0.1d0 0.2d0 0.3d0 0.4d0 0.5d0 0.6d0 0.7d0 0.8d0 0.9d0 0.95d0 0.99d0)
JOB_TAGS=(00 01 02 03 04 05 06 07 08 09 095 099)

test -f "${STAGE}"
for chi in "${CHIS[@]}"; do
    target=${ROOT}/chinu${chi}_t100_latest
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_source() {
    local src=$1
    local chi=$2

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=2049, ny=2049/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi}/" \
        -e 's/integer(kind=4), parameter :: loadInitField=1/integer(kind=4), parameter :: loadInitField=0/' \
        -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=100.0d0/' \
        "${src}"
    sed -i "/^#ifdef FLOW_ODD_FIXED_SQ$/,/^#endif$/ s|^[[:space:]]*real(kind=8), parameter :: Sq=.*$|        real(kind=8), parameter :: Sq=1.0d0 ! 固定奇矩松弛率对照|" "${src}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq 'parameter :: nx=2049, ny=2049' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d10' "${src}"
    grep -Fq "parameter :: chi_nu=${chi}" "${src}"
    grep -Fq 'parameter :: Sq=1.0d0' "${src}"
    grep -Fq 'parameter :: loadInitField=0' "${src}"
    grep -Fq 'parameter :: unsteadyRunDuration=100.0d0' "${src}"
    grep -Fq 'parameter :: nonFiniteCheckIntervalTf=1.0d0' "${src}"
    grep -Fq 'parameter :: flowOddMomentDiagnosticIntervalTf=0.1d0' "${src}"
}

write_pbs() {
    local case_dir=$1
    local job_name=$2
    local node_name=$3
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
        -e "s|__NODE__|${node_name}|g" \
        -e "s|__SOURCE__|${src}|g" \
        -e "s|__RESULT__|${result}|g" \
        "${pbs}"
    chmod 700 "${pbs}"
}

make_case() {
    local chi_label=$1
    local chi_value=$2
    local job_tag=$3
    local node_name=$4
    local case_dir=${ROOT}/chinu${chi_label}_t100_latest
    local job_name=S10C${job_tag}T100

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    configure_source "${case_dir}/source/${SOURCE_NAME}" "${chi_value}"
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${case_dir}/source/${SOURCE_NAME}" > "${case_dir}/source/case_source.sha256"
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
node=${node_name}
purpose=fixed_sq1_chinu_short_stability_screen
EOF
    write_pbs "${case_dir}" "${job_name}" "${node_name}"
    printf '%s\n' "${case_dir}"
}

node07_cases=()
node05_cases=()
for index in "${!CHIS[@]}"; do
    if [[ "${index}" -le 5 ]]; then
        node07_cases+=("$(make_case "${CHIS[index]}" "${CHI_FORTRAN[index]}" "${JOB_TAGS[index]}" node07)")
    else
        node05_cases+=("$(make_case "${CHIS[index]}" "${CHI_FORTRAN[index]}" "${JOB_TAGS[index]}" node05)")
    fi
done

previous=${NODE07_DEPEND_AFTER}
for index in "${!node07_cases[@]}"; do
    jid=$(qsub -W depend=afterany:${previous} "${node07_cases[index]}/pbs/run.pbs")
    printf 'node07 chi=%s job=%s depends_on=%s\n' "${CHIS[index]}" "${jid}" "${previous}"
    previous=${jid}
done

previous=${NODE05_DEPEND_AFTER}
for local_index in "${!node05_cases[@]}"; do
    global_index=$((local_index+6))
    jid=$(qsub -W depend=afterany:${previous} "${node05_cases[local_index]}/pbs/run.pbs")
    printf 'node05 chi=%s job=%s depends_on=%s\n' "${CHIS[global_index]}" "${jid}" "${previous}"
    previous=${jid}
done
