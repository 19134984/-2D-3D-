#!/usr/bin/env bash
set -euo pipefail

# Ra=1e9 第一阶段短筛选：固定 chi_nu=0.4，仅改变流场奇矩松弛率 Sq。
# 两个算例均从静止初场运行至 100 t_ff；检测到非数时由求解器提前结束。
# 调用前必须把最新本地源文件上传到 STAGE；已有算例目录不会被覆盖。

STAGE=/data2/XLLi/LBMCDE/FLOW-TEST/_staging/20260827_ra1e9_sq_stage1/latest_local.F90
SQROOT=/data2/XLLi/LBMCDE/FLOW-TEST/SQ_SCAN/Ra1e9/chinu0.4
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
DEPEND_AFTER=6313.master

C50=${SQROOT}/sq0.5000000_t100
C10=${SQROOT}/sq1.0000000_t100

test -f "${STAGE}"

for target in "${C50}" "${C10}"; do
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_source() {
    local src=$1
    local sq=$2

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        -e "s|^[[:space:]]*real(kind=8), parameter :: Sq=1.0d0 ! 固定奇模态松弛率对照：tau_q=1$|        real(kind=8), parameter :: Sq=${sq} ! 固定奇矩松弛率诊断|" \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=1025, ny=1025/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d9/' \
        -e 's/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=0.4d0/' \
        -e 's/integer(kind=4), parameter :: loadInitField=1/integer(kind=4), parameter :: loadInitField=0/' \
        -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=100.0d0/' \
        "${src}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq "parameter :: Sq=${sq}" "${src}"
    grep -Fq 'parameter :: nx=1025, ny=1025' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d9' "${src}"
    grep -Fq 'parameter :: chi_nu=0.4d0' "${src}"
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
#PBS -l walltime=12:00:00
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

if [[ ! -x "${NVFORTRAN}" ]]; then
    printf '127\n' > compile.status
    printf '125\n' > run.status
    exit 127
fi

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
    local sq=$2
    local job=$3

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    configure_source "${case_dir}/source/${SOURCE_NAME}" "${sq}"
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${case_dir}/source/${SOURCE_NAME}" > "${case_dir}/source/case_source.sha256"
    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job}
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=1.0d9
Pr=0.7d0
nx=1025
ny=1025
chi_nu=0.4d0
flow_policy=FLOW_ODD_FIXED_SQ
Sq=${sq}
loadInitField=0
target_time_tf=100.0
nonfinite_macro_interval_tff=1.0
population_and_odd_moment_interval_tff=0.1
purpose=short_stability_screen_not_statistical_convergence
EOF
    write_pbs "${case_dir}" "${job}"
}

make_case "${C50}" 0.5d0 S9C4Q50
make_case "${C10}" 1.0d0 S9C4Q10

jid50=$(qsub -W depend=afterany:${DEPEND_AFTER} "${C50}/pbs/run.pbs")
jid10=$(qsub -W depend=afterany:${jid50} "${C10}/pbs/run.pbs")

printf 'sq050=%s\nsq100=%s\n' "${jid50}" "${jid10}"
