#!/usr/bin/env bash
set -euo pipefail

# Ra=1e8, chi_nu=0.4 固定 Sq 机制诊断：
# 先执行 5 t_ff GPU 试跑，再将四个 300 t_ff 算例分配到 node05/node07。
# 每个算例保存独立源代码快照、源文件哈希、PBS 文件和运行结果。

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST/SQ_SCAN/Ra1e8/chinu0.4
STAGED_SOURCE="${ROOT}/_staging/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90"
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
NODE05_TAIL=6254.master
NODE07_TAIL=6256.master

if [[ ! -f "${STAGED_SOURCE}" ]]; then
    printf 'Missing staged source: %s\n' "${STAGED_SOURCE}" >&2
    exit 2
fi

mkdir -p "${ROOT}"

make_case() {
    local case_name=$1
    local job_name=$2
    local node_name=$3
    local sq_value=$4
    local run_duration=$5
    local case_dir="${ROOT}/${case_name}"
    local source_file="${case_dir}/source/${SOURCE_NAME}"
    local pbs_file="${case_dir}/pbs/run.pbs"

    if [[ -e "${case_dir}" ]]; then
        printf 'Refusing to overwrite existing case: %s\n' "${case_dir}" >&2
        exit 3
    fi

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    cp "${STAGED_SOURCE}" "${source_file}"
    sed -i 's/\r$//' "${source_file}"

    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e "s/real(kind=8), parameter :: Sq=1.0d0 ! 固定奇模态松弛率对照：tau_q=1/real(kind=8), parameter :: Sq=${sq_value} ! 固定奇模态松弛率机制诊断/" \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=513, ny=513/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d8/' \
        -e 's/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=0.4d0/' \
        -e "s/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=${run_duration}/" \
        "${source_file}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${source_file}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${source_file}"
    grep -Fq "parameter :: Sq=${sq_value}" "${source_file}"
    grep -Fq 'parameter :: nx=513, ny=513' "${source_file}"
    grep -Fq 'parameter :: Rayleigh=1.0d8' "${source_file}"
    grep -Fq 'parameter :: chi_nu=0.4d0' "${source_file}"
    grep -Fq "parameter :: unsteadyRunDuration=${run_duration}" "${source_file}"
    grep -Fq 'flowOddMomentDiagnosticIntervalTf=0.1d0' "${source_file}"
    grep -Fq 'FlowOddMomentHistory_2DOpenaccLBMCDE_D2Q5.dat' "${source_file}"

    sha256sum "${STAGED_SOURCE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${source_file}" > "${case_dir}/source/case_source.sha256"
    {
        printf 'case=%s\n' "${case_name}"
        printf 'node=%s\n' "${node_name}"
        printf 'Ra=1.0d8\nPr=0.7d0\nnx=513\nny=513\n'
        printf 'chi_nu=0.4d0\nchi_b=0.0d0\n'
        printf 'flow_policy=FLOW_ODD_FIXED_SQ\nSq=%s\n' "${sq_value}"
        printf 'unsteadyRunDuration=%s\n' "${run_duration}"
        printf 'odd_moment_diagnostic_interval_tff=0.1\n'
    } > "${case_dir}/case_settings.txt"

    cat > "${pbs_file}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=${node_name}:ppn=1
#PBS -l walltime=96:00:00
#PBS -o ${case_dir}/results/pbs.stdout
#PBS -e ${case_dir}/results/pbs.stderr

set -u
SRC=${source_file}
RESULT_DIR=${case_dir}/results
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=\${NVHPC_ROOT}/compilers/bin/nvfortran

mkdir -p "\${RESULT_DIR}"
cd "\${RESULT_DIR}" || exit 90
export PATH="\${NVHPC_ROOT}/compilers/bin:\${PATH}"
export LD_LIBRARY_PATH="\${NVHPC_ROOT}/compilers/lib:\${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1
export ACC_DEVICE_TYPE=nvidia

printf 'PBS_JOBID=%s\n' "\${PBS_JOBID:-unknown}" > job_identity.txt
printf 'HOSTNAME=%s\n' "\$(hostname)" >> job_identity.txt
date '+CASE_START %F %T %z' > timing.txt
sha256sum "\${SRC}" > runtime_source.sha256
printf '%s\n' "\${NVFORTRAN}" > compiler_path.txt
"\${NVFORTRAN}" --version > compiler_version.txt 2>&1
nvidia-smi > nvidia_smi_before.log 2>&1 || true

if [[ ! -x "\${NVFORTRAN}" ]]; then
    printf '127\n' > compile.status
    printf '125\n' > run.status
    exit 127
fi

BUILD_FLAGS=(-cpp -O3 -acc -gpu=cc60 -Minfo=accel)
printf '%q ' "\${BUILD_FLAGS[@]}" > build_flags.txt
printf '\n' >> build_flags.txt
"\${NVFORTRAN}" "\${BUILD_FLAGS[@]}" "\${SRC}" -o solver.exe > compile.log 2>&1
compile_status=\$?
printf '%s\n' "\${compile_status}" > compile.status
if [[ "\${compile_status}" -ne 0 ]]; then
    printf '125\n' > run.status
    exit "\${compile_status}"
fi

set +e
./solver.exe > solver.stdout 2>&1
run_status=\$?
set -e
if [[ "\${run_status}" -eq 0 ]] && grep -Eqi '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)' \
    solver.stdout Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat \
    Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat \
    DissipationHistory_2DOpenaccLBMCDE.dat \
    FlowOddMomentHistory_2DOpenaccLBMCDE_D2Q5.dat 2>/dev/null; then
    printf 'non-finite value found in solver output/history\n' >> solver.stdout
    run_status=86
fi
printf '%s\n' "\${run_status}" > run.status
nvidia-smi > nvidia_smi_after.log 2>&1 || true
date '+CASE_END %F %T %z' >> timing.txt
exit "\${run_status}"
EOF
    chmod 750 "${pbs_file}"
}

make_case gpu_smoke_sq0.0639672 SQ8SMOKE node05 0.0639671806308998d0 5.0d0
make_case sq0.0388777 SQ8C0389 node05 0.0388776875914752d0 300.0d0
make_case sq0.0500000 SQ8C0500 node05 0.05d0 300.0d0
make_case sq0.0639672 SQ8C0640 node07 0.0639671806308998d0 300.0d0
make_case sq0.0800000 SQ8C0800 node07 0.08d0 300.0d0

smoke_job=$(qsub -W depend=afterany:${NODE05_TAIL} "${ROOT}/gpu_smoke_sq0.0639672/pbs/run.pbs")
sq0389_job=$(qsub -W depend=afterok:${smoke_job} "${ROOT}/sq0.0388777/pbs/run.pbs")
sq0500_job=$(qsub -W depend=afterok:${sq0389_job} "${ROOT}/sq0.0500000/pbs/run.pbs")
sq0640_job=$(qsub -W depend=afterany:${NODE07_TAIL},afterok:${smoke_job} "${ROOT}/sq0.0639672/pbs/run.pbs")
sq0800_job=$(qsub -W depend=afterok:${sq0640_job} "${ROOT}/sq0.0800000/pbs/run.pbs")

{
    printf 'gpu_smoke_sq0.0639672=%s\n' "${smoke_job}"
    printf 'sq0.0388777=%s\n' "${sq0389_job}"
    printf 'sq0.0500000=%s\n' "${sq0500_job}"
    printf 'sq0.0639672=%s\n' "${sq0640_job}"
    printf 'sq0.0800000=%s\n' "${sq0800_job}"
} | tee "${ROOT}/submission_manifest.txt"
