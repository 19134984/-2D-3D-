#!/usr/bin/env bash
set -euo pipefail

# Ra=1e8, chi_nu=0.4 固定 Sq 长时间统计：
# Sq=0.05/0.08 从有效的 300 t_ff 重启，Sq=0.10 从零计算，目标均为 1000 t_ff。

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST/SQ_SCAN/Ra1e8/chinu0.4
STAGED_SOURCE="${ROOT}/_staging/first_solver_restart_tolerance_fixed.F90"
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
SEED_005="${ROOT}/sq0.0500000_v2/results"
SEED_008="${ROOT}/sq0.0800000/results"

if [[ ! -f "${STAGED_SOURCE}" ]]; then
    printf 'Missing staged source: %s\n' "${STAGED_SOURCE}" >&2
    exit 2
fi
for seed_dir in "${SEED_005}" "${SEED_008}"; do
    if [[ ! -f "${seed_dir}/reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta" ]]; then
        printf 'Missing restart metadata: %s\n' "${seed_dir}" >&2
        exit 2
    fi
done

configure_source() {
    local source_file=$1
    local sq_value=$2
    local load_init=$3

    cp "${STAGED_SOURCE}" "${source_file}"
    sed -i 's/\r$//' "${source_file}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e "s/real(kind=8), parameter :: Sq=1.0d0 ! 固定奇模态松弛率对照：tau_q=1/real(kind=8), parameter :: Sq=${sq_value} ! 固定奇模态松弛率长时间测试/" \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=513, ny=513/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d8/' \
        -e 's/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=0.4d0/' \
        -e "s/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=${load_init}/" \
        "${source_file}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${source_file}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${source_file}"
    grep -Fq "parameter :: Sq=${sq_value}" "${source_file}"
    grep -Fq 'parameter :: nx=513, ny=513' "${source_file}"
    grep -Fq 'parameter :: Rayleigh=1.0d8' "${source_file}"
    grep -Fq 'parameter :: chi_nu=0.4d0' "${source_file}"
    grep -Fq "parameter :: loadInitField=${load_init}" "${source_file}"
    grep -Fq 'parameter :: unsteadyRunDuration=1000.0d0' "${source_file}"
    grep -Fq '0.5d0/timeUnit+64.0d0*epsilon(1.0d0)' "${source_file}"
}

write_metadata() {
    local case_dir=$1
    local case_name=$2
    local sq_value=$3
    local load_init=$4
    local seed_dir=$5
    local source_file="${case_dir}/source/${SOURCE_NAME}"

    sha256sum "${STAGED_SOURCE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${source_file}" > "${case_dir}/source/case_source.sha256"
    {
        printf 'case=%s\nnode=node07\n' "${case_name}"
        printf 'Ra=1.0d8\nPr=0.7d0\nnx=513\nny=513\nchi_nu=0.4d0\nchi_b=0.0d0\n'
        printf 'flow_policy=FLOW_ODD_FIXED_SQ\nSq=%s\n' "${sq_value}"
        printf 'unsteadyRunDuration=1000.0d0\nloadInitField=%s\n' "${load_init}"
        printf 'nonfinite_interval_tff=1.0d0\nflow_odd_interval_tff=0.1d0\n'
        printf 'restart_seed=%s\n' "${seed_dir}"
    } > "${case_dir}/case_settings.txt"
}

write_pbs() {
    local case_dir=$1
    local job_name=$2
    local source_file="${case_dir}/source/${SOURCE_NAME}"
    local pbs_file="${case_dir}/pbs/run.pbs"

    cat > "${pbs_file}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=node07:ppn=1
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
NVFLAGS=(-cpp -O3 -acc -gpu=cc60 -Minfo=accel)
printf '%q ' "\${NVFLAGS[@]}" > build_flags.txt
printf '\n' >> build_flags.txt
"\${NVFORTRAN}" --version > compiler_version.txt 2>&1
"\${NVFORTRAN}" "\${NVFLAGS[@]}" "\${SRC}" -o solver.exe > compile.log 2>&1
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
    Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat DissipationHistory_2DOpenaccLBMCDE.dat \
    PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q5.dat \
    FlowOddMomentHistory_2DOpenaccLBMCDE_D2Q5.dat 2>/dev/null; then
    printf 'non-finite value found in solver output/history\n' >> solver.stdout
    run_status=86
fi
printf '%s\n' "\${run_status}" > run.status
date '+CASE_END %F %T %z' >> timing.txt
exit "\${run_status}"
EOF
    chmod 750 "${pbs_file}"
}

make_restart_case() {
    local case_name=$1
    local job_name=$2
    local sq_value=$3
    local seed_dir=$4
    local keep_seed_diagnostics=$5
    local case_dir="${ROOT}/${case_name}"

    if [[ -e "${case_dir}" ]]; then
        printf 'Refusing to overwrite existing case: %s\n' "${case_dir}" >&2
        exit 3
    fi
    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    cp -a "${seed_dir}/." "${case_dir}/results/"
    cp "${seed_dir}/run.status" "${case_dir}/results/seed_300tff_run.status"
    cp "${seed_dir}/timing.txt" "${case_dir}/results/seed_300tff_timing.txt"
    cp "${seed_dir}/runtime_source.sha256" "${case_dir}/results/seed_300tff_runtime_source.sha256"
    if [[ "${keep_seed_diagnostics}" -eq 0 ]]; then
        # 0.08 的旧诊断由 nvfortran 折成多行；原目录保留，新目录从 300 t_ff 建立单行诊断分段。
        rm -f "${case_dir}/results/PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q5.dat"
        rm -f "${case_dir}/results/FlowOddMomentHistory_2DOpenaccLBMCDE_D2Q5.dat"
    fi
    configure_source "${case_dir}/source/${SOURCE_NAME}" "${sq_value}" 1
    write_metadata "${case_dir}" "${case_name}" "${sq_value}" 1 "${seed_dir}"
    write_pbs "${case_dir}" "${job_name}"
}

make_fresh_010() {
    local case_name=sq0.1000000_t1000
    local case_dir="${ROOT}/${case_name}"

    if [[ -e "${case_dir}" ]]; then
        printf 'Refusing to overwrite existing case: %s\n' "${case_dir}" >&2
        exit 3
    fi
    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    configure_source "${case_dir}/source/${SOURCE_NAME}" 0.10d0 0
    write_metadata "${case_dir}" "${case_name}" 0.10d0 0 none
    write_pbs "${case_dir}" SQ8C10T1000
}

make_restart_case sq0.0500000_t1000 SQ8C05T1000 0.05d0 "${SEED_005}" 1
make_restart_case sq0.0800000_t1000_v2 SQ8C08T1V2 0.08d0 "${SEED_008}" 0
make_fresh_010

sq005_job=$(qsub "${ROOT}/sq0.0500000_t1000/pbs/run.pbs")
sq008_job=$(qsub -W depend=afterany:${sq005_job} "${ROOT}/sq0.0800000_t1000_v2/pbs/run.pbs")
sq010_job=$(qsub -W depend=afterany:${sq008_job} "${ROOT}/sq0.1000000_t1000/pbs/run.pbs")

{
    printf 'sq0.0500000_t1000=%s\n' "${sq005_job}"
    printf 'sq0.0800000_t1000_v2=%s\n' "${sq008_job}"
    printf 'sq0.1000000_t1000=%s\n' "${sq010_job}"
} | tee "${ROOT}/long_run_submission_manifest.txt"
