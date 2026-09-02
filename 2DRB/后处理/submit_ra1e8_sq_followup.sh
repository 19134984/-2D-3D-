#!/usr/bin/env bash
set -euo pipefail

# 第一套 D2Q9-TRT/D2Q5-TRT 在 Ra=1e8、chi_nu=0.4 下的后续固定 Sq 测试。
# node07 顺序：0.0388777 的 60 t_ff 加密诊断 -> 0.05 的 300 t_ff -> 0.08 从 300 续算至 1000 t_ff。
# 使用 afterany 串联，确保预期会发散的诊断任务不会阻断后续独立算例。

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST/SQ_SCAN/Ra1e8/chinu0.4
STAGED_SOURCE="${ROOT}/_staging/first_solver_diag_coord_fixed.F90"
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
SEED_008="${ROOT}/sq0.0800000/results"

if [[ ! -f "${STAGED_SOURCE}" ]]; then
    printf 'Missing staged source: %s\n' "${STAGED_SOURCE}" >&2
    exit 2
fi
if [[ ! -f "${SEED_008}/reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta" ]]; then
    printf 'Missing Sq=0.08 restart metadata: %s\n' "${SEED_008}" >&2
    exit 2
fi

configure_source() {
    local source_file=$1
    local sq_value=$2
    local run_duration=$3
    local load_init=$4
    local nonfinite_interval=$5

    cp "${STAGED_SOURCE}" "${source_file}"
    sed -i 's/\r$//' "${source_file}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e "s/real(kind=8), parameter :: Sq=1.0d0 ! 固定奇模态松弛率对照：tau_q=1/real(kind=8), parameter :: Sq=${sq_value} ! 固定奇模态松弛率机制诊断/" \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=513, ny=513/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d8/' \
        -e 's/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=0.4d0/' \
        -e "s/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=${load_init}/" \
        -e "s/real(kind=8), parameter :: nonFiniteCheckIntervalTf=1.0d0/real(kind=8), parameter :: nonFiniteCheckIntervalTf=${nonfinite_interval}/" \
        -e "s/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=${run_duration}/" \
        "${source_file}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${source_file}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${source_file}"
    grep -Fq "parameter :: Sq=${sq_value}" "${source_file}"
    grep -Fq 'parameter :: nx=513, ny=513' "${source_file}"
    grep -Fq 'parameter :: Rayleigh=1.0d8' "${source_file}"
    grep -Fq 'parameter :: chi_nu=0.4d0' "${source_file}"
    grep -Fq "parameter :: loadInitField=${load_init}" "${source_file}"
    grep -Fq "parameter :: nonFiniteCheckIntervalTf=${nonfinite_interval}" "${source_file}"
    grep -Fq "parameter :: unsteadyRunDuration=${run_duration}" "${source_file}"
    grep -Fq 'xQOverH = xp(iMaxQ)' "${source_file}"
    grep -Fq 'xFOverH = xp(iMaxF)' "${source_file}"
    grep -Fq "write(diagnosticUnit,'(ES24.16E3,1X,I0" "${source_file}"
}

write_case_metadata() {
    local case_dir=$1
    local case_name=$2
    local sq_value=$3
    local run_duration=$4
    local load_init=$5
    local nonfinite_interval=$6
    local source_file="${case_dir}/source/${SOURCE_NAME}"

    sha256sum "${STAGED_SOURCE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${source_file}" > "${case_dir}/source/case_source.sha256"
    {
        printf 'case=%s\nnode=node07\n' "${case_name}"
        printf 'Ra=1.0d8\nPr=0.7d0\nnx=513\nny=513\nchi_nu=0.4d0\nchi_b=0.0d0\n'
        printf 'flow_policy=FLOW_ODD_FIXED_SQ\nSq=%s\n' "${sq_value}"
        printf 'unsteadyRunDuration=%s\nloadInitField=%s\n' "${run_duration}" "${load_init}"
        printf 'population_and_macro_nonfinite_interval_tff=%s\n' "${nonfinite_interval}"
        printf 'flow_odd_interval_tff=0.1d0\n'
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

make_fresh_case() {
    local case_name=$1
    local job_name=$2
    local sq_value=$3
    local run_duration=$4
    local nonfinite_interval=$5
    local case_dir="${ROOT}/${case_name}"

    if [[ -e "${case_dir}" ]]; then
        printf 'Refusing to overwrite existing case: %s\n' "${case_dir}" >&2
        exit 3
    fi
    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    configure_source "${case_dir}/source/${SOURCE_NAME}" "${sq_value}" "${run_duration}" 0 "${nonfinite_interval}"
    write_case_metadata "${case_dir}" "${case_name}" "${sq_value}" "${run_duration}" 0 "${nonfinite_interval}"
    write_pbs "${case_dir}" "${job_name}"
}

make_restart_008() {
    local case_name=sq0.0800000_t1000
    local case_dir="${ROOT}/${case_name}"

    if [[ -e "${case_dir}" ]]; then
        printf 'Refusing to overwrite existing case: %s\n' "${case_dir}" >&2
        exit 3
    fi
    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    cp -a "${SEED_008}/." "${case_dir}/results/"
    cp "${SEED_008}/run.status" "${case_dir}/results/seed_300tff_run.status"
    cp "${SEED_008}/timing.txt" "${case_dir}/results/seed_300tff_timing.txt"
    cp "${SEED_008}/runtime_source.sha256" "${case_dir}/results/seed_300tff_runtime_source.sha256"
    # 旧诊断由 nvfortran 自动折成多行；原目录完整保留，新续算目录从 300 t_ff 开启单行诊断新分段。
    rm -f "${case_dir}/results/PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q5.dat"
    rm -f "${case_dir}/results/FlowOddMomentHistory_2DOpenaccLBMCDE_D2Q5.dat"
    configure_source "${case_dir}/source/${SOURCE_NAME}" 0.08d0 1000.0d0 1 1.0d0
    write_case_metadata "${case_dir}" "${case_name}" 0.08d0 1000.0d0 1 1.0d0
    printf 'restart_seed=%s\nrestart_seed_time_tf=299.99997829607514\n' "${SEED_008}" >> "${case_dir}/case_settings.txt"
    write_pbs "${case_dir}" SQ8C08T1000
}

make_fresh_case sq0.0388777_diag_t60 SQ8D0389 0.0388776875914752d0 60.0d0 0.1d0
make_fresh_case sq0.0500000_v2 SQ8C050V2 0.05d0 300.0d0 1.0d0
make_restart_008

diag_job=$(qsub "${ROOT}/sq0.0388777_diag_t60/pbs/run.pbs")
sq005_job=$(qsub -W depend=afterany:${diag_job} "${ROOT}/sq0.0500000_v2/pbs/run.pbs")
sq008_job=$(qsub -W depend=afterany:${sq005_job} "${ROOT}/sq0.0800000_t1000/pbs/run.pbs")

{
    printf 'sq0.0388777_diag_t60=%s\n' "${diag_job}"
    printf 'sq0.0500000_v2=%s\n' "${sq005_job}"
    printf 'sq0.0800000_t1000=%s\n' "${sq008_job}"
} | tee "${ROOT}/followup_submission_manifest.txt"
