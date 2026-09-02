#!/usr/bin/env bash
set -euo pipefail

# Ra=1e8、513^2 网格的固定 Sq 稳定性与统计收敛测试。
# node07 串行执行；每个后续任务使用 afterany，因此前一任务因非数退出后仍会继续下一个参数。
ROOT=/data2/XLLi/LBMCDE/FLOW-TEST/SQ_SCAN/Ra1e8
STAGED_SOURCE="${ROOT}/_staging/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT_latest_local_20260827.F90"
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
SEED_005="${ROOT}/chinu0.4/sq0.0500000_t1000_v2/results"

if [[ ! -f "${STAGED_SOURCE}" ]]; then
    printf 'Missing latest local source: %s\n' "${STAGED_SOURCE}" >&2
    exit 2
fi

for required in \
    compile.status run.status \
    reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta \
    NuRe_InstantaneousVolAvg_2DOpenaccLBMCDE_D2Q5.plt \
    DissipationHistory_2DOpenaccLBMCDE.dat \
    TemperatureProfileHistory_2DOpenaccLBMCDE.bin; do
    if [[ ! -f "${SEED_005}/${required}" ]]; then
        printf 'Missing Sq=0.05 restart input: %s\n' "${SEED_005}/${required}" >&2
        exit 2
    fi
done

[[ "$(<"${SEED_005}/compile.status")" == 0 ]]
[[ "$(<"${SEED_005}/run.status")" == 0 ]]
grep -Fq 'time_tf  1.0000000401974960E+003' \
    "${SEED_005}/reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta"
grep -Fq 'cumulativeStatisticSampleCount 1001' \
    "${SEED_005}/reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta"

configure_source() {
    local source_file=$1
    local chi_value=$2
    local sq_value=$3
    local target_tf=$4
    local load_init=$5

    cp "${STAGED_SOURCE}" "${source_file}"
    sed -i 's/\r$//' "${source_file}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        -e "s/real(kind=8), parameter :: Sq=1.0d0 ! 固定奇模态松弛率对照：tau_q=1/real(kind=8), parameter :: Sq=${sq_value} ! 固定奇模态松弛率稳定性测试/" \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=513, ny=513/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d8/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi_value}/" \
        -e "s/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=${load_init}/" \
        -e "s/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=${target_tf}/" \
        "${source_file}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${source_file}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${source_file}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${source_file}"
    grep -Fxq '!#define FLOW_BGK' "${source_file}"
    grep -Fq "parameter :: Sq=${sq_value}" "${source_file}"
    grep -Fq 'parameter :: nx=513, ny=513' "${source_file}"
    grep -Fq 'parameter :: Rayleigh=1.0d8' "${source_file}"
    grep -Fq "parameter :: chi_nu=${chi_value}" "${source_file}"
    grep -Fq "parameter :: loadInitField=${load_init}" "${source_file}"
    grep -Fq "parameter :: unsteadyRunDuration=${target_tf}" "${source_file}"
    grep -Fq 'parameter :: nonFiniteCheckIntervalTf=1.0d0' "${source_file}"
    grep -Fq 'parameter :: flowOddMomentDiagnosticIntervalTf=0.1d0' "${source_file}"
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

make_case() {
    local chi_tag=$1
    local chi_value=$2
    local sq_tag=$3
    local sq_value=$4
    local target_tf=$5
    local load_init=$6
    local job_name=$7
    local seed_dir=$8
    local case_dir="${ROOT}/${chi_tag}/${sq_tag}_t${target_tf%%.*}"

    if [[ -e "${case_dir}" ]]; then
        printf 'Refusing to overwrite existing case: %s\n' "${case_dir}" >&2
        exit 3
    fi
    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"

    if [[ "${load_init}" -eq 1 ]]; then
        cp -a "${seed_dir}/." "${case_dir}/results/"
        cp "${seed_dir}/run.status" "${case_dir}/results/seed_1000tff_run.status"
        cp "${seed_dir}/timing.txt" "${case_dir}/results/seed_1000tff_timing.txt"
        cp "${seed_dir}/runtime_source.sha256" \
            "${case_dir}/results/seed_1000tff_runtime_source.sha256"
    fi

    configure_source "${case_dir}/source/${SOURCE_NAME}" \
        "${chi_value}" "${sq_value}" "${target_tf}" "${load_init}"
    sha256sum "${STAGED_SOURCE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${case_dir}/source/${SOURCE_NAME}" > "${case_dir}/source/case_source.sha256"
    {
        printf 'node=node07\nRa=1.0d8\nPr=0.7d0\nnx=513\nny=513\n'
        printf 'chi_nu=%s\nchi_b=0.0d0\n' "${chi_value}"
        printf 'flow_policy=FLOW_ODD_FIXED_SQ\nSq=%s\n' "${sq_value}"
        printf 'unsteadyRunDuration=%s\nloadInitField=%s\n' "${target_tf}" "${load_init}"
        printf 'restart_seed=%s\n' "${seed_dir:-none}"
        printf 'nonfinite_interval_tff=1.0d0\nflow_odd_interval_tff=0.1d0\n'
    } > "${case_dir}/case_settings.txt"
    write_pbs "${case_dir}" "${job_name}"
    printf '%s\n' "${case_dir}"
}

CASE_005=$(make_case chinu0.4 0.4d0 sq0.0500000 0.05d0 1400.0d0 1 SQ8C04Q05R "${SEED_005}")
CASE_04_05=$(make_case chinu0.4 0.4d0 sq0.5000000 0.5d0 1000.0d0 0 SQ8C04Q50 '')
CASE_04_10=$(make_case chinu0.4 0.4d0 sq1.0000000 1.0d0 1000.0d0 0 SQ8C04Q100 '')
CASE_05_05=$(make_case chinu0.5 0.5d0 sq0.5000000 0.5d0 1000.0d0 0 SQ8C05Q50 '')
CASE_05_10=$(make_case chinu0.5 0.5d0 sq1.0000000 1.0d0 1000.0d0 0 SQ8C05Q100 '')
CASE_09_05=$(make_case chinu0.9 0.9d0 sq0.5000000 0.5d0 1000.0d0 0 SQ8C09Q50 '')
CASE_09_10=$(make_case chinu0.9 0.9d0 sq1.0000000 1.0d0 1000.0d0 0 SQ8C09Q100 '')

tail_job=''
submit_afterany() {
    local case_dir=$1
    local job_id
    if [[ -n "${tail_job}" ]]; then
        job_id=$(qsub -W depend=afterany:${tail_job} "${case_dir}/pbs/run.pbs")
    else
        job_id=$(qsub "${case_dir}/pbs/run.pbs")
    fi
    tail_job=${job_id}
    printf '%s=%s\n' "$(basename "${case_dir}")" "${job_id}"
}

{
    submit_afterany "${CASE_005}"
    submit_afterany "${CASE_04_05}"
    submit_afterany "${CASE_04_10}"
    submit_afterany "${CASE_05_05}"
    submit_afterany "${CASE_05_10}"
    submit_afterany "${CASE_09_05}"
    submit_afterany "${CASE_09_10}"
} | tee "${ROOT}/high_sq_submission_manifest.txt"
