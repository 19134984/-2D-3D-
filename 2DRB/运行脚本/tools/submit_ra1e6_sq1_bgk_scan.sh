#!/usr/bin/env bash
# Submit fresh Ra=1e6, 1000-t_ff scans for fixed Sq=1 and flow BGK.
set -euo pipefail

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
TEMPLATE=${ROOT}/analysis/current_d2q5_luo_trt_stats750_1000.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NODE=node07

write_pbs() {
    local source_file=$1
    local result_dir=$2
    local pbs_file=$3
    local job_name=$4

    cat > "${pbs_file}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=${NODE}:ppn=1
#PBS -l walltime=04:00:00
#PBS -o ${result_dir}/pbs.stdout
#PBS -e ${result_dir}/pbs.stderr

set -u
SRC=${source_file}
RESULT_DIR=${result_dir}
NVHPC_ROOT=${NVHPC_ROOT}
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
    solver.stdout Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat \
    DissipationHistory_2DOpenaccLBMCDE.dat 2>/dev/null; then
    printf 'non-finite value found in solver output/history\n' >> solver.stdout
    run_status=86
fi
printf '%s\n' "\${run_status}" > run.status
nvidia-smi > nvidia_smi_after.log 2>&1 || true
date '+CASE_END %F %T %z' >> timing.txt
exit "\${run_status}"
EOF
    chmod 700 "${pbs_file}"
}

configure_source() {
    local source_file=$1
    local chi_value=$2
    local policy=$3

    cp "${TEMPLATE}" "${source_file}"
    sed -i \
        -e 's@integer(kind=4), parameter :: nx=2048, ny=2048@integer(kind=4), parameter :: nx=129, ny=129@' \
        -e 's@real(kind=8), parameter :: Rayleigh=1.0d10@real(kind=8), parameter :: Rayleigh=1.0d6@' \
        -e "s@real(kind=8), parameter :: chi_nu=0.5d0@real(kind=8), parameter :: chi_nu=${chi_value}d0@" \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        "${source_file}"

    case "${policy}" in
        FIXED_SQ1)
            sed -i 's/^!#define FLOW_ODD_FIXED_SQ/#define FLOW_ODD_FIXED_SQ/' "${source_file}"
            grep -q 'real(kind=8), parameter :: Sq=1.0d0' "${source_file}"
            ;;
        FLOW_BGK)
            sed -i 's/^!#define FLOW_BGK/#define FLOW_BGK/' "${source_file}"
            grep -q 'real(kind=8), parameter :: Sq=Snu' "${source_file}"
            ;;
        *)
            printf 'unknown flow policy: %s\n' "${policy}" >&2
            exit 80
            ;;
    esac

    grep -q 'nx=129, ny=129' "${source_file}"
    grep -q 'Rayleigh=1.0d6' "${source_file}"
    grep -q "chi_nu=${chi_value}d0" "${source_file}"
    grep -q 'unsteadyRunDuration=1000.0d0' "${source_file}"
    grep -q 'subroutine check_nonfinite_state' "${source_file}"
    grep -q 'error stop 86' "${source_file}"
    if grep -qE 'rbInitPerturb|rbPerturbation' "${source_file}"; then
        printf 'unexpected initial perturbation code in %s\n' "${source_file}" >&2
        exit 81
    fi
}

prepare_case() {
    local policy=$1
    local chi_label=$2
    local case_dir=${ROOT}/${policy}/Ra1e6/chinu${chi_label}
    local source_dir=${case_dir}/source
    local result_dir=${case_dir}/results
    local pbs_dir=${case_dir}/pbs
    local source_file=${source_dir}/${SOURCE_NAME}
    local pbs_file=${pbs_dir}/run.pbs
    local job_name

    if find "${case_dir}" -mindepth 1 -type f -print -quit | grep -q .; then
        printf 'refusing to overwrite non-empty case: %s\n' "${case_dir}" >&2
        exit 82
    fi
    mkdir -p "${source_dir}" "${result_dir}" "${pbs_dir}"
    configure_source "${source_file}" "${chi_label}" "${policy}"
    sha256sum "${TEMPLATE}" > "${source_dir}/source_template.sha256"
    sha256sum "${source_file}" > "${source_dir}/source.sha256"

    if [[ "${policy}" == FIXED_SQ1 ]]; then
        job_name="SQ1_R6_c${chi_label}"
    else
        job_name="BGK_R6_c${chi_label}"
    fi
    {
        printf 'policy=%s\n' "${policy}"
        printf 'Ra=1.0d6\n'
        printf 'grid=129x129\n'
        printf 'chi_nu=%sd0\n' "${chi_label}"
        printf 'unsteady_run_duration_tff=1000\n'
        printf 'complete_history_window_tff=500_to_1000\n'
        printf 'final_statistics_window_tff=750_to_1000\n'
        printf 'nonfinite_guard=macro_fields_every_1_tff_and_statistics_before_write\n'
        printf 'node=%s\n' "${NODE}"
    } > "${case_dir}/case_manifest.txt"
    write_pbs "${source_file}" "${result_dir}" "${pbs_file}" "${job_name}"
    printf '%s\n' "${pbs_file}"
}

submit_case() {
    local pbs_file=$1
    local previous_job=${2:-}
    if [[ -n "${previous_job}" ]]; then
        qsub -W "depend=afterany:${previous_job}" "${pbs_file}"
    else
        qsub "${pbs_file}"
    fi
}

main() {
    local previous_job=''
    local chi policy pbs_file job_id case_dir

    [[ -f "${TEMPLATE}" ]]
    for chi in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
        for policy in FIXED_SQ1 FLOW_BGK; do
            pbs_file=$(prepare_case "${policy}" "${chi}")
            job_id=$(submit_case "${pbs_file}" "${previous_job}")
            case_dir=${ROOT}/${policy}/Ra1e6/chinu${chi}
            printf '%s\n' "${job_id}" > "${case_dir}/pbs/job_id.txt"
            printf '%s chi=%s %s\n' "${policy}" "${chi}" "${job_id}"
            previous_job=${job_id}
        done
    done
}

main "$@"
