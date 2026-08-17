#!/usr/bin/env bash
set -euo pipefail

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
BASE_SOURCE="${ROOT}/ORIGINAL/Ra1e6/chinu0.0/source/current_source_template.F90"
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3

if [ ! -f "${BASE_SOURCE}" ]; then
    printf 'missing current source template: %s\n' "${BASE_SOURCE}" >&2
    exit 80
fi

setup_case() {
    local branch="$1"
    local ra_label="$2"
    local rayleigh="$3"
    local grid="$4"
    local node="$5"
    local job_name="$6"
    local case_dir="${ROOT}/${branch}/${ra_label}/chinu0.0"
    local source_dir="${case_dir}/source"
    local result_dir="${case_dir}/results"
    local pbs_dir="${case_dir}/pbs"
    local source_file="${source_dir}/${SOURCE_NAME}"
    local pbs_file="${pbs_dir}/run.pbs"

    if [ -e "${source_file}" ] || [ -e "${pbs_file}" ] || [ -e "${result_dir}/run.status" ]; then
        printf 'refusing to overwrite existing prepared or completed case: %s\n' "${case_dir}" >&2
        exit 81
    fi

    mkdir -p "${source_dir}" "${result_dir}" "${pbs_dir}"
    cp "${BASE_SOURCE}" "${source_file}"

    sed -i \
        -e "s@integer(kind=4), parameter :: nx=2048, ny=2048@integer(kind=4), parameter :: nx=${grid}, ny=${grid}@" \
        -e "s@real(kind=8), parameter :: Rayleigh=1.0d10@real(kind=8), parameter :: Rayleigh=${rayleigh}@" \
        -e "s@real(kind=8), parameter :: chi_nu=0.5d0@real(kind=8), parameter :: chi_nu=0.0d0@" \
        "${source_file}"

    case "${branch}" in
        ORIGINAL)
            ;;
        EFFECTIVE)
            sed -i \
                -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
                -e 's/^!#define FLOW_ODD_EFFECTIVE_MAGIC/#define FLOW_ODD_EFFECTIVE_MAGIC/' \
                "${source_file}"
            ;;
        *)
            printf 'unknown flow branch: %s\n' "${branch}" >&2
            exit 82
            ;;
    esac

    grep -q "nx=${grid}, ny=${grid}" "${source_file}"
    grep -q "Rayleigh=${rayleigh}" "${source_file}"
    grep -q "chi_nu=0.0d0" "${source_file}"
    if [ "${branch}" = EFFECTIVE ]; then
        grep -q '^#define FLOW_ODD_EFFECTIVE_MAGIC' "${source_file}"
    else
        grep -q '^#define FLOW_ODD_ORIGINAL_MAGIC' "${source_file}"
    fi

    sha256sum "${BASE_SOURCE}" > "${source_dir}/source_template.sha256"
    sha256sum "${source_file}" > "${source_dir}/source.sha256"
    {
        printf 'branch=%s\n' "${branch}"
        printf 'Ra=%s\n' "${rayleigh}"
        printf 'grid=%sx%s\n' "${grid}" "${grid}"
        printf 'chi_nu=0.0d0\n'
        printf 'flow_policy=%s\n' "${branch}"
        printf 'unsteady_run_duration_tff=1000.0d0\n'
        printf 'statistics_window_tff=500.0d0_to_1000.0d0\n'
        printf 'final_NuRe_window_tff=750.0d0_to_1000.0d0\n'
        printf 'node=%s\n' "${node}"
        printf 'source_file=%s\n' "${source_file}"
    } > "${case_dir}/case_manifest.txt"

    cat > "${pbs_file}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=${node}:ppn=1
#PBS -l walltime=72:00:00
#PBS -o ${result_dir}/pbs.stdout
#PBS -e ${result_dir}/pbs.stderr

set -u

CASE_DIR=${case_dir}
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
export ACC_NOTIFY=1

printf 'PBS_JOBID=%s\n' "\${PBS_JOBID:-unknown}" > job_identity.txt
printf 'HOSTNAME=%s\n' "\$(hostname)" >> job_identity.txt
date '+CASE_START %F %T %z' > timing.txt
sha256sum "\${SRC}" > runtime_source.sha256
printf '%s\n' "\${NVFORTRAN}" > compiler_path.txt
"\${NVFORTRAN}" --version > compiler_version.txt 2>&1
nvidia-smi > nvidia_smi_before.log 2>&1 || true

if [ ! -x "\${NVFORTRAN}" ]; then
    printf 'nvfortran unavailable: %s\n' "\${NVFORTRAN}" >&2
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
if [ "\${compile_status}" -ne 0 ]; then
    printf '125\n' > run.status
    exit "\${compile_status}"
fi

set +e
./solver.exe > solver.stdout 2>&1
run_status=\$?
set -e
printf '%s\n' "\${run_status}" > run.status
nvidia-smi > nvidia_smi_after.log 2>&1 || true
date '+CASE_END %F %T %z' >> timing.txt
exit "\${run_status}"
EOF
    chmod 700 "${pbs_file}"
}

setup_case ORIGINAL  Ra1e6 1.0d6 129 node05 LuoO_Ra1e6
setup_case EFFECTIVE Ra1e6 1.0d6 129 node07 LuoE_Ra1e6
setup_case ORIGINAL  Ra1e7 1.0d7 257 node05 LuoO_Ra1e7
setup_case EFFECTIVE Ra1e7 1.0d7 257 node07 LuoE_Ra1e7
setup_case ORIGINAL  Ra1e8 1.0d8 513 node05 LuoO_Ra1e8
setup_case EFFECTIVE Ra1e8 1.0d8 513 node07 LuoE_Ra1e8

job_o6=$(qsub "${ROOT}/ORIGINAL/Ra1e6/chinu0.0/pbs/run.pbs")
job_e6=$(qsub "${ROOT}/EFFECTIVE/Ra1e6/chinu0.0/pbs/run.pbs")
job_o7=$(qsub -W "depend=afterok:${job_o6}" "${ROOT}/ORIGINAL/Ra1e7/chinu0.0/pbs/run.pbs")
job_e7=$(qsub -W "depend=afterok:${job_e6}" "${ROOT}/EFFECTIVE/Ra1e7/chinu0.0/pbs/run.pbs")
job_o8=$(qsub -W "depend=afterok:${job_o7}" "${ROOT}/ORIGINAL/Ra1e8/chinu0.0/pbs/run.pbs")
job_e8=$(qsub -W "depend=afterok:${job_e7}" "${ROOT}/EFFECTIVE/Ra1e8/chinu0.0/pbs/run.pbs")

printf '%s\n' "${job_o6}" > "${ROOT}/ORIGINAL/Ra1e6/chinu0.0/pbs/job_id.txt"
printf '%s\n' "${job_e6}" > "${ROOT}/EFFECTIVE/Ra1e6/chinu0.0/pbs/job_id.txt"
printf '%s\n' "${job_o7}" > "${ROOT}/ORIGINAL/Ra1e7/chinu0.0/pbs/job_id.txt"
printf '%s\n' "${job_e7}" > "${ROOT}/EFFECTIVE/Ra1e7/chinu0.0/pbs/job_id.txt"
printf '%s\n' "${job_o8}" > "${ROOT}/ORIGINAL/Ra1e8/chinu0.0/pbs/job_id.txt"
printf '%s\n' "${job_e8}" > "${ROOT}/EFFECTIVE/Ra1e8/chinu0.0/pbs/job_id.txt"

printf 'ORIGINAL_Ra1e6 %s\n' "${job_o6}"
printf 'EFFECTIVE_Ra1e6 %s\n' "${job_e6}"
printf 'ORIGINAL_Ra1e7 %s afterok:%s\n' "${job_o7}" "${job_o6}"
printf 'EFFECTIVE_Ra1e7 %s afterok:%s\n' "${job_e7}" "${job_e6}"
printf 'ORIGINAL_Ra1e8 %s afterok:%s\n' "${job_o8}" "${job_o7}"
printf 'EFFECTIVE_Ra1e8 %s afterok:%s\n' "${job_e8}" "${job_e7}"
