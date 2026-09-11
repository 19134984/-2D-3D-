#!/usr/bin/env bash
set -euo pipefail

# Third solver: D2Q9-TRT flow + D2Q9-TRT temperature.
# Ra=1e10, chi_nu=chi_kappa, fresh 0--100 t_ff screens.

ROOT=/data2/XLLi/LBMCDE/THRERMAL-D2Q9TRT-TEST
STAGING=${ROOT}/_staging/20260904_ra1e10_bgk_fixedsq_equal_t100
MASTER=${STAGING}/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.local-master.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90

test -f "${MASTER}"
LOCAL_MASTER_HASH=$(sha256sum "${MASTER}" | awk '{print $1}')

case_dir_for() {
    local policy=$1
    local chi=$2
    local campaign
    case "${policy}" in
        BGK) campaign=FLOW_BGK_CHI_EQUAL_SCAN ;;
        FIXED_SQ1) campaign=FIXED_SQ1_CHI_EQUAL_SCAN ;;
        *) printf 'Unknown policy: %s\n' "${policy}" >&2; return 2 ;;
    esac
    printf '%s/%s/Ra1e10_t100_latest/chinu%s_chikappa%s' \
        "${ROOT}" "${campaign}" "${chi}" "${chi}"
}

for chi in 0.0 0.2 0.7 0.9; do
    test ! -e "$(case_dir_for BGK "${chi}")"
done
for chi in 0.0 0.2 0.5 0.7 0.9; do
    test ! -e "$(case_dir_for FIXED_SQ1 "${chi}")"
done

submit_case() {
    local policy=$1
    local chi=$2
    local node=$3
    local dependency=$4
    local job_name=$5
    local case_dir source result pbs case_hash job_id

    case_dir=$(case_dir_for "${policy}" "${chi}")
    source=${case_dir}/source/${SOURCE_NAME}
    result=${case_dir}/results
    pbs=${case_dir}/pbs/run.pbs

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${result}"
    cp "${MASTER}" "${source}"
    sed -i 's/\r$//' "${source}"

    case "${policy}" in
        BGK)
            sed -i \
                -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
                -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
                -e 's/^#define FLOW_ODD_FIXED_SQ$/!#define FLOW_ODD_FIXED_SQ/' \
                -e 's/^!#define FLOW_BGK$/#define FLOW_BGK/' \
                "${source}"
            ;;
        FIXED_SQ1)
            sed -i \
                -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
                -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
                -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
                -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
                "${source}"
            ;;
    esac

    sed -i \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi}d0/" \
        -e "s/real(kind=8), parameter :: chi_kappa=0.5d0/real(kind=8), parameter :: chi_kappa=${chi}d0/" \
        -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=100.0d0/' \
        "${source}"

    grep -Fq 'parameter :: loadInitField=0' "${source}"
    grep -Fq 'parameter :: nx=2048, ny=2048' "${source}"
    grep -Fq 'parameter :: Rayleigh=1.0d10' "${source}"
    grep -Fq "parameter :: chi_nu=${chi}d0" "${source}"
    grep -Fq "parameter :: chi_kappa=${chi}d0" "${source}"
    grep -Fq 'parameter :: unsteadyRunDuration=100.0d0' "${source}"
    if [[ "${policy}" == BGK ]]; then
        grep -Fxq '#define FLOW_BGK' "${source}"
        grep -Fxq '!#define FLOW_ODD_FIXED_SQ' "${source}"
    else
        grep -Fxq '!#define FLOW_BGK' "${source}"
        grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${source}"
    fi
    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${source}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${source}"

    case_hash=$(sha256sum "${source}" | awk '{print $1}')
    {
        printf 'local_master_sha256 %s\n' "${LOCAL_MASTER_HASH}"
        printf 'case_source_sha256 %s\n' "${case_hash}"
        printf 'case Ra1e10 nx2048 ny2048 chi_nu=%s chi_kappa=%s policy=%s duration_tf=100\n' \
            "${chi}" "${chi}" "${policy}"
    } > "${case_dir}/source/source_provenance.txt"

    cat > "${case_dir}/case_settings.txt" <<EOF
algorithm=D2Q9TRT_flow_D2Q9TRT_temperature
Ra=1.0d10
Pr=0.7d0
nx=2048
ny=2048
chi_nu=${chi}
chi_kappa=${chi}
flow_policy=${policy}
loadInitField=0
target_time_tf=100.0d0
statistics_window_tf=50--100
half_windows_tf=50--75,75--100
node=${node}
dependency=afterany:${dependency}
EOF

    cat > "${pbs}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=${node}:ppn=1
#PBS -l walltime=96:00:00
#PBS -o ${result}/pbs.stdout
#PBS -e ${result}/pbs.stderr

set -u
SRC=${source}
RESULT_DIR=${result}
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=\${NVHPC_ROOT}/compilers/bin/nvfortran

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
if [[ "\${run_status}" -eq 0 ]] && grep -Eqi \
    '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)' \
    solver.stdout \
    Nu_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat \
    Re_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat \
    DissipationHistory_2DOpenaccLBMCDE.dat 2>/dev/null; then
    printf 'non-finite value found in solver output/history\n' >> solver.stdout
    run_status=86
fi
printf '%s\n' "\${run_status}" > run.status
nvidia-smi > nvidia_smi_after.log 2>&1 || true
date '+CASE_END %F %T %z' >> timing.txt
exit "\${run_status}"
EOF
    chmod 0755 "${pbs}"

    job_id=$(qsub -W "depend=afterany:${dependency}" "${pbs}")
    printf '%s\n' "${job_id}" > "${case_dir}/pbs/job.id"
    printf '%s|%s|%s|afterany:%s|%s\n' \
        "${job_id}" "${policy}" "${chi}" "${dependency}" "${case_dir}" \
        >> "${STAGING}/submitted_jobs.txt"
    printf '%s\n' "${job_id}"
}

: > "${STAGING}/submitted_jobs.txt"

# node05: five new cases after the current node05 tail.
node05_tail=6668.master
node05_tail=$(submit_case BGK 0.0 node05 "${node05_tail}" T9B10C00)
node05_tail=$(submit_case BGK 0.7 node05 "${node05_tail}" T9B10C07)
node05_tail=$(submit_case FIXED_SQ1 0.2 node05 "${node05_tail}" T9Q10C02)
node05_tail=$(submit_case FIXED_SQ1 0.7 node05 "${node05_tail}" T9Q10C07)
node05_tail=$(submit_case FIXED_SQ1 0.9 node05 "${node05_tail}" T9Q10C09)

# node07: four new cases after the current node07 tail.
node07_tail=6670.master
node07_tail=$(submit_case BGK 0.2 node07 "${node07_tail}" T9B10C02)
node07_tail=$(submit_case BGK 0.9 node07 "${node07_tail}" T9B10C09)
node07_tail=$(submit_case FIXED_SQ1 0.0 node07 "${node07_tail}" T9Q10C00)
node07_tail=$(submit_case FIXED_SQ1 0.5 node07 "${node07_tail}" T9Q10C05)

printf 'LOCAL_MASTER_SHA256=%s\n' "${LOCAL_MASTER_HASH}"
printf 'NODE05_FINAL_TAIL=%s\n' "${node05_tail}"
printf 'NODE07_FINAL_TAIL=%s\n' "${node07_tail}"
cat "${STAGING}/submitted_jobs.txt"
