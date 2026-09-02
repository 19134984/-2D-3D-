#!/usr/bin/env bash
set -euo pipefail

ROOT=/data2/XLLi/LBMCDE/THRERMAL-D2Q9TRT-TEST
STAGING=${ROOT}/_staging/20260830_ra1e10_original_equal_scan_t100_latest
MASTER=${STAGING}/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.local-master.F90
CAMPAIGN=${ROOT}/ORIGINAL_CHI_EQUAL_SCAN/Ra1e10_t100_latest

if [[ ! -f "${MASTER}" ]]; then
    echo "Missing uploaded local master: ${MASTER}" >&2
    exit 2
fi
if ! grep -q '^#define FLOW_ODD_ORIGINAL_MAGIC' "${MASTER}"; then
    echo 'Local master does not enable FLOW_ODD_ORIGINAL_MAGIC' >&2
    exit 3
fi
if grep -q '^#define FLOW_ODD_EFFECTIVE_MAGIC' "${MASTER}"; then
    echo 'Local master unexpectedly enables FLOW_ODD_EFFECTIVE_MAGIC' >&2
    exit 4
fi

mkdir -p "${CAMPAIGN}"
LOCAL_MASTER_HASH=$(sha256sum "${MASTER}" | awk '{print $1}')

submit_case() {
    local chi=$1
    local node=$2
    local dependency=$3
    local tag=${chi/./}
    local case_dir=${CAMPAIGN}/chinu${chi}_chikappa${chi}
    local src=${case_dir}/source/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90
    local pbs=${case_dir}/pbs/run.pbs

    if [[ -e "${case_dir}" ]]; then
        echo "Refusing to overwrite existing case: ${case_dir}" >&2
        exit 5
    fi

    mkdir -p "${case_dir}/source" "${case_dir}/results" "${case_dir}/pbs"
    cp "${MASTER}" "${src}"

    perl -0pi -e 's/integer\(kind=4\), parameter :: loadInitField=\d+/integer(kind=4), parameter :: loadInitField=0/' "${src}"
    perl -0pi -e 's/integer\(kind=4\), parameter :: nx=\d+, ny=\d+/integer(kind=4), parameter :: nx=2048, ny=2048/' "${src}"
    perl -0pi -e 's/real\(kind=8\), parameter :: Rayleigh=[^\r\n]*/real(kind=8), parameter :: Rayleigh=1.0d10/' "${src}"
    perl -0pi -e "s/(real\\(kind=8\\), parameter :: chi_nu=)[^!\\r\\n]*/\${1}${chi}d0          /" "${src}"
    perl -0pi -e "s/(real\\(kind=8\\), parameter :: chi_kappa=)[^!\\r\\n]*/\${1}${chi}d0       /" "${src}"
    perl -0pi -e 's/(real\(kind=8\), parameter :: unsteadyRunDuration=)[^!\r\n]*/${1}100.0d0  /' "${src}"

    grep -q '^#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    ! grep -q '^#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -q 'parameter :: loadInitField=0' "${src}"
    grep -q 'parameter :: nx=2048, ny=2048' "${src}"
    grep -q 'parameter :: Rayleigh=1.0d10' "${src}"
    grep -q "parameter :: chi_nu=${chi}d0" "${src}"
    grep -q "parameter :: chi_kappa=${chi}d0" "${src}"
    grep -q 'parameter :: unsteadyRunDuration=100.0d0' "${src}"

    {
        printf 'local_master_sha256 %s\n' "${LOCAL_MASTER_HASH}"
        printf 'case_source_sha256 '
        sha256sum "${src}" | awk '{print $1}'
        printf 'case Ra1e10 nx2048 ny2048 chi_nu=%s chi_kappa=%s policy=ORIGINAL duration_tf=100\n' "${chi}" "${chi}"
    } > "${case_dir}/source/source_provenance.txt"

    cat > "${pbs}" <<EOF
#!/usr/bin/env bash
#PBS -N T9O10C${tag}
#PBS -q batch
#PBS -l nodes=${node}:ppn=1
#PBS -l walltime=96:00:00
#PBS -o ${case_dir}/results/pbs_original_t100.stdout
#PBS -e ${case_dir}/results/pbs_original_t100.stderr

set -u
SRC=${src}
RESULT_DIR=${case_dir}/results
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
if [[ "\${run_status}" -eq 0 ]] && grep -Eqi '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)' \
    solver.stdout Nu_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat Re_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat \
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

    local job_id
    job_id=$(qsub -W "depend=afterany:${dependency}" "${pbs}")
    printf '%s|%s|%s|%s\n' "${job_id}" "${case_dir}" "${node}" "${dependency}" | tee -a "${STAGING}/submitted_jobs.txt"
    SUBMITTED_JOB_ID=${job_id}
}

: > "${STAGING}/submitted_jobs.txt"

dep05=6527.master
for chi in 0.0 0.2 0.4 0.6 0.8 0.99; do
    submit_case "${chi}" node05 "${dep05}"
    dep05=${SUBMITTED_JOB_ID}
done

dep07=6532.master
for chi in 0.1 0.3 0.5 0.7 0.9; do
    submit_case "${chi}" node07 "${dep07}"
    dep07=${SUBMITTED_JOB_ID}
done

echo "LOCAL_MASTER_SHA256=${LOCAL_MASTER_HASH}"
echo "NODE05_CHAIN_LAST=${dep05}"
echo "NODE07_CHAIN_LAST=${dep07}"
