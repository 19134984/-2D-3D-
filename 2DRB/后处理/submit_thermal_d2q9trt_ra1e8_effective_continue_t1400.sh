#!/usr/bin/env bash
set -euo pipefail

ROOT=/data2/XLLi/LBMCDE/THRERMAL-D2Q9TRT-TEST
STAGING=${ROOT}/_staging/20260830_ra1e8_effective_continue_t1400
MASTER=${STAGING}/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.local-master.F90
CAMPAIGN=${ROOT}/EFFECTIVE_CHI_EQUAL_SCAN/Ra1e8

if [[ ! -f "${MASTER}" ]]; then
    echo "Missing uploaded local master: ${MASTER}" >&2
    exit 2
fi
if ! grep -q '^#define FLOW_ODD_ORIGINAL_MAGIC' "${MASTER}"; then
    echo 'Uploaded local master is not the expected latest source snapshot' >&2
    exit 3
fi

LOCAL_MASTER_HASH=$(sha256sum "${MASTER}" | awk '{print $1}')

submit_continuation() {
    local chi=$1
    local node=$2
    local dependency=$3
    local tag=${chi/./}
    local base=${CAMPAIGN}/chinu${chi}_chikappa${chi}
    local old_results=${base}/results
    local case_dir=${base}/continuation_t1400
    local src=${case_dir}/source/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90
    local results=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs
    local restart_name

    if [[ -e "${case_dir}" ]]; then
        echo "Refusing to overwrite existing continuation: ${case_dir}" >&2
        exit 4
    fi
    if [[ "$(cat "${old_results}/compile.status")" != 0 ]] || \
       [[ "$(cat "${old_results}/run.status")" != 0 ]]; then
        echo "Base run is not a completed calculation: ${base}" >&2
        exit 5
    fi

    test "$(wc -l < "${old_results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat")" -eq 1001
    test "$(wc -l < "${old_results}/Re_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat")" -eq 1001
    test "$(wc -l < "${old_results}/DissipationHistory_2DOpenaccLBMCDE.dat")" -eq 1002
    test "$(stat -c %s "${old_results}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin")" -eq 8224216
    grep -q '^cumulativeStatisticSampleCount 1000' \
        "${old_results}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta"
    restart_name=$(awk '$1=="reloadFileName" {print $2}' \
        "${old_results}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta")
    test -n "${restart_name}"
    test -f "${old_results}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-${restart_name}.bin"

    mkdir -p "${case_dir}/source" "${results}" "${case_dir}/pbs"
    cp "${MASTER}" "${src}"
    cp "${old_results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat" "${results}/"
    cp "${old_results}/Re_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat" "${results}/"
    cp "${old_results}/DissipationHistory_2DOpenaccLBMCDE.dat" "${results}/"
    cp "${old_results}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin" "${results}/"
    cp "${old_results}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta" "${results}/"
    cp "${old_results}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-${restart_name}.bin" "${results}/"

    # The old final metadata missed the t=1000 sample; all four histories contain t=0..1000.
    perl -0pi -e 's/^cumulativeStatisticSampleCount\s+1000\r?$/cumulativeStatisticSampleCount 1001/m' \
        "${results}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta"

    perl -0pi -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC/!#define FLOW_ODD_ORIGINAL_MAGIC/m; s/^!#define FLOW_ODD_EFFECTIVE_MAGIC/#define FLOW_ODD_EFFECTIVE_MAGIC/m' "${src}"
    perl -0pi -e 's/integer\(kind=4\), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=1/' "${src}"
    perl -0pi -e 's/integer\(kind=4\), parameter :: nx=\d+, ny=\d+/integer(kind=4), parameter :: nx=513, ny=513/' "${src}"
    perl -0pi -e 's/real\(kind=8\), parameter :: Rayleigh=[^\r\n]*/real(kind=8), parameter :: Rayleigh=1.0d8/' "${src}"
    perl -0pi -e "s/(real\\(kind=8\\), parameter :: chi_nu=)[^!\\r\\n]*/\${1}${chi}d0          /" "${src}"
    perl -0pi -e "s/(real\\(kind=8\\), parameter :: chi_kappa=)[^!\\r\\n]*/\${1}${chi}d0       /" "${src}"
    perl -0pi -e 's/(real\(kind=8\), parameter :: unsteadyRunDuration=)[^!\r\n]*/${1}1400.0d0  /' "${src}"

    grep -q '^!#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -q '^#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -q 'parameter :: loadInitField=1' "${src}"
    grep -q 'parameter :: nx=513, ny=513' "${src}"
    grep -q 'parameter :: Rayleigh=1.0d8' "${src}"
    grep -q "parameter :: chi_nu=${chi}d0" "${src}"
    grep -q "parameter :: chi_kappa=${chi}d0" "${src}"
    grep -q 'parameter :: unsteadyRunDuration=1400.0d0' "${src}"
    grep -q '^cumulativeStatisticSampleCount 1001' \
        "${results}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta"

    {
        printf 'local_master_sha256 %s\n' "${LOCAL_MASTER_HASH}"
        printf 'case_source_sha256 '
        sha256sum "${src}" | awk '{print $1}'
        printf 'restart_source %s\n' "${old_results}"
        printf 'restart_file %s\n' "${restart_name}"
        printf 'metadata_repair cumulativeStatisticSampleCount=1000_to_1001\n'
        printf 'case Ra1e8 nx513 ny513 chi_nu=%s chi_kappa=%s policy=EFFECTIVE target_tf=1400\n' "${chi}" "${chi}"
    } > "${case_dir}/source/source_provenance.txt"

    cat > "${pbs}" <<EOF
#!/usr/bin/env bash
#PBS -N T9E8C${tag}R
#PBS -q batch
#PBS -l nodes=${node}:ppn=1
#PBS -l walltime=96:00:00
#PBS -o ${results}/pbs_continue.stdout
#PBS -e ${results}/pbs_continue.stderr

set -u
SRC=${src}
RESULT_DIR=${results}
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

submit_continuation 0.0 node05 6511.master
node05_last=${SUBMITTED_JOB_ID}
submit_continuation 0.2 node05 "${node05_last}"
node05_last=${SUBMITTED_JOB_ID}

submit_continuation 0.1 node07 6515.master
node07_last=${SUBMITTED_JOB_ID}
submit_continuation 0.3 node07 "${node07_last}"
node07_last=${SUBMITTED_JOB_ID}

# Give the user's newest continuation request priority without cancelling the queued Ra1e7 scan.
qalter -W "depend=afterany:${node05_last}" 6523.master
qalter -W "depend=afterany:${node07_last}" 6528.master

echo "LOCAL_MASTER_SHA256=${LOCAL_MASTER_HASH}"
echo "NODE05_CONTINUATION_LAST=${node05_last}"
echo "NODE07_CONTINUATION_LAST=${node07_last}"
