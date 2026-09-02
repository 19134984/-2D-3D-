#!/usr/bin/env bash
set -euo pipefail

# Third solver: D2Q9-TRT flow + D2Q9-TRT temperature.
# Ra=1e10, original flow magic, chi_nu=chi_kappa=0.99.
# Continue the completed 0--100 t_ff case to an absolute target of 1000 t_ff.

ROOT=/data2/XLLi/LBMCDE/THRERMAL-D2Q9TRT-TEST
STAGING=${ROOT}/_staging/20260901_ra1e10_original_chi099_restart100_t1000_retry1
MASTER=${STAGING}/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.local-master.F90
SEED=${ROOT}/ORIGINAL_CHI_EQUAL_SCAN/Ra1e10_t100_latest/chinu0.99_chikappa0.99/results
CASE=${ROOT}/ORIGINAL_CHI_EQUAL_SCAN/Ra1e10_t100_latest/chinu0.99_chikappa0.99/continuation_t1000_latest_retry1
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90
NODE=node05

test -f "${MASTER}"
test "$(cat "${SEED}/compile.status")" = 0
test "$(cat "${SEED}/run.status")" = 0
if [[ -e "${CASE}" ]]; then
    printf 'Refusing to overwrite existing case: %s\n' "${CASE}" >&2
    exit 3
fi

META=${SEED}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta
test -f "${META}"
grep -Fxq 'reload_meta_version 4' "${META}"
grep -Fxq 'flowMode unsteadyFlow' "${META}"
grep -Fxq 'nx 2048' "${META}"
grep -Fxq 'ny 2048' "${META}"
grep -Fxq 'cumulativeStatisticSampleCount 101' "${META}"
restart_name=$(awk '$1=="reloadFileName" {print $2}' "${META}")
test -n "${restart_name}"
RESTART=${SEED}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-${restart_name}.bin
test -f "${RESTART}"

NU_HISTORY=Nu_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat
RE_HISTORY=Re_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat
DISS_HISTORY=DissipationHistory_2DOpenaccLBMCDE.dat
PROFILE_HISTORY=TemperatureProfileHistory_2DOpenaccLBMCDE.bin
POP_HISTORY=PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q9TRT.dat
test "$(wc -l < "${SEED}/${NU_HISTORY}")" -eq 101
test "$(wc -l < "${SEED}/${RE_HISTORY}")" -eq 101
test "$(wc -l < "${SEED}/${DISS_HISTORY}")" -eq 102
test -s "${SEED}/${PROFILE_HISTORY}"
test -s "${SEED}/${POP_HISTORY}"

mkdir -p "${CASE}/source" "${CASE}/pbs" "${CASE}/results"
SRC=${CASE}/source/${SOURCE_NAME}
RESULT=${CASE}/results
PBS=${CASE}/pbs/run.pbs

cp "${MASTER}" "${SRC}"
sed -i 's/\r$//' "${SRC}"
sed -i \
    -e 's/^!#define FLOW_ODD_ORIGINAL_MAGIC$/#define FLOW_ODD_ORIGINAL_MAGIC/' \
    -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
    -e 's/^#define FLOW_ODD_FIXED_SQ$/!#define FLOW_ODD_FIXED_SQ/' \
    -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
    -e 's/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=1/' \
    -e 's/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=0.99d0/' \
    -e 's/real(kind=8), parameter :: chi_kappa=0.5d0/real(kind=8), parameter :: chi_kappa=0.99d0/' \
    -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/' \
    "${SRC}"

grep -Fxq '#define FLOW_ODD_ORIGINAL_MAGIC' "${SRC}"
grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${SRC}"
grep -Fxq '!#define FLOW_ODD_FIXED_SQ' "${SRC}"
grep -Fxq '!#define FLOW_BGK' "${SRC}"
grep -Fq 'parameter :: loadInitField=1' "${SRC}"
grep -Fq 'parameter :: nx=2048, ny=2048' "${SRC}"
grep -Fq 'parameter :: Rayleigh=1.0d10' "${SRC}"
grep -Fq 'parameter :: chi_nu=0.99d0' "${SRC}"
grep -Fq 'parameter :: chi_kappa=0.99d0' "${SRC}"
grep -Fq 'parameter :: unsteadyRunDuration=1000.0d0' "${SRC}"

cp "${SEED}/${NU_HISTORY}" "${RESULT}/"
cp "${SEED}/${RE_HISTORY}" "${RESULT}/"
cp "${SEED}/${DISS_HISTORY}" "${RESULT}/"
cp "${SEED}/${PROFILE_HISTORY}" "${RESULT}/"
cp "${SEED}/${POP_HISTORY}" "${RESULT}/"
cp "${META}" "${RESULT}/"
cp "${RESTART}" "${RESULT}/"

local_master_hash=$(sha256sum "${MASTER}" | awk '{print $1}')
case_source_hash=$(sha256sum "${SRC}" | awk '{print $1}')
restart_hash=$(sha256sum "${RESTART}" | awk '{print $1}')
{
    printf 'local_master_sha256 %s\n' "${local_master_hash}"
    printf 'case_source_sha256 %s\n' "${case_source_hash}"
    printf 'restart_seed %s\n' "${SEED}"
    printf 'restart_file %s\n' "${restart_name}"
    printf 'restart_sha256 %s\n' "${restart_hash}"
    printf 'restart_time_tf 100\n'
    printf 'restart_sample_count 101\n'
    printf 'case Ra1e10 nx2048 ny2048 chi_nu=0.99 chi_kappa=0.99 policy=ORIGINAL target_tf=1000\n'
} > "${CASE}/source/source_provenance.txt"

cat > "${CASE}/case_settings.txt" <<EOF
case=T9O10C099R1K
algorithm=D2Q9TRT_flow_D2Q9TRT_temperature
Ra=1.0d10
Pr=0.7d0
nx=2048
ny=2048
chi_nu=0.99d0
chi_kappa=0.99d0
flow_policy=FLOW_ODD_ORIGINAL_MAGIC
loadInitField=1
target_time_tf=1000.0d0
restart_seed=${SEED}
restart_seed_time_tf=100
restart_seed_sample_count=101
statistics_window_tf=500--1000
half_windows_tf=500--750,750--1000
node=${NODE}
dependency=none
EOF

cat > "${PBS}" <<EOF
#!/usr/bin/env bash
#PBS -N T9O10C099R1K
#PBS -q batch
#PBS -l nodes=${NODE}:ppn=1
#PBS -l walltime=96:00:00
#PBS -o ${RESULT}/pbs_continue.stdout
#PBS -e ${RESULT}/pbs_continue.stderr

set -u
SRC=${SRC}
RESULT_DIR=${RESULT}
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
    solver.stdout "\${NU_HISTORY}" "\${RE_HISTORY}" "\${DISS_HISTORY}" 2>/dev/null; then
    printf 'non-finite value found in solver output/history\n' >> solver.stdout
    run_status=86
fi
printf '%s\n' "\${run_status}" > run.status
nvidia-smi > nvidia_smi_after.log 2>&1 || true
date '+CASE_END %F %T %z' >> timing.txt
exit "\${run_status}"
EOF
chmod 0755 "${PBS}"

job_id=$(qsub "${PBS}")
printf '%s|%s|%s|dependency=none\n' "${job_id}" "${CASE}" "${NODE}" | tee "${STAGING}/submitted_job.txt"
printf 'LOCAL_MASTER_SHA256=%s\n' "${local_master_hash}"
printf 'CASE_SOURCE_SHA256=%s\n' "${case_source_hash}"
printf 'JOB_ID=%s\n' "${job_id}"
