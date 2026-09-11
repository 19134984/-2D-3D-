#!/usr/bin/env bash
set -euo pipefail

# Third solver: D2Q9-TRT flow + D2Q9-TRT temperature.
# Continue valid FLOW_BGK Ra=1e10 equal-chi screens from 100 to 300 t_ff.
# Restart files are copied and validated by the PBS job immediately before use.

ROOT=/data2/XLLi/LBMCDE/THRERMAL-D2Q9TRT-TEST
STAGING=${ROOT}/_staging/20260905_ra1e10_bgk_continue_t300_latest
MASTER=${STAGING}/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.local-master.F90
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90
NODE=node07

test -f "${MASTER}"
LOCAL_MASTER_HASH=$(sha256sum "${MASTER}" | awk '{print $1}')

parent_dir_for() {
    local chi=$1
    printf '%s/FLOW_BGK_CHI_EQUAL_SCAN/Ra1e10_t100_latest/chinu%s_chikappa%s' \
        "${ROOT}" "${chi}" "${chi}"
}

case_dir_for() {
    local chi=$1
    printf '%s/continuation_t300_latest' "$(parent_dir_for "${chi}")"
}

# The first two restart seeds are already complete.  The chi=0.7 seed is
# validated later by its PBS job because its 0--100 t_ff parent is still active.
for chi in 0.0 0.2; do
    seed=$(parent_dir_for "${chi}")/results
    test "$(cat "${seed}/compile.status")" = 0
    test "$(cat "${seed}/run.status")" = 0
    grep -Fxq 'reload_meta_version 4' \
        "${seed}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta"
    grep -Fxq 'cumulativeStatisticSampleCount 101' \
        "${seed}/reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta"
done

for chi in 0.0 0.2 0.7; do
    test ! -e "$(case_dir_for "${chi}")"
done

prepare_case() {
    local chi=$1
    local dependency=$2
    local job_name=$3
    local parent seed case_dir source result pbs case_hash

    parent=$(parent_dir_for "${chi}")
    seed=${parent}/results
    case_dir=$(case_dir_for "${chi}")
    source=${case_dir}/source/${SOURCE_NAME}
    result=${case_dir}/results
    pbs=${case_dir}/pbs/run.pbs

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${result}"
    cp "${MASTER}" "${source}"
    sed -i 's/\r$//' "${source}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^#define FLOW_ODD_FIXED_SQ$/!#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^!#define FLOW_BGK$/#define FLOW_BGK/' \
        -e 's/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=1/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi}d0/" \
        -e "s/real(kind=8), parameter :: chi_kappa=0.5d0/real(kind=8), parameter :: chi_kappa=${chi}d0/" \
        -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=300.0d0/' \
        "${source}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${source}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${source}"
    grep -Fxq '!#define FLOW_ODD_FIXED_SQ' "${source}"
    grep -Fxq '#define FLOW_BGK' "${source}"
    test "$(grep -Ec '^#define FLOW_(ODD_ORIGINAL_MAGIC|ODD_EFFECTIVE_MAGIC|ODD_FIXED_SQ|BGK)$' "${source}")" -eq 1
    grep -Fq 'parameter :: loadInitField=1' "${source}"
    grep -Fq 'parameter :: nx=2048, ny=2048' "${source}"
    grep -Fq 'parameter :: Rayleigh=1.0d10' "${source}"
    grep -Fq "parameter :: chi_nu=${chi}d0" "${source}"
    grep -Fq "parameter :: chi_kappa=${chi}d0" "${source}"
    grep -Fq 'parameter :: unsteadyRunDuration=300.0d0' "${source}"

    case_hash=$(sha256sum "${source}" | awk '{print $1}')
    {
        printf 'local_master_sha256 %s\n' "${LOCAL_MASTER_HASH}"
        printf 'case_source_sha256 %s\n' "${case_hash}"
        printf 'restart_seed %s\n' "${seed}"
        printf 'restart_expected_time_tf 100\n'
        printf 'restart_expected_sample_count 101\n'
        printf 'case Ra1e10 nx2048 ny2048 chi_nu=%s chi_kappa=%s policy=BGK target_tf=300\n' \
            "${chi}" "${chi}"
    } > "${case_dir}/source/source_provenance.txt"

    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job_name}
algorithm=D2Q9TRT_flow_D2Q9TRT_temperature
Ra=1.0d10
Pr=0.7d0
nx=2048
ny=2048
chi_nu=${chi}d0
chi_kappa=${chi}d0
flow_policy=FLOW_BGK
loadInitField=1
target_time_tf=300.0d0
restart_seed=${seed}
restart_expected_time_tf=100
restart_expected_sample_count=101
statistics_window_tf=150--300
half_windows_tf=150--225,225--300
node=${NODE}
dependency=${dependency}
EOF

    cat > "${pbs}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=${NODE}:ppn=1
#PBS -l walltime=96:00:00
#PBS -o ${result}/pbs_continue.stdout
#PBS -e ${result}/pbs_continue.stderr

set -u
SRC=${source}
SEED=${seed}
RESULT_DIR=${result}
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=\${NVHPC_ROOT}/compilers/bin/nvfortran
META_NAME=reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta
RELOAD_PREFIX=reloadFile2DOpenaccLBMCDE_D2Q9TRT
NU_HISTORY=Nu_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat
RE_HISTORY=Re_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat
DISS_HISTORY=DissipationHistory_2DOpenaccLBMCDE.dat
PROFILE_HISTORY=TemperatureProfileHistory_2DOpenaccLBMCDE.bin
POP_HISTORY=PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q9TRT.dat

cd "\${RESULT_DIR}" || exit 90
export PATH="\${NVHPC_ROOT}/compilers/bin:\${PATH}"
export LD_LIBRARY_PATH="\${NVHPC_ROOT}/compilers/lib:\${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1
export ACC_DEVICE_TYPE=nvidia

seed_fail() {
    printf '%s\n' "\$1" > seed.status
    printf 'restart seed validation failed with status %s\n' "\$1" > seed_error.log
    exit "\$1"
}

test -f "\${SEED}/compile.status" || seed_fail 91
test -f "\${SEED}/run.status" || seed_fail 91
test "\$(cat "\${SEED}/compile.status")" = 0 || seed_fail 91
test "\$(cat "\${SEED}/run.status")" = 0 || seed_fail 91
META=\${SEED}/\${META_NAME}
test -f "\${META}" || seed_fail 92
grep -Fxq 'reload_meta_version 4' "\${META}" || seed_fail 92
grep -Fxq 'flowMode unsteadyFlow' "\${META}" || seed_fail 92
grep -Fxq 'nx 2048' "\${META}" || seed_fail 92
grep -Fxq 'ny 2048' "\${META}" || seed_fail 92
grep -Fxq 'cumulativeStatisticSampleCount 101' "\${META}" || seed_fail 92
awk '\$1=="time_tf" && \$2>99.99 && \$2<100.01 {ok=1} END{exit(ok?0:1)}' "\${META}" || seed_fail 92

restart_name=\$(awk '\$1=="reloadFileName" {print \$2}' "\${META}")
test -n "\${restart_name}" || seed_fail 92
RESTART=\${SEED}/\${RELOAD_PREFIX}-\${restart_name}.bin
test -s "\${RESTART}" || seed_fail 92
test "\$(wc -l < "\${SEED}/\${NU_HISTORY}")" -eq 101 || seed_fail 93
test "\$(wc -l < "\${SEED}/\${RE_HISTORY}")" -eq 101 || seed_fail 93
test "\$(wc -l < "\${SEED}/\${DISS_HISTORY}")" -eq 102 || seed_fail 93
test -s "\${SEED}/\${PROFILE_HISTORY}" || seed_fail 93
test -s "\${SEED}/\${POP_HISTORY}" || seed_fail 93

for history in "\${NU_HISTORY}" "\${RE_HISTORY}" "\${DISS_HISTORY}" \
               "\${PROFILE_HISTORY}" "\${POP_HISTORY}"; do
    cp "\${SEED}/\${history}" . || seed_fail 94
done
cp "\${META}" . || seed_fail 94
cp "\${RESTART}" . || seed_fail 94
printf '0\n' > seed.status
{
    printf 'restart_seed %s\n' "\${SEED}"
    printf 'restart_file %s\n' "\${restart_name}"
    sha256sum "\${RESTART}"
    awk '\$1=="time_tf" || \$1=="cumulativeStatisticSampleCount" {print}' "\${META}"
} > restart_runtime_provenance.txt

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
    chmod 0755 "${pbs}"
}

prepare_case 0.0 none T9B10C00R3
prepare_case 0.2 pending T9B10C02R3
prepare_case 0.7 pending T9B10C07R3

: > "${STAGING}/submitted_jobs.txt"

case0=$(case_dir_for 0.0)
job0=$(qsub "${case0}/pbs/run.pbs")
printf '%s\n' "${job0}" > "${case0}/pbs/job.id"
printf '%s|BGK|0.0|dependency=none|%s\n' "${job0}" "${case0}" >> "${STAGING}/submitted_jobs.txt"

case2=$(case_dir_for 0.2)
job2=$(qsub -W "depend=afterany:${job0}" "${case2}/pbs/run.pbs")
printf '%s\n' "${job2}" > "${case2}/pbs/job.id"
sed -i "s/^dependency=pending$/dependency=afterany:${job0}/" "${case2}/case_settings.txt"
printf '%s|BGK|0.2|dependency=afterany:%s|%s\n' "${job2}" "${job0}" "${case2}" >> "${STAGING}/submitted_jobs.txt"

case7=$(case_dir_for 0.7)
job7=$(qsub -W "depend=afterany:${job2}" "${case7}/pbs/run.pbs")
printf '%s\n' "${job7}" > "${case7}/pbs/job.id"
sed -i "s/^dependency=pending$/dependency=afterany:${job2}/" "${case7}/case_settings.txt"
printf '%s|BGK|0.7|dependency=afterany:%s|%s\n' "${job7}" "${job2}" "${case7}" >> "${STAGING}/submitted_jobs.txt"

printf 'LOCAL_MASTER_SHA256=%s\n' "${LOCAL_MASTER_HASH}"
cat "${STAGING}/submitted_jobs.txt"
