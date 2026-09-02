#!/usr/bin/env bash
set -euo pipefail

# First solver (D2Q9-TRT flow + D2Q5 Luo-TRT temperature):
# 1) reproduce the late Sq=1 failure by restarting the committed t=2700 field;
# 2) screen Sq=0.05/0.08/0.20 at Ra=1e9, chi_nu=0.4 for 100 t_ff.
#
# The caller must upload the latest local source to STAGE before running this
# script. Existing case directories are never overwritten.

STAGE=/data2/XLLi/LBMCDE/FLOW-TEST/_staging/20260827_ra1e9_sq_diag/latest_local.F90
OLD=/data2/XLLi/LBMCDE/FLOW-TEST/FIXED_SQ1/Ra1e9/chinu0.0_t2800/results
RESTART=/data2/XLLi/LBMCDE/FLOW-TEST/FIXED_SQ1/Ra1e9/chinu0.0_restart2700_latest_t2800
SQROOT=/data2/XLLi/LBMCDE/FLOW-TEST/SQ_SCAN/Ra1e9/chinu0.4
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90

C05=${SQROOT}/sq0.0500000_t100
C08=${SQROOT}/sq0.0800000_t100
C20=${SQROOT}/sq0.2000000_t100

test -f "${STAGE}"
test -f "${OLD}/reloadFile2DOpenaccLBMCDE_D2Q5-000000000027.bin"
test -f "${OLD}/reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta"

for target in "${RESTART}" "${C05}" "${C08}" "${C20}"; do
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_source() {
    local src=$1
    local chi=$2
    local sq=$3
    local load=$4
    local duration=$5
    local nonfinite_interval=$6

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC$/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^!#define FLOW_ODD_FIXED_SQ$/#define FLOW_ODD_FIXED_SQ/' \
        -e "s|^[[:space:]]*real(kind=8), parameter :: Sq=1.0d0.*$|        real(kind=8), parameter :: Sq=${sq} ! fixed odd relaxation diagnostic|" \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=1025, ny=1025/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d9/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi}/" \
        -e "s/integer(kind=4), parameter :: loadInitField=0/integer(kind=4), parameter :: loadInitField=${load}/" \
        -e "s/real(kind=8), parameter :: nonFiniteCheckIntervalTf=1.0d0/real(kind=8), parameter :: nonFiniteCheckIntervalTf=${nonfinite_interval}/" \
        -e "s/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=${duration}/" \
        "${src}"

    grep -Fxq '!#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fq "parameter :: Sq=${sq}" "${src}"
    grep -Fq 'parameter :: nx=1025, ny=1025' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d9' "${src}"
    grep -Fq "parameter :: chi_nu=${chi}" "${src}"
    grep -Fq "parameter :: loadInitField=${load}" "${src}"
    grep -Fq "parameter :: nonFiniteCheckIntervalTf=${nonfinite_interval}" "${src}"
    grep -Fq "parameter :: unsteadyRunDuration=${duration}" "${src}"
    grep -Fq 'flowOddMomentDiagnosticIntervalTf=0.1d0' "${src}"
}

write_pbs() {
    local case_dir=$1
    local job_name=$2
    local src=${case_dir}/source/${SOURCE_NAME}
    local result=${case_dir}/results
    local pbs=${case_dir}/pbs/run.pbs

    cat > "${pbs}" <<'PBSEOF'
#!/usr/bin/env bash
#PBS -N __JOB__
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=04:00:00
#PBS -o __RESULT__/pbs.stdout
#PBS -e __RESULT__/pbs.stderr

set -u
SRC=__SOURCE__
RESULT_DIR=__RESULT__
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
NVFORTRAN=${NVHPC_ROOT}/compilers/bin/nvfortran

mkdir -p "${RESULT_DIR}"
cd "${RESULT_DIR}" || exit 90
export PATH="${NVHPC_ROOT}/compilers/bin:${PATH}"
export LD_LIBRARY_PATH="${NVHPC_ROOT}/compilers/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1
export ACC_DEVICE_TYPE=nvidia

printf 'PBS_JOBID=%s\n' "${PBS_JOBID:-unknown}" > job_identity.txt
printf 'HOSTNAME=%s\n' "$(hostname)" >> job_identity.txt
date '+CASE_START %F %T %z' > timing.txt
sha256sum "${SRC}" > runtime_source.sha256
"${NVFORTRAN}" --version > compiler_version.txt 2>&1
nvidia-smi > nvidia_smi_before.log 2>&1 || true

if [[ ! -x "${NVFORTRAN}" ]]; then
    printf '127\n' > compile.status
    printf '125\n' > run.status
    exit 127
fi

BUILD_FLAGS=(-cpp -O3 -acc -gpu=cc60 -Minfo=accel)
printf '%q ' "${BUILD_FLAGS[@]}" > build_flags.txt
printf '\n' >> build_flags.txt
"${NVFORTRAN}" "${BUILD_FLAGS[@]}" "${SRC}" -o solver.exe > compile.log 2>&1
compile_status=$?
printf '%s\n' "${compile_status}" > compile.status
if [[ "${compile_status}" -ne 0 ]]; then
    printf '125\n' > run.status
    date '+CASE_END %F %T %z' >> timing.txt
    exit "${compile_status}"
fi

set +e
./solver.exe > solver.stdout 2>&1
run_status=$?
set -e
printf '%s\n' "${run_status}" > run.status
nvidia-smi > nvidia_smi_after.log 2>&1 || true
date '+CASE_END %F %T %z' >> timing.txt
exit "${run_status}"
PBSEOF

    sed -i \
        -e "s|__JOB__|${job_name}|g" \
        -e "s|__SOURCE__|${src}|g" \
        -e "s|__RESULT__|${result}|g" \
        "${pbs}"
    chmod 700 "${pbs}"
}

mkdir -p "${RESTART}/source" "${RESTART}/pbs" "${RESTART}/results"
configure_source "${RESTART}/source/${SOURCE_NAME}" 0.0d0 1.0d0 1 2800.0d0 0.1d0
cp "${OLD}/reloadFile2DOpenaccLBMCDE_D2Q5-000000000027.bin" "${RESTART}/results/"
cp "${OLD}/reloadFile2DOpenaccLBMCDE_D2Q5-latest.meta" "${RESTART}/results/"
cp "${OLD}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat" "${RESTART}/results/"
cp "${OLD}/Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat" "${RESTART}/results/"
cp "${OLD}/DissipationHistory_2DOpenaccLBMCDE.dat" "${RESTART}/results/"
cp "${OLD}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin" "${RESTART}/results/"
sha256sum "${STAGE}" > "${RESTART}/source/local_master.sha256"
sha256sum "${RESTART}/source/${SOURCE_NAME}" > "${RESTART}/source/case_source.sha256"
cat > "${RESTART}/case_settings.txt" <<EOF
case=restart_from_t2700_latest_source
Ra=1.0d9
Pr=0.7d0
nx=1025
ny=1025
chi_nu=0.0d0
flow_policy=FLOW_ODD_FIXED_SQ
Sq=1.0d0
loadInitField=1
restart_seed=${OLD}/reloadFile2DOpenaccLBMCDE_D2Q5-000000000027.bin
restart_meta_time_tf=2699.9998817435312
target_time_tf=2800.0
nonfinite_interval_tff=0.1
odd_moment_interval_tff=0.1
EOF
write_pbs "${RESTART}" R9Q1R27

make_sq_case() {
    local case_dir=$1
    local sq=$2
    local job=$3
    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    configure_source "${case_dir}/source/${SOURCE_NAME}" 0.4d0 "${sq}" 0 100.0d0 0.1d0
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${case_dir}/source/${SOURCE_NAME}" > "${case_dir}/source/case_source.sha256"
    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job}
Ra=1.0d9
Pr=0.7d0
nx=1025
ny=1025
chi_nu=0.4d0
flow_policy=FLOW_ODD_FIXED_SQ
Sq=${sq}
loadInitField=0
target_time_tf=100.0
nonfinite_interval_tff=0.1
odd_moment_interval_tff=0.1
EOF
    write_pbs "${case_dir}" "${job}"
}

make_sq_case "${C05}" 0.05d0 S9C4Q05
make_sq_case "${C08}" 0.08d0 S9C4Q08
make_sq_case "${C20}" 0.20d0 S9C4Q20

jid_restart=$(qsub "${RESTART}/pbs/run.pbs")
jid05=$(qsub -W depend=afterany:${jid_restart} "${C05}/pbs/run.pbs")
jid08=$(qsub -W depend=afterany:${jid05} "${C08}/pbs/run.pbs")
jid20=$(qsub -W depend=afterany:${jid08} "${C20}/pbs/run.pbs")

printf 'restart=%s\nsq005=%s\nsq008=%s\nsq020=%s\n' \
    "${jid_restart}" "${jid05}" "${jid08}" "${jid20}"
