#!/usr/bin/env bash
set -euo pipefail

# 第一套算法 Ra=1e9 原始魔法参数 chi_nu 短扫描。
# chi_nu=0 已有长期结果，本脚本只提交 0.1--0.9、0.95、0.99。
# 所有算例均从静止初场运行到 100 t_ff；检测到非数时当前算例提前结束，
# 后续算例通过 afterany 依赖继续执行。已有目录绝不覆盖。

STAGE=/data2/XLLi/LBMCDE/FLOW-TEST/_staging/20260827_ra1e9_original_chi_t100/latest_local.F90
ROOT=/data2/XLLi/LBMCDE/FLOW-TEST/ORIGINAL/Ra1e9
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
DEPEND_AFTER=6317.master

CHIS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99)
CHI_FORTRAN=(0.1d0 0.2d0 0.3d0 0.4d0 0.5d0 0.6d0 0.7d0 0.8d0 0.9d0 0.95d0 0.99d0)
JOB_TAGS=(01 02 03 04 05 06 07 08 09 095 099)

test -f "${STAGE}"

for chi in "${CHIS[@]}"; do
    target=${ROOT}/chinu${chi}_t100_latest
    if [[ -e "${target}" ]]; then
        printf 'Refusing to overwrite %s\n' "${target}" >&2
        exit 3
    fi
done

configure_source() {
    local src=$1
    local chi=$2

    cp "${STAGE}" "${src}"
    sed -i 's/\r$//' "${src}"
    sed -i \
        -e 's/^!#define FLOW_ODD_ORIGINAL_MAGIC$/#define FLOW_ODD_ORIGINAL_MAGIC/' \
        -e 's/^#define FLOW_ODD_EFFECTIVE_MAGIC$/!#define FLOW_ODD_EFFECTIVE_MAGIC/' \
        -e 's/^#define FLOW_ODD_FIXED_SQ$/!#define FLOW_ODD_FIXED_SQ/' \
        -e 's/^#define FLOW_BGK$/!#define FLOW_BGK/' \
        -e 's/integer(kind=4), parameter :: nx=2048, ny=2048/integer(kind=4), parameter :: nx=1025, ny=1025/' \
        -e 's/real(kind=8), parameter :: Rayleigh=1.0d10/real(kind=8), parameter :: Rayleigh=1.0d9/' \
        -e "s/real(kind=8), parameter :: chi_nu=0.5d0/real(kind=8), parameter :: chi_nu=${chi}/" \
        -e 's/integer(kind=4), parameter :: loadInitField=1/integer(kind=4), parameter :: loadInitField=0/' \
        -e 's/real(kind=8), parameter :: unsteadyRunDuration=1000.0d0/real(kind=8), parameter :: unsteadyRunDuration=100.0d0/' \
        "${src}"

    grep -Fxq '#define FLOW_ODD_ORIGINAL_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_EFFECTIVE_MAGIC' "${src}"
    grep -Fxq '!#define FLOW_ODD_FIXED_SQ' "${src}"
    grep -Fxq '!#define FLOW_BGK' "${src}"
    grep -Fq 'parameter :: nx=1025, ny=1025' "${src}"
    grep -Fq 'parameter :: Rayleigh=1.0d9' "${src}"
    grep -Fq "parameter :: chi_nu=${chi}" "${src}"
    grep -Fq 'parameter :: loadInitField=0' "${src}"
    grep -Fq 'parameter :: unsteadyRunDuration=100.0d0' "${src}"
    grep -Fq 'parameter :: nonFiniteCheckIntervalTf=1.0d0' "${src}"
    grep -Fq 'parameter :: flowOddMomentDiagnosticIntervalTf=0.1d0' "${src}"
    grep -Fq 'Sq=1.0d0/(0.5d0+flowMagicParameter/(tauf-0.5d0))' "${src}"
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
#PBS -l walltime=12:00:00
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

make_case() {
    local chi_label=$1
    local chi_value=$2
    local job_tag=$3
    local case_dir=${ROOT}/chinu${chi_label}_t100_latest
    local job_name=O9C${job_tag}T100

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${case_dir}/results"
    configure_source "${case_dir}/source/${SOURCE_NAME}" "${chi_value}"
    sha256sum "${STAGE}" > "${case_dir}/source/local_master.sha256"
    sha256sum "${case_dir}/source/${SOURCE_NAME}" > "${case_dir}/source/case_source.sha256"
    cat > "${case_dir}/case_settings.txt" <<EOF
case=${job_name}
algorithm=D2Q9TRT_flow_D2Q5LuoTRT_temperature
Ra=1.0d9
Pr=0.7d0
nx=1025
ny=1025
chi_nu=${chi_value}
flow_policy=FLOW_ODD_ORIGINAL_MAGIC
flow_magic_parameter=3/16
Sq_formula=1/(0.5+(3/16)/(tau_0-0.5))
loadInitField=0
target_time_tf=100.0
nonfinite_macro_interval_tff=1.0
population_and_odd_moment_interval_tff=0.1
purpose=short_stability_screen_not_statistical_convergence
EOF
    write_pbs "${case_dir}" "${job_name}"
    printf '%s\n' "${case_dir}"
}

case_dirs=()
for index in "${!CHIS[@]}"; do
    case_dirs+=("$(make_case "${CHIS[index]}" "${CHI_FORTRAN[index]}" "${JOB_TAGS[index]}")")
done

previous=${DEPEND_AFTER}
for index in "${!case_dirs[@]}"; do
    jid=$(qsub -W depend=afterany:${previous} "${case_dirs[index]}/pbs/run.pbs")
    printf 'chi=%s job=%s depends_on=%s\n' "${CHIS[index]}" "${jid}" "${previous}"
    previous=${jid}
done
