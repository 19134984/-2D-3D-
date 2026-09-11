#!/usr/bin/env bash
# Submit the current D2Q5 Luo-TRT FLOW-TEST campaign without overwriting prior data.
set -euo pipefail

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
TEMPLATE="${ROOT}/analysis/current_d2q5_luo_trt_stats750_1000.F90"
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
MODE=${1:-all}

require_template() {
    if [[ ! -f "${TEMPLATE}" ]]; then
        printf 'missing current source template: %s\n' "${TEMPLATE}" >&2
        exit 80
    fi
}

write_pbs() {
    local source_file=$1
    local result_dir=$2
    local pbs_file=$3
    local node=$4
    local job_name=$5
    local walltime=$6

    cat > "${pbs_file}" <<EOF
#!/usr/bin/env bash
#PBS -N ${job_name}
#PBS -q batch
#PBS -l nodes=${node}:ppn=1
#PBS -l walltime=${walltime}
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
    local grid=$2
    local rayleigh=$3
    local chinu=$4
    local branch=$5
    local duration=$6
    local restart=$7

    cp "${TEMPLATE}" "${source_file}"
    sed -i \
        -e "s@integer(kind=4), parameter :: loadInitField=0@integer(kind=4), parameter :: loadInitField=${restart}@" \
        -e "s@integer(kind=4), parameter :: nx=2048, ny=2048@integer(kind=4), parameter :: nx=${grid}, ny=${grid}@" \
        -e "s@real(kind=8), parameter :: Rayleigh=1.0d10@real(kind=8), parameter :: Rayleigh=${rayleigh}@" \
        -e "s@real(kind=8), parameter :: chi_nu=0.5d0@real(kind=8), parameter :: chi_nu=${chinu}@" \
        -e "s@real(kind=8), parameter :: unsteadyRunDuration=1000.0d0@real(kind=8), parameter :: unsteadyRunDuration=${duration}@" \
        "${source_file}"
    if [[ "${branch}" == EFFECTIVE ]]; then
        sed -i \
            -e 's/^#define FLOW_ODD_ORIGINAL_MAGIC/!#define FLOW_ODD_ORIGINAL_MAGIC/' \
            -e 's/^!#define FLOW_ODD_EFFECTIVE_MAGIC/#define FLOW_ODD_EFFECTIVE_MAGIC/' \
            "${source_file}"
    fi
    grep -q "nx=${grid}, ny=${grid}" "${source_file}"
    grep -q "Rayleigh=${rayleigh}" "${source_file}"
    grep -q "chi_nu=${chinu}" "${source_file}"
    grep -q "unsteadyRunDuration=${duration}" "${source_file}"
    grep -q "loadInitField=${restart}" "${source_file}"
    if [[ "${branch}" == EFFECTIVE ]]; then
        grep -q '^#define FLOW_ODD_EFFECTIVE_MAGIC' "${source_file}"
    else
        grep -q '^#define FLOW_ODD_ORIGINAL_MAGIC' "${source_file}"
    fi
}

prepare_fresh_case() {
    local branch=$1
    local ra_label=$2
    local rayleigh=$3
    local grid=$4
    local chinu_label=$5
    local chinu_value=$6
    local node=$7
    local job_name=$8
    local walltime=$9
    local case_dir="${ROOT}/${branch}/${ra_label}/chinu${chinu_label}"
    local source_dir="${case_dir}/source"
    local result_dir="${case_dir}/results"
    local pbs_dir="${case_dir}/pbs"
    local source_file="${source_dir}/${SOURCE_NAME}"
    local pbs_file="${pbs_dir}/run.pbs"

    if [[ -e "${source_file}" || -e "${pbs_file}" || -e "${result_dir}/run.status" ]]; then
        printf 'refusing to overwrite prepared or completed case: %s\n' "${case_dir}" >&2
        exit 81
    fi
    mkdir -p "${source_dir}" "${result_dir}" "${pbs_dir}"
    configure_source "${source_file}" "${grid}" "${rayleigh}" "${chinu_value}" "${branch}" 1000.0d0 0
    sha256sum "${TEMPLATE}" > "${source_dir}/source_template.sha256"
    sha256sum "${source_file}" > "${source_dir}/source.sha256"
    {
        printf 'branch=%s\n' "${branch}"
        printf 'Ra=%s\n' "${rayleigh}"
        printf 'grid=%sx%s\n' "${grid}" "${grid}"
        printf 'chi_nu=%s\n' "${chinu_value}"
        printf 'flow_policy=%s\n' "${branch}"
        printf 'unsteady_run_duration_tff=1000\n'
        printf 'complete_history_window_tff=500_to_1000\n'
        printf 'final_statistics_window_tff=750_to_1000\n'
        printf 'node=%s\n' "${node}"
    } > "${case_dir}/case_manifest.txt"
    write_pbs "${source_file}" "${result_dir}" "${pbs_file}" "${node}" "${job_name}" "${walltime}"
    printf '%s\n' "${pbs_file}"
}

prepare_ra7_continuation() {
    local branch=$1
    local case_dir="${ROOT}/${branch}/Ra1e7/chinu0.0_t1400"
    local source_dir="${case_dir}/source"
    local result_dir="${case_dir}/results"
    local pbs_dir="${case_dir}/pbs"
    local source_file="${source_dir}/${SOURCE_NAME}"
    local pbs_file="${pbs_dir}/run.pbs"
    local old_results="${ROOT}/${branch}/Ra1e7/chinu0.0/results"
    local reload_prefix=reloadFile2DOpenaccLBMCDE_D2Q5
    local profile_record_bytes=$(( (1 + 2 * 257) * 8 ))

    if [[ -e "${source_file}" || -e "${pbs_file}" || -e "${result_dir}/run.status" ]]; then
        printf 'refusing to overwrite prepared or completed continuation: %s\n' "${case_dir}" >&2
        exit 82
    fi
    for path in \
        "${old_results}/${reload_prefix}-000000000009.bin" \
        "${old_results}/${reload_prefix}-latest.meta" \
        "${old_results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat" \
        "${old_results}/Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat" \
        "${old_results}/DissipationHistory_2DOpenaccLBMCDE.dat" \
        "${old_results}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin"; do
        [[ -f "${path}" ]] || { printf 'missing continuation input: %s\n' "${path}" >&2; exit 83; }
    done
    mkdir -p "${source_dir}" "${result_dir}" "${pbs_dir}"
    configure_source "${source_file}" 257 1.0d7 0.0d0 "${branch}" 1400.0d0 1
    cp "${old_results}/${reload_prefix}-000000000009.bin" "${result_dir}/"
    sed 's/^cumulativeStatisticSampleCount .*/cumulativeStatisticSampleCount 201/' \
        "${old_results}/${reload_prefix}-latest.meta" > "${result_dir}/${reload_prefix}-latest.meta"
    sed -n '201,401p' "${old_results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat" > \
        "${result_dir}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat"
    sed -n '201,401p' "${old_results}/Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat" > \
        "${result_dir}/Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat"
    sed -n '1p;202,402p' "${old_results}/DissipationHistory_2DOpenaccLBMCDE.dat" > \
        "${result_dir}/DissipationHistory_2DOpenaccLBMCDE.dat"
    dd if="${old_results}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin" \
        of="${result_dir}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin" \
        bs="${profile_record_bytes}" skip=200 count=201 status=none
    [[ $(wc -l < "${result_dir}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat") -eq 201 ]]
    [[ $(wc -l < "${result_dir}/Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat") -eq 201 ]]
    [[ $(wc -l < "${result_dir}/DissipationHistory_2DOpenaccLBMCDE.dat") -eq 202 ]]
    [[ $(stat -c %s "${result_dir}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin") -eq $((201 * profile_record_bytes)) ]]
    sha256sum "${TEMPLATE}" > "${source_dir}/source_template.sha256"
    sha256sum "${source_file}" > "${source_dir}/source.sha256"
    {
        printf 'branch=%s\n' "${branch}"
        printf 'Ra=1.0d7\n'
        printf 'grid=257x257\n'
        printf 'chi_nu=0.0d0\n'
        printf 'restart_field_origin=%s\n' "${old_results}/${reload_prefix}-000000000009.bin"
        printf 'strict_restart_time_tf=900.00179578872587\n'
        printf 'retained_history_tff=700_to_900\n'
        printf 'new_complete_history_tff=700_to_1400\n'
        printf 'new_final_statistics_tff=1050_to_1400\n'
        printf 'retained_history_samples=201\n'
    } > "${case_dir}/continuation_manifest.txt"
    write_pbs "${source_file}" "${result_dir}" "${pbs_file}" node07 "Luo${branch:0:1}_Ra1e7T14" 24:00:00
    printf '%s\n' "${pbs_file}"
}

submit_one() {
    local pbs_file=$1
    local parent_job=${2:-}
    local job_id
    if [[ -n "${parent_job}" ]]; then
        job_id=$(qsub -W "depend=afterany:${parent_job}" "${pbs_file}")
    else
        job_id=$(qsub "${pbs_file}")
    fi
    printf '%s\n' "${job_id}"
}

submit_high_ra() {
    local previous=""
    local spec branch ra_label rayleigh grid job_name pbs_file job_id
    local specs=(
        'ORIGINAL Ra1e9 1.0d9 1025 LuoO_Ra1e9'
        'EFFECTIVE Ra1e9 1.0d9 1025 LuoE_Ra1e9'
        'ORIGINAL Ra1e10 1.0d10 2049 LuoO_Ra1e10'
        'EFFECTIVE Ra1e10 1.0d10 2049 LuoE_Ra1e10'
    )
    for spec in "${specs[@]}"; do
        read -r branch ra_label rayleigh grid job_name <<< "${spec}"
        pbs_file=$(prepare_fresh_case "${branch}" "${ra_label}" "${rayleigh}" "${grid}" 0.0 0.0d0 node05 "${job_name}" 96:00:00)
        job_id=$(submit_one "${pbs_file}" "${previous}")
        printf '%s %s %s\n' "${branch}" "${ra_label}" "${job_id}"
        printf '%s\n' "${job_id}" > "${ROOT}/${branch}/${ra_label}/chinu0.0/pbs/job_id.txt"
        previous="${job_id}"
    done
}

submit_scan_ra6_ra8() {
    local previous=""
    local pbs_file job_id branch ra_label rayleigh grid chi
    for branch in ORIGINAL EFFECTIVE; do
        pbs_file=$(prepare_ra7_continuation "${branch}")
        job_id=$(submit_one "${pbs_file}" "${previous}")
        printf '%s Ra1e7_t1400 %s\n' "${branch}" "${job_id}"
        printf '%s\n' "${job_id}" > "${ROOT}/${branch}/Ra1e7/chinu0.0_t1400/pbs/job_id.txt"
        previous="${job_id}"
    done
    for ra_label in Ra1e6 Ra1e8; do
        if [[ "${ra_label}" == Ra1e6 ]]; then
            rayleigh=1.0d6
            grid=129
        else
            rayleigh=1.0d8
            grid=513
        fi
        for branch in ORIGINAL EFFECTIVE; do
            for chi in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
                pbs_file=$(prepare_fresh_case "${branch}" "${ra_label}" "${rayleigh}" "${grid}" "${chi}" "${chi}d0" node07 "Luo${branch:0:1}_${ra_label}_c${chi}" 24:00:00)
                job_id=$(submit_one "${pbs_file}" "${previous}")
                printf '%s %s chinu%s %s\n' "${branch}" "${ra_label}" "${chi}" "${job_id}"
                printf '%s\n' "${job_id}" > "${ROOT}/${branch}/${ra_label}/chinu${chi}/pbs/job_id.txt"
                previous="${job_id}"
            done
        done
    done
    printf '%s\n' "${previous}" > "${ROOT}/analysis/node07_ra6_ra8_scan_tail_job.txt"
}

require_template
case "${MODE}" in
    high)
        submit_high_ra
        ;;
    scan)
        submit_scan_ra6_ra8
        ;;
    all)
        submit_high_ra
        submit_scan_ra6_ra8
        ;;
    *)
        printf 'usage: %s [high|scan|all]\n' "$0" >&2
        exit 64
        ;;
esac
