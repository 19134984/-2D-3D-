#!/usr/bin/env bash
# Repair the strict 1000 t_ff restart inputs, then requeue node07 serially.
set -euo pipefail

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST
PREFIX=reloadFile2DOpenaccLBMCDE_D2Q5
SOURCE_NAME=2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
RETAINED_SAMPLES=300  # 700 <= t/t_ff <= 999; the restart field writes t=1000 anew.
PROFILE_RECORD_BYTES=$(( (1 + 2 * 257) * 8 ))

archive_failed_attempt() {
    local results=$1
    local archive=$2
    mkdir -p "${archive}"
    local name
    for name in build_flags.txt commondata.mod compile.log compile.status compiler_path.txt compiler_version.txt \
        job_identity.txt nvidia_smi_before.log nvidia_smi_after.log pbs.stderr pbs.stdout run.status solver.exe \
        solver.stdout timing.txt SimulationSettings2DOpenaccLBMCDE_D2Q5.txt; do
        if [[ -e "${results}/${name}" ]]; then
            mv "${results}/${name}" "${archive}/"
        fi
    done
}

repair_continuation() {
    local branch=$1
    local old_results="${ROOT}/${branch}/Ra1e7/chinu0.0/results"
    local case_dir="${ROOT}/${branch}/Ra1e7/chinu0.0_t1400"
    local results="${case_dir}/results"
    local archive="${case_dir}/failed_restart_attempt_20260817"

    [[ -f "${old_results}/${PREFIX}-000000000010.bin" ]]
    [[ -f "${old_results}/${PREFIX}-latest.meta" ]]
    [[ -f "${case_dir}/source/${SOURCE_NAME}" ]]
    [[ -f "${case_dir}/pbs/run.pbs" ]]

    archive_failed_attempt "${results}" "${archive}"
    cp "${old_results}/${PREFIX}-000000000010.bin" "${results}/${PREFIX}-000000000010.bin"
    sed "s/^cumulativeStatisticSampleCount .*/cumulativeStatisticSampleCount ${RETAINED_SAMPLES}/" \
        "${old_results}/${PREFIX}-latest.meta" > "${results}/${PREFIX}-latest.meta"
    sed -n '201,500p' "${old_results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat" > \
        "${results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat"
    sed -n '201,500p' "${old_results}/Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat" > \
        "${results}/Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat"
    sed -n '1p;202,501p' "${old_results}/DissipationHistory_2DOpenaccLBMCDE.dat" > \
        "${results}/DissipationHistory_2DOpenaccLBMCDE.dat"
    dd if="${old_results}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin" \
        of="${results}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin" \
        bs="${PROFILE_RECORD_BYTES}" skip=200 count="${RETAINED_SAMPLES}" status=none

    [[ $(wc -l < "${results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat") -eq ${RETAINED_SAMPLES} ]]
    [[ $(wc -l < "${results}/Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat") -eq ${RETAINED_SAMPLES} ]]
    [[ $(wc -l < "${results}/DissipationHistory_2DOpenaccLBMCDE.dat") -eq $((RETAINED_SAMPLES + 1)) ]]
    [[ $(stat -c %s "${results}/TemperatureProfileHistory_2DOpenaccLBMCDE.bin") -eq $((RETAINED_SAMPLES * PROFILE_RECORD_BYTES)) ]]
    grep -q '^cumulativeStatisticSampleCount 300$' "${results}/${PREFIX}-latest.meta"
    grep -q '^reloadFileName 000000000010$' "${results}/${PREFIX}-latest.meta"
    grep -q '^ 7.0000000000000000E+002' "${results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat"
    grep -q '^ 9.9900000000000000E+002' "${results}/Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat"
    sed -i \
        -e 's/^strict_restart_time_tf=.*/strict_restart_time_tf=999.99987072475665/' \
        -e 's/^retained_history_tff=.*/retained_history_tff=700_to_999/' \
        -e 's/^retained_history_samples=.*/retained_history_samples=300/' \
        "${case_dir}/continuation_manifest.txt"
    printf 'restart_repair=20260817: copied strict checkpoint 000000000010 and retained 700_to_999 history\n' \
        >> "${case_dir}/continuation_manifest.txt"
}

submit_after_any() {
    local pbs_file=$1
    local parent=$2
    if [[ -n "${parent}" ]]; then
        qsub -W "depend=afterany:${parent}" "${pbs_file}"
    else
        qsub "${pbs_file}"
    fi
}

repair_continuation ORIGINAL
repair_continuation EFFECTIVE

previous=$(submit_after_any "${ROOT}/ORIGINAL/Ra1e7/chinu0.0_t1400/pbs/run.pbs" '')
printf 'ORIGINAL Ra1e7_t1400 %s\n' "${previous}"
printf '%s\n' "${previous}" > "${ROOT}/ORIGINAL/Ra1e7/chinu0.0_t1400/pbs/job_id.txt"
job_id=$(submit_after_any "${ROOT}/EFFECTIVE/Ra1e7/chinu0.0_t1400/pbs/run.pbs" "${previous}")
printf 'EFFECTIVE Ra1e7_t1400 %s\n' "${job_id}"
printf '%s\n' "${job_id}" > "${ROOT}/EFFECTIVE/Ra1e7/chinu0.0_t1400/pbs/job_id.txt"
previous=${job_id}

for ra_label in Ra1e6 Ra1e8; do
    for branch in ORIGINAL EFFECTIVE; do
        for chi in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
            pbs_file="${ROOT}/${branch}/${ra_label}/chinu${chi}/pbs/run.pbs"
            [[ -f "${pbs_file}" ]]
            job_id=$(submit_after_any "${pbs_file}" "${previous}")
            printf '%s %s chinu%s %s\n' "${branch}" "${ra_label}" "${chi}" "${job_id}"
            printf '%s\n' "${job_id}" > "${ROOT}/${branch}/${ra_label}/chinu${chi}/pbs/job_id.txt"
            previous=${job_id}
        done
    done
done
printf '%s\n' "${previous}" > "${ROOT}/analysis/node07_ra6_ra8_scan_tail_job.txt"
