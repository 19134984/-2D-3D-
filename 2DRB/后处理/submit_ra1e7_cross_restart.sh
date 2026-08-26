#!/usr/bin/env bash
set -euo pipefail

ROOT=/data2/XLLi/LBMCDE/FLOW-TEST/ORIGINAL/Ra1e7
CROSS_ROOT=${ROOT}/cross_restart
STAGING_SOURCE=${CROSS_ROOT}/_staging/latest_local.F90
PBS_TEMPLATE=${ROOT}/chinu0.5_fullhistory_t1400/pbs/run.pbs
RELOAD_PREFIX=reloadFile2DOpenaccLBMCDE_D2Q5

SOURCE_C05=${ROOT}/chinu0.5/results
SOURCE_C06=${ROOT}/chinu0.6/results
SOURCE_BIN_C05=${SOURCE_C05}/${RELOAD_PREFIX}-000000000028.bin
SOURCE_BIN_C06=${SOURCE_C06}/${RELOAD_PREFIX}-000000000028.bin
SOURCE_META_C05=${SOURCE_C05}/${RELOAD_PREFIX}-latest.meta
SOURCE_META_C06=${SOURCE_C06}/${RELOAD_PREFIX}-latest.meta

for required in \
    "${STAGING_SOURCE}" "${PBS_TEMPLATE}" \
    "${SOURCE_BIN_C05}" "${SOURCE_BIN_C06}" \
    "${SOURCE_META_C05}" "${SOURCE_META_C06}"; do
    test -f "${required}"
done

make_case() {
    local source_chi=$1
    local target_chi=$2
    local node=$3
    local job_name=$4
    local source_dir source_bin source_meta case_dir source_file result_dir

    case "${source_chi}" in
        0.5)
            source_dir=${SOURCE_C05}
            source_bin=${SOURCE_BIN_C05}
            source_meta=${SOURCE_META_C05}
            ;;
        0.6)
            source_dir=${SOURCE_C06}
            source_bin=${SOURCE_BIN_C06}
            source_meta=${SOURCE_META_C06}
            ;;
        *)
            printf 'Unsupported source chi: %s\n' "${source_chi}" >&2
            exit 2
            ;;
    esac

    case_dir=${CROSS_ROOT}/from_chinu${source_chi}_to_chinu${target_chi}
    source_file=${case_dir}/source/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90
    result_dir=${case_dir}/results

    if [[ -e "${case_dir}" ]]; then
        printf 'Refusing to overwrite existing case: %s\n' "${case_dir}" >&2
        exit 3
    fi

    mkdir -p "${case_dir}/source" "${case_dir}/pbs" "${result_dir}"
    cp "${STAGING_SOURCE}" "${source_file}"
    dos2unix "${source_file}" >/dev/null 2>&1 || true

    # 交叉重启使用相同的第一套算法，只修改网格、Ra、目标 chi_nu 和重启开关。
    sed -i \
        -e 's/loadInitField=0/loadInitField=1/' \
        -e 's/nx=2048, ny=2048/nx=257, ny=257/' \
        -e 's/Rayleigh=1.0d10/Rayleigh=1.0d7/' \
        -e "s/chi_nu=0.5d0/chi_nu=${target_chi}d0/" \
        "${source_file}"

    grep -q '#define FLOW_ODD_ORIGINAL_MAGIC' "${source_file}"
    grep -q 'loadInitField=1' "${source_file}"
    grep -q 'nx=257, ny=257' "${source_file}"
    grep -q 'Rayleigh=1.0d7' "${source_file}"
    grep -q "chi_nu=${target_chi}d0" "${source_file}"
    grep -q 'unsteadyRunDuration=1000.0d0' "${source_file}"

    cp "${source_bin}" "${result_dir}/${RELOAD_PREFIX}-000000000000.bin"
    cp "${source_meta}" "${result_dir}/input_restart_original.meta"
    cp "${source_meta}" "${result_dir}/${RELOAD_PREFIX}-latest.meta"

    # 这是“以成熟场为新初场”的参数突变实验：f/g 原样保留，实验时间轴从 0 重新计数。
    sed -i \
        -e 's/^reloadFileName .*/reloadFileName 000000000000/' \
        -e 's/^itc_total .*/itc_total 0/' \
        -e 's/^time_tf .*/time_tf 0.0000000000000000E+000/' \
        -e 's/^cumulativeStatisticSampleCount .*/cumulativeStatisticSampleCount 0/' \
        -e 's/^snapshotFileNum .*/snapshotFileNum 0/' \
        -e 's/^pltFileNum .*/pltFileNum 0/' \
        -e 's/^reloadFileNum .*/reloadFileNum 0/' \
        "${result_dir}/${RELOAD_PREFIX}-latest.meta"

    sha256sum "${source_file}" > "${case_dir}/source/source.sha256"
    sha256sum "${source_bin}" > "${result_dir}/input_restart.sha256"
    cp "${source_dir}/runtime_source_t2800.sha256" "${result_dir}/input_state_runtime_source.sha256"
    printf '%s\n' \
        "source_state_chi_nu=${source_chi}" \
        "target_chi_nu=${target_chi}" \
        "source_state_time_tf=2799.9996380293187" \
        "experiment_time_tf=0_to_1000" \
        "convergence_windows_tf=500_to_750,750_to_1000" \
        "final_statistics_tf=750_to_1000" \
        "restart_distribution_rescaled=no" \
        "source_restart=${source_bin}" \
        > "${result_dir}/cross_restart_provenance.txt"

    cp "${PBS_TEMPLATE}" "${case_dir}/pbs/run.pbs"
    sed -i \
        -e "s/#PBS -N L7O05FULL/#PBS -N ${job_name}/" \
        -e "s/nodes=node05:ppn=1/nodes=${node}:ppn=1/" \
        -e "s|${ROOT}/chinu0.5_fullhistory_t1400|${case_dir}|g" \
        "${case_dir}/pbs/run.pbs"
}

make_case 0.6 0.6 node05 L7R06T06
make_case 0.6 0.5 node05 L7R06T05
make_case 0.5 0.5 node07 L7R05T05
make_case 0.5 0.6 node07 L7R05T06

# 两条节点链都等待 chi=0.5、0.7 的全历史从零启动任务成功完成。
dependency=afterok:6251.master:6252.master
job_06_to_06=$(qsub -W depend="${dependency}" "${CROSS_ROOT}/from_chinu0.6_to_chinu0.6/pbs/run.pbs")
job_06_to_05=$(qsub -W depend="afterok:${job_06_to_06}" "${CROSS_ROOT}/from_chinu0.6_to_chinu0.5/pbs/run.pbs")
job_05_to_05=$(qsub -W depend="${dependency}" "${CROSS_ROOT}/from_chinu0.5_to_chinu0.5/pbs/run.pbs")
job_05_to_06=$(qsub -W depend="afterok:${job_05_to_05}" "${CROSS_ROOT}/from_chinu0.5_to_chinu0.6/pbs/run.pbs")

printf '%s\n' \
    "from_0.6_to_0.6=${job_06_to_06}" \
    "from_0.6_to_0.5=${job_06_to_05}" \
    "from_0.5_to_0.5=${job_05_to_05}" \
    "from_0.5_to_0.6=${job_05_to_06}" \
    | tee "${CROSS_ROOT}/submitted_jobs.txt"
