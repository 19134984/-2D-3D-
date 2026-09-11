#!/usr/bin/env bash
set -u

root=/data2/XLLi/LBMCDE/THRERMAL-D2Q9TRT-TEST

printf 'REL\tCOMPILE\tSEED\tRUN\tNCOUNT\tMAX_TF\tLAST_NU\tLAST_RE\tTARGET_TF\tPOLICY\tGRID\tABORT\tFORMAL\tDISS\n'

find "${root}" -type d -name results | grep 'Ra1e10' | sort | while IFS= read -r result; do
    rel=${result#${root}/}
    compile=NA
    seed=NA
    run=NA
    [[ -f "${result}/compile.status" ]] && compile=$(tr -d '\r\n ' < "${result}/compile.status")
    [[ -f "${result}/seed.status" ]] && seed=$(tr -d '\r\n ' < "${result}/seed.status")
    [[ -f "${result}/run.status" ]] && run=$(tr -d '\r\n ' < "${result}/run.status")

    nufile=$(find "${result}" -maxdepth 1 -type f -name 'Nu_VolAvg_*.dat' | head -n 1)
    refile=$(find "${result}" -maxdepth 1 -type f -name 'Re_VolAvg_*.dat' | head -n 1)
    ncount=0
    max_tf=NA
    last_nu=NA
    last_re=NA
    if [[ -n "${nufile}" && -s "${nufile}" ]]; then
        ncount=$(wc -l < "${nufile}" | tr -d ' ')
        read -r max_tf last_nu < <(awk 'NF>=2{t=$1;v=$2}END{print t,v}' "${nufile}")
    fi
    if [[ -n "${refile}" && -s "${refile}" ]]; then
        last_re=$(awk 'NF>=2{v=$2}END{print v}' "${refile}")
    fi

    case_dir=$(dirname "${result}")
    src=$(find "${case_dir}/source" -maxdepth 1 -type f -iname '*.F90' 2>/dev/null | head -n 1)
    target_tf=NA
    policy=NA
    grid=NA
    if [[ -n "${src}" && -s "${src}" ]]; then
        target_tf=$(grep -m1 -E 'parameter *:: *unsteadyRunDuration=' "${src}" | sed -E 's/.*unsteadyRunDuration=([^ !,]+).*/\1/' | tr -d '\r' || true)
        policy=$(grep -E '^#define FLOW_(ODD_ORIGINAL_MAGIC|ODD_EFFECTIVE_MAGIC|ODD_FIXED_SQ|BGK)' "${src}" | tr -d '\r' | tr '\n' ',' | sed 's/,$//' || true)
        nx=$(grep -m1 -E 'parameter *:: *nx=' "${src}" | sed -E 's/.*nx=([0-9]+).*/\1/' | tr -d '\r' || true)
        ny=$(grep -m1 -E 'parameter *:: *nx=[0-9]+, *ny=' "${src}" | sed -E 's/.*ny=([0-9]+).*/\1/' | tr -d '\r' || true)
        [[ -n "${nx:-}" && -n "${ny:-}" ]] && grid="${nx}x${ny}"
    fi

    abort=NONE
    if [[ -f "${result}/solver.stdout" ]]; then
        abort=$(grep 'NONFINITE_ABORT' "${result}/solver.stdout" | tail -n 1 | sed -E 's/[[:space:]]+/ /g;s/^ //;s/ /_/g' || true)
        [[ -z "${abort}" ]] && abort=NONE
    fi

    formal=NONE
    statfile=$(find "${result}" -maxdepth 1 -type f -name 'NuRe_TimeAverage_*.txt' | head -n 1)
    if [[ -n "${statfile}" && -s "${statfile}" ]]; then
        formal=$(awk '!/^#/ && NF>=10{line=$0}END{gsub(/[[:space:]]+/,"_",line);sub(/^_/,"",line);print line}' "${statfile}")
        [[ -z "${formal}" ]] && formal=NONE
    fi

    diss=NONE
    dissfile="${result}/NuReDissStatistics_2DOpenaccLBMCDE.dat"
    if [[ -s "${dissfile}" ]]; then
        diss=$(awk '
            $1=="Nu_volume_from_heat_flux"{nu=$2}
            $1=="Re_rms_space_time"{re=$2}
            $1=="Nu_wall_relative_difference_percent"{wall=$2}
            $1=="Nu_kinetic_relative_difference_percent"{kin=$2}
            $1=="Nu_thermal_relative_difference_percent"{therm=$2}
            $1=="R_u_calculated_over_exact"{ru=$2}
            $1=="R_T_calculated_over_exact"{rt=$2}
            END{printf("Nu=%s,Re=%s,wall=%s,kin=%s,therm=%s,Ru=%s,RT=%s",nu,re,wall,kin,therm,ru,rt)}
        ' "${dissfile}" | tr ' ' '_')
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${rel}" "${compile}" "${seed}" "${run}" "${ncount}" "${max_tf}" \
        "${last_nu}" "${last_re}" "${target_tf}" "${policy}" "${grid}" \
        "${abort}" "${formal}" "${diss}"
done
