#!/usr/bin/env bash
set -u

root300=/data2/XLLi/LBMCDE/THRERMAL-D2Q9BGK-TEST/FLOW_BGK_CHI_EQUAL_SCAN/Ra1e10_restart100_latest_t300_grid2049_historyfix
root1000=/data2/XLLi/LBMCDE/THRERMAL-D2Q9BGK-TEST/FLOW_BGK_CHI_EQUAL_SCAN/Ra1e10_restart300_latest_t1000_grid2049_historyfix

printf 'stage\tchi\tstatus\tlast_nu_record\tlast_re_record\tstatistics_record\n'
for stage in 300 1000; do
    if [[ "$stage" == 300 ]]; then root="$root300"; else root="$root1000"; fi
    for chi in 0.0 0.2 0.5 0.7 0.9 0.99; do
        d="$root/chinu${chi}_chikappa${chi}/results"
        status=missing
        [[ -f "$d/run.status" ]] && status=$(tr -d '[:space:]' < "$d/run.status")
        nu=NA
        re=NA
        stat=NA
        if [[ -f "$d/Nu_VolAvg_2DOpenaccLBMCDE_D2Q9BGK.dat" ]]; then
            nu=$(awk 'NF && $1 !~ /^#/{line=$0} END{gsub(/^[[:space:]]+|[[:space:]]+$/,"",line); print line}' "$d/Nu_VolAvg_2DOpenaccLBMCDE_D2Q9BGK.dat")
        fi
        if [[ -f "$d/Re_VolAvg_2DOpenaccLBMCDE_D2Q9BGK.dat" ]]; then
            re=$(awk 'NF && $1 !~ /^#/{line=$0} END{gsub(/^[[:space:]]+|[[:space:]]+$/,"",line); print line}' "$d/Re_VolAvg_2DOpenaccLBMCDE_D2Q9BGK.dat")
        fi
        if [[ -f "$d/NuRe_TimeAverage_2DOpenaccLBMCDE_D2Q9BGK.txt" ]]; then
            stat=$(awk 'NF && $1 !~ /^#/{line=$0} END{gsub(/^[[:space:]]+|[[:space:]]+$/,"",line); print line}' "$d/NuRe_TimeAverage_2DOpenaccLBMCDE_D2Q9BGK.txt")
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$stage" "$chi" "$status" "$nu" "$re" "$stat"
    done
done
