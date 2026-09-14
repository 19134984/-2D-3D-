#!/usr/bin/env bash
set -u

emit_result_dir() {
    local family="$1"
    local d="$2"
    local status="missing"
    local nu_file=""
    local re_file=""
    local settings=""
    local meta=""
    local last_nu="NA"
    local last_re="NA"
    local meta_tf="NA"
    local chi_nu="NA"
    local chi_kappa="NA"
    local tau0="NA"
    local snu="NA"
    local sq="NA"
    local policy="NA"
    local fail="none"

    [[ -f "$d/run.status" ]] && status=$(tr -d '[:space:]' < "$d/run.status")
    nu_file=$(find "$d" -maxdepth 1 -type f -name 'Nu_VolAvg_*.dat' -print -quit 2>/dev/null || true)
    re_file=$(find "$d" -maxdepth 1 -type f -name 'Re_VolAvg_*.dat' -print -quit 2>/dev/null || true)
    settings=$(find "$d" -maxdepth 1 -type f -name 'SimulationSettings*.txt' -print -quit 2>/dev/null || true)
    meta=$(find "$d" -maxdepth 1 -type f -name 'reloadFile*-latest.meta' -print -quit 2>/dev/null || true)

    [[ -n "$nu_file" ]] && last_nu=$(awk 'NF && $1 !~ /^#/ && $1+0==$1 {v=$1} END{if(v!="") print v; else print "NA"}' "$nu_file")
    [[ -n "$re_file" ]] && last_re=$(awk 'NF && $1 !~ /^#/ && $1+0==$1 {v=$1} END{if(v!="") print v; else print "NA"}' "$re_file")
    [[ -n "$meta" ]] && meta_tf=$(awk '$1=="time_tf"{v=$2} END{if(v!="") print v; else print "NA"}' "$meta")

    if [[ -n "$settings" ]]; then
        chi_nu=$(awk -F'[=; ]+' '/chi_nu=/{for(i=1;i<=NF;i++) if($i=="chi_nu"){print $(i+1); exit}}' "$settings")
        chi_kappa=$(awk -F'[=; ]+' '/chi_kappa=/{for(i=1;i<=NF;i++) if($i=="chi_kappa"){print $(i+1); exit}}' "$settings")
        tau0=$(awk -F'[=; ]+' '/tau_0\/tauf=/{for(i=1;i<=NF;i++) if($i=="tau_0/tauf"){print $(i+1); exit}}' "$settings")
        snu=$(awk -F'[=; ]+' '/tau_0\/tauf=/{for(i=1;i<=NF;i++) if($i=="Snu"){print $(i+1); exit}}' "$settings")
        sq=$(awk '/tau_0\/tauf=/{getline; print $1; exit}' "$settings")
        policy=$(grep -m1 -E 'Flow (odd magic|collision) policy' "$settings" | sed 's/^[[:space:]]*//; s/[[:space:]]\+/ /g' || true)
    fi
    [[ -z "$chi_nu" ]] && chi_nu="NA"
    [[ -z "$chi_kappa" ]] && chi_kappa="NA"
    [[ -z "$tau0" ]] && tau0="NA"
    [[ -z "$snu" ]] && snu="NA"
    [[ -z "$sq" ]] && sq="NA"
    [[ -z "$policy" ]] && policy="NA"
    if [[ -f "$d/solver.stdout" ]]; then
        fail=$(grep -m1 'NONFINITE_ABORT' "$d/solver.stdout" | sed 's/^[[:space:]]*//; s/[[:space:]]\+/ /g' || true)
        [[ -z "$fail" ]] && fail="none"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$family" "$d" "$status" "$last_nu" "$last_re" "$meta_tf" \
        "$chi_nu" "$chi_kappa" "$tau0" "$snu" "$sq" "$policy" "$fail"
}

printf 'family\tresult_dir\trun_status\tlast_nu_tf\tlast_re_tf\tmeta_tf\tchi_nu\tchi_kappa\ttau0\tSnu\tSq\tpolicy\tnonfinite\n'

for branch in ORIGINAL EFFECTIVE FIXED_SQ1 FLOW_BGK; do
    root="/data2/XLLi/LBMCDE/FLOW-TEST/${branch}/Ra1e10"
    [[ -d "$root" ]] || continue
    while IFS= read -r d; do
        emit_result_dir "FIRST_${branch}" "$d"
    done < <(find "$root" -type d -name results | sort)
done

root2=/data2/XLLi/LBMCDE/THRERMAL-D2Q9BGK-TEST/FLOW_BGK_CHI_EQUAL_SCAN
while IFS= read -r d; do
    emit_result_dir "SECOND_FLOW_BGK" "$d"
done < <(find "$root2" -type d -path '*Ra1e10*' -name results | sort)
