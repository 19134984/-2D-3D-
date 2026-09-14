#!/usr/bin/env bash
set -euo pipefail
SEED=/data2/XLLi/LBMCDE/THRERMAL-D2Q9TRT-TEST/FLOW_BGK_CHI_EQUAL_SCAN/Ra1e10_t100_latest/chinu0.2_chikappa0.2/continuation_t300_latest/continuation_t1000_latest/results
CASE=/data2/XLLi/LBMCDE/THRERMAL-D2Q9TRT-TEST/FLOW_BGK_CHI_EQUAL_SCAN/Ra1e10_t100_latest/chinu0.2_chikappa0.2/continuation_t300_latest/continuation_t1000_latest/resume_t900_to1000_20260912
META=reloadFile2DOpenaccLBMCDE_D2Q9TRT-latest.meta
grep -q 'walltime .*exceeded limit' "$SEED/pbs_continue.stderr"
! grep -Eq 'NONFINITE_ABORT|^[[:space:]]*Error:' "$SEED/solver.stdout"
grep -Fxq 'nx 2048' "$SEED/$META"
grep -Fxq 'ny 2048' "$SEED/$META"
grep -Fxq 'cumulativeStatisticSampleCount 901' "$SEED/$META"
awk '$1=="time_tf" {ok=($2>899.99 && $2<900.01)} END{exit !ok}' "$SEED/$META"
test ! -e "$CASE"
mkdir -p "$CASE/source" "$CASE/results"
cp /data2/XLLi/LBMCDE/_resume_T9_20260912/solver.F90 "$CASE/source/"
cp /data2/XLLi/LBMCDE/_resume_T9_20260912/run.pbs "$CASE/"
dos2unix "$CASE/source/solver.F90" "$CASE/run.pbs" >/dev/null 2>&1
for name in Nu_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat Re_VolAvg_2DOpenaccLBMCDE_D2Q9TRT.dat DissipationHistory_2DOpenaccLBMCDE.dat; do
    awk '/^#/ {print;next} NF && $1+0<=900.00001 {print}' "$SEED/$name" > "$CASE/results/$name"
    awk '/^#/ {next} NF {if(tolower($0) ~ /nan|inf/) exit 1; if(($1-(n))>0.00001 || (($1-(n)) < -0.00001)) exit 2; n++} END {if(n!=901) exit 3}' "$CASE/results/$name"
done
# 原输出均保留；仅新目录截取前901条完整温度剖面，每条为(1+2*ny)个float64。
dd if="$SEED/TemperatureProfileHistory_2DOpenaccLBMCDE.bin" of="$CASE/results/TemperatureProfileHistory_2DOpenaccLBMCDE.bin" bs=32776 count=901 status=none
test "$(stat -c%s "$CASE/results/TemperatureProfileHistory_2DOpenaccLBMCDE.bin")" -eq 29531176
# 分布函数诊断历史由求解器按restart时刻恢复；保留完整副本作为输入。
cp "$SEED/PopulationNonequilibriumHistory_2DOpenaccLBMCDE_D2Q9TRT.dat" "$CASE/results/"
cp "$SEED/$META" "$CASE/results/"
BIN=reloadFile2DOpenaccLBMCDE_D2Q9TRT-000000000009.bin
cp "$SEED/$BIN" "$CASE/results/"
cmp "$SEED/$BIN" "$CASE/results/$BIN"
sha256sum "$SEED/$BIN" "$CASE/results/$BIN" > "$CASE/results/input_restart.sha256"
printf 'local_master_sha256 a2845f4920b680bc2fae33ce467de0d0f7a5821fabe9a05bfd69448596bd4051\nseed_tf 900\ntarget_tf 1000\n' > "$CASE/source/provenance.txt"
sha256sum "$CASE/source/solver.F90" >> "$CASE/source/provenance.txt"
bash -n "$CASE/run.pbs"
printf '%s\n' "$CASE"

