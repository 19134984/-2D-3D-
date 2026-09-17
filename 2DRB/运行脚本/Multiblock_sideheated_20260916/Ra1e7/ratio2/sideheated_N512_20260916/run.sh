#!/bin/bash
set -eu
CASE=/data2/XLLi/Multiblock/Ra1e7/ratio2/sideheated_N512_20260916
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia
cd "$CASE/build"
sha256sum "$CASE/source/solver.F90" > compiled_source.sha256
test "$(cut -d ' ' -f1 compiled_source.sha256)" = '3ec302a8e1102331a1f4833248604d5881330d9bf36644f370519a69d59c9b73'
nvfortran --version > compiler.version.txt 2>&1
printf '%s
' '-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend' > build_flags.txt
set +e
nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend "$CASE/source/solver.F90" -o solver.exe > compile.log 2>&1
rc=$?
set -e
echo "$rc" > compile.status
test "$rc" -eq 0 || exit "$rc"
if [ "${1:-}" = compile-only ]; then exit 0; fi
cd "$CASE/results"
test ! -e run.started
date -Is > run.started
printf 'job=%s
host=%s
' "${PBS_JOBID:-unknown}" "$(hostname)" > job_identity.txt
nvidia-smi > nvidia_smi_before.log
set +e
"$CASE/build/solver.exe" > solver.stdout 2>&1
rc=$?
set -e
echo "$rc" > run.status
date -Is > run.finished
if [ "$rc" -ne 0 ]; then echo runtime_failed > validation.status; exit "$rc"; fi
if grep -Eqi '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)' solver.stdout ./*.dat; then
    echo nonfinite > validation.status; exit 86
fi
awk 'END {if(NF==3 && $2<=1e-7 && $3<=1e-7) print "steady_threshold_passed"; else print "steady_threshold_not_met"}' Convergence_2DOpenaccMultiblock.dat > validation.status
nvidia-smi > nvidia_smi_after.log
