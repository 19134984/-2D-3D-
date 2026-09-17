#!/bin/bash
set -eu
CASE=/data2/XLLi/Multiblock/Ra1e6/ratio2/sideheated_N256_20260916_v2
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia
cd "$CASE/results"
sha256sum "$CASE/source/solver.F90" > runtime_source.sha256
test "$(cut -d ' ' -f1 runtime_source.sha256)" = 'e31b2d317b901a481c647c65d249ea818a31b974d9a0726b181a5942e0671114'
nvfortran --version > compiler_version.txt 2>&1
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
date -Is > timing.txt
printf 'job=%s
host=%s
' "${PBS_JOBID:-unknown}" "$(hostname)" > job_identity.txt
nvidia-smi > nvidia_smi_before.log
set +e
"$CASE/results/solver.exe" > solver.stdout 2>&1
rc=$?
set -e
echo "$rc" > run.status
date -Is > run.finished
date -Is >> timing.txt
if [ "$rc" -ne 0 ]; then echo runtime_failed > validation.status; exit "$rc"; fi
if grep -Eqi '^[[:space:]]*(Error:|ERROR STOP|FATAL)' solver.stdout; then
    echo 87 > run.status; echo fatal_error > validation.status; exit 87
fi
if grep -Eqi '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)' solver.stdout ./*.dat; then
    echo 86 > run.status; echo nonfinite > validation.status; exit 86
fi
awk 'END {if(NF==3 && $2<=1e-7 && $3<=1e-7) print "steady_threshold_passed"; else print "steady_threshold_not_met"}' Convergence_2DOpenaccMultiblock.dat > validation.status
nvidia-smi > nvidia_smi_after.log
