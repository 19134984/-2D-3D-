#!/bin/bash
set -uo pipefail
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/bin:$PATH
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia CUDA_VISIBLE_DEVICES=0
cd "/data2/XLLi/Multiblock/data2D/RB/diagnosis_20261020/Ra8_replay/results" || exit 90
date -Is > run.started
sha256sum ../source/solver.F90 > runtime_source.sha256
nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend ../source/solver.F90 -o solver.exe > compile.log 2>&1
rc=$?; printf '%s\n' "$rc" > compile.status
if [ "$rc" -eq 0 ]; then ./solver.exe > solver.stdout 2>&1; rc=$?; fi
printf '%s\n' "$rc" > run.status
date -Is > run.finished
exit "$rc"
