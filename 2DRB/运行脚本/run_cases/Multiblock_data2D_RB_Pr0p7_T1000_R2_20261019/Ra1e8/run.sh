#!/bin/bash
set -euo pipefail
CASE=/data2/XLLi/Multiblock/data2D/RB/Ra1e8
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia CUDA_VISIBLE_DEVICES=0
cd "$CASE/results"
test "$(hostname -s)" = node07
test ! -e run.started
date -Is > run.started
finish() {
    rc=$?
    printf '%s\n' "$rc" > run.status
    date -Is > run.finished
    date -Is >> timing.txt
    trap - EXIT
    exit "$rc"
}
trap finish EXIT
date -Is > timing.txt
printf 'job=%s\nhost=%s\ncase=%s\n' "${PBS_JOBID:-unknown}" "$(hostname)" "$CASE" > job_identity.txt
sha256sum "$CASE/source/solver.F90" > runtime_source.sha256
test "$(cut -d ' ' -f1 runtime_source.sha256)" = 431e77597829829d3dd0d2b18df07bdc1a1cdff97acb435d365f65fa8496e185
nvfortran --version > compiler_version.txt 2>&1
printf '%s\n' '-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend' > build_flags.txt
set +e
nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend "$CASE/source/solver.F90" -o solver.exe > compile.log 2>&1
rc=$?
set -e
printf '%s\n' "$rc" > compile.status
if [ "$rc" -ne 0 ]; then printf 'compile_failed\n' > validation.status; exit "$rc"; fi
nvidia-smi > nvidia_smi_before.log 2>&1 || true
set +e
./solver.exe > solver.stdout 2>&1
rc=$?
set -e
if [ "$rc" -ne 0 ]; then printf 'runtime_failed\n' > validation.status; exit "$rc"; fi
if grep -Eqi '^[[:space:]]*(Error:|ERROR STOP|FATAL)' solver.stdout; then
    printf 'fatal_error\n' > validation.status; exit 87
fi
if ! python3 "$CASE/../serial_node07/validate_completion.py" > completion_validation.log 2>&1; then
    printf 'duration_validation_failed\n' > validation.status; exit 88
fi
nvidia-smi > nvidia_smi_after.log 2>&1 || true
