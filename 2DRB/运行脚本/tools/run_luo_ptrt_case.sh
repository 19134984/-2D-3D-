#!/bin/bash
set -euo pipefail
ROOT=/data2/XLLi/PTRT/LuoPTRT_Ra9_Ra10_20261025_r2
name=${1:?case required}
case "$name" in GPU_validation|P1_Ra1e9|P2_Ra1e9|P1_Ra1e10|P2_Ra1e10) ;; *) exit 2;; esac
test "$(hostname -s)" = node05
cd "$ROOT/$name/results"
RESULT_DIR="$ROOT/$name/results"
test ! -e run.started
date -Is > run.started
trap 'rc=$?; printf "%s\n" "$rc" > "$RESULT_DIR/run.status"; date -Is > "$RESULT_DIR/run.finished"' EXIT
printf 'job=%s\nhost=%s\ncase=%s\n' "${PBS_JOBID:-unknown}" "$(hostname)" "$name" > job_identity.txt
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia CUDA_VISIBLE_DEVICES=0
nvfortran --version > compiler_version.txt 2>&1
nvidia-smi > nvidia_smi_before.log 2>&1 || true
# Preserve other users' processes; start only when measured compute demand is low.
for attempt in $(seq 1 120); do
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits | head -1)
    test "$util" -lt 10 && break
    if test "$attempt" -eq 120; then echo 'GPU remained busy; no solver launched' >&2; exit 75; fi
    sleep 30
done
FLAGS=(-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend)
printf '%s\n' "${FLAGS[*]}" > build_flags.txt

# node05 only has Python 2.7. Lightweight assessment reads shared files on master.
# The solver and every GPU kernel continue to run exclusively inside this PBS node.
control_python() {
    ssh -o BatchMode=yes -o ConnectTimeout=10 master /usr/bin/python3 "$@"
}

build_and_run() {
    local label=$1
    cd "$ROOT/$label/results"
    sha256sum "$ROOT/$label/source/solver.F90" > runtime_source.sha256
    (cd "$ROOT" && sha256sum -c upload.sha256) > input_hash_verification.log
    set +e
    nvfortran "${FLAGS[@]}" "$ROOT/$label/source/solver.F90" -o solver.exe > compile.log 2>&1
    rc=$?
    printf '%s\n' "$rc" > compile.status
    set -e
    test "$rc" -eq 0 || return "$rc"
    set +e
    /usr/bin/time -p -o timing.txt ./solver.exe > solver.stdout 2>&1
    rc=$?
    printf '%s\n' "$rc" > solver.status
    set -e
    return "$rc"
}

if test "$name" = GPU_validation; then
    build_and_run GPU_P1
    build_and_run GPU_P2
    cd "$ROOT/GPU_validation/results"
    control_python "$ROOT/check_gpu_smoke.py" > gpu_comparison.log 2>&1
    exit 0
fi
test "$(cat "$ROOT/GPU_validation/results/run.status")" = 0
if [[ "$name" == *_Ra1e10 ]]; then
    previous="${name:0:2}_Ra1e9"
    if ! control_python "$ROOT/assess_case.py" "$previous" --gate > Ra9_gate.log 2>&1; then
        echo "Ra10 not launched: $previous has not passed the full Ra9 assessment" > NOT_LAUNCHED.txt
        exit 91
    fi
fi
expected=$(awk -v p="$name/source/solver.F90" '$2 == p {print $1}' "$ROOT/upload.sha256")
test "$(sha256sum "$ROOT/$name/source/solver.F90" | cut -d ' ' -f1)" = "$expected"
set +e
( build_and_run "$name" )
solver_rc=$?
set -e
cd "$ROOT/$name/results"
set +e
control_python "$ROOT/assess_case.py" "$name" > assessment.log 2>&1
assess_rc=$?
set -e
nvidia-smi > nvidia_smi_after.log 2>&1 || true
test "$solver_rc" -eq 0 || exit "$solver_rc"
exit "$assess_rc"
