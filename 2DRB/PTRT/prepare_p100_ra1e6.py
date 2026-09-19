"""Prepare a new immutable steady side-heated P100 case; does not submit."""
from pathlib import Path
import hashlib
import json
import re
import tarfile
import sys

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / '2DRBOpenaccPTRT.F90'
exponent = int(sys.argv[1]) if len(sys.argv) > 1 else 6
assert exponent in (6, 7, 8)
n = 256 if exponent == 6 else 512
NAME = f'Ra1e{exponent}/sideheated_N{n}_20260918'
REMOTE = '/data2/XLLi/PTRT/' + NAME
CASE = ROOT / 'p100' / NAME
CASE.mkdir(parents=True, exist_ok=False)
for folder in ('source', 'build', 'results'):
    (CASE / folder).mkdir()
raw = SOURCE.read_bytes()
text = raw.decode('utf-8-sig').replace('\r\n', '\n')
for macro, enabled in [('steadyFlow', True), ('unsteadyFlow', False),
        ('RayleighBenardCell', False), ('HorizontalWallsConstT', False),
        ('VerticalWallsAdiabatic', False), ('SideHeatedCell', True),
        ('HorizontalWallsAdiabatic', True), ('VerticalWallsConstT', True)]:
    text, count = re.subn(r'^!?#define ' + macro + r'[^\S\n]*$',
                         ('#define ' if enabled else '!#define ') + macro,
                         text, flags=re.M)
    assert count == 1, macro
for old, new in [('#define NX_OVERRIDE 1024', f'#define NX_OVERRIDE {n}'),
                 ('#define NY_OVERRIDE 1024', f'#define NY_OVERRIDE {n}'),
                 ('#define RAYLEIGH_OVERRIDE 10000000', f'#define RAYLEIGH_OVERRIDE {10**exponent}'),
                 ('Prandtl=0.7d0', 'Prandtl=0.71d0')]:
    assert text.count(old) == 1, old
    text = text.replace(old, new)
# Diagnostics only, needed to match the earlier multiblock acceptance workflow.
old = '        if(MOD(restartItcOffset+itc,2000).EQ.0) call check()'
assert text.count(old) == 1
text = text.replace(old, '''        if(MOD(restartItcOffset+itc,2000).EQ.0) then
            call check()
            call ptrt_case_check()
        endif''')
anchor = '''        if( (outputSnapshotFile.EQ.1).AND.(MOD(restartItcOffset+itc, outputSnapshotIntervalItc).EQ.0) ) then
            call update_host_snapshot_2d_openacc()'''
assert text.count(anchor) == 1
text = text.replace(anchor, anchor.replace('            call update_host_snapshot_2d_openacc()',
                                          '            call calNuRe()\n            call update_host_snapshot_2d_openacc()'))
anchor = '''#ifdef steadyFlow
    call output_Tecplot()'''
assert text.count(anchor) == 1
text = text.replace(anchor, '''#ifdef steadyFlow
    ! Preserve the actual stop state for independent final-checkpoint acceptance.
    call ptrt_case_check()
    call update_host_reload_2d_openacc()
    call output_ReloadFile()
    call output_Tecplot()''')
text += '''
! Case-only diagnostics: does not alter collision, streaming or wall treatment.
subroutine ptrt_case_check()
    use commondata
    use, intrinsic :: ieee_arithmetic
    implicit none
    !$acc wait(1)
    !$acc update self(u,v,T,rho)
    if (.not.all(ieee_is_finite(u))) error stop 86
    if (.not.all(ieee_is_finite(v))) error stop 86
    if (.not.all(ieee_is_finite(T))) error stop 86
    if (.not.all(ieee_is_finite(rho))) error stop 86
    if (.not.ieee_is_finite(errorU)) error stop 86
    if (.not.ieee_is_finite(errorT)) error stop 86
    if (minval(rho).le.0.0d0) error stop 86
    open(unit=89,file='Convergence_PTRT.dat',status='unknown',position='append')
    write(89,'(i12,3(1x,es24.16))') restartItcOffset+itc,errorU,errorT, &
        dble(restartItcOffset+itc)/timeUnit
    close(89)
end subroutine ptrt_case_check
'''
data = text.encode('utf-8')
(CASE / 'source/solver.F90').write_bytes(data)
(CASE / 'source/local_master.F90').write_bytes(raw)
master_hash = hashlib.sha256(raw).hexdigest()
case_hash = hashlib.sha256(data).hexdigest()
manifest = dict(local_source=str(SOURCE), remote_case=REMOTE,
    local_master_sha256=master_hash, case_source_sha256=case_hash,
    case='steady side-heated square cavity', flow='P-TRT conserved-moment Guo/MRT',
    thermal='Luo/Wang D2Q5', nx=n, ny=n, Ra=10**exponent, Pr=0.71, Ma=0.1,
    start='from_zero', epsU=1e-7, epsT=1e-7, check_steps=2000, max_steps=20000000,
    snapshot_tff=10, checkpoint_tff=100, tecplot_tff=100,
    node='node07', walltime='24:00:00', dependency=None,
    comparison=f'Multiblock Ra1e{exponent} sideheated_N{n}_20260916_v3; uniform grid here',
    kernel_changes=False, diagnostics=['finite macro fields and positive density every 2000 steps',
        'Nu/Re at snapshot times', 'final f/g checkpoint at actual stopping step',
        'Convergence_PTRT.dat (step,errorU,errorT,t_ff); final row may repeat'],
    acceptance='both errors <= 1e-7; final checkpoint finite/metrics audit remains separate')
(CASE / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
pbs = '''#!/bin/bash
#PBS -N PTRT_S6_N256
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=24:00:00
#PBS -o __CASE__/results/pbs.stdout
#PBS -e __CASE__/results/pbs.stderr
set -eu
CASE=__CASE__
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia
export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}
cd "$CASE/results"
test ! -e run.started
printf 'job=%s\nhost=%s\n' "${PBS_JOBID:-unknown}" "$(hostname)" > job_identity.txt
date -Is > timing.txt
sha256sum "$CASE/source/solver.F90" > runtime_source.sha256
test "$(cut -d ' ' -f1 runtime_source.sha256)" = '__HASH__'
nvfortran --version > compiler_version.txt 2>&1
printf '%s\n' '-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend' > build_flags.txt
cd "$CASE/build"
set +e
nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend "$CASE/source/solver.F90" -o solver.exe > "$CASE/results/compile.log" 2>&1
rc=$?
set -e
echo "$rc" > "$CASE/results/compile.status"
test "$rc" -eq 0 || exit "$rc"
cd "$CASE/results"
nvidia-smi > nvidia_smi_before.log
# Do not start GPU work over an unexpected existing GPU process.
nvidia-smi --query-compute-apps=pid --format=csv,noheader > gpu_pids_before.txt
if grep -Eq '[0-9]' gpu_pids_before.txt; then
    echo gpu_occupied > validation.status
    echo 88 > run.status
    exit 88
fi
date -Is > run.started
set +e
"$CASE/build/solver.exe" > solver.stdout 2>&1
rc=$?
set -e
echo "$rc" > run.status
date -Is > run.finished
date -Is >> timing.txt
if [ "$rc" -ne 0 ]; then echo runtime_failed > validation.status; exit "$rc"; fi
if grep -Eqi '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)|CUDA_ERROR|ERROR STOP|FATAL|Error:' solver.stdout ./*.dat SimulationSettings2DOpenacc.txt; then
    echo 86 > run.status; echo nonfinite_or_error > validation.status; exit 86
fi
if awk 'END {exit !(NF==4 && $2<=1e-7 && $3<=1e-7)}' Convergence_PTRT.dat; then
    echo steady_threshold_passed_pending_checkpoint_audit > validation.status
else
    echo steady_threshold_not_met > validation.status
    exit 89
fi
nvidia-smi > nvidia_smi_after.log
'''.replace('__CASE__', REMOTE).replace('__HASH__', case_hash).replace('PTRT_S6_N256', f'PTRT_S{exponent}_N{n}')
(CASE / 'run.pbs').write_text(pbs, encoding='utf-8', newline='\n')
with tarfile.open(ROOT / f'p100/ra1e{exponent}_submission.tar', 'w') as archive:
    archive.add(CASE, arcname=NAME)
print(json.dumps(manifest, indent=2))
