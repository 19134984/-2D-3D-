"""Create immutable P100 side-heated cases from the current multiblock source."""
from pathlib import Path
import hashlib
import json
import re
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '运行脚本' / 'Multiblock_sideheated_20260916_v3'
BASE = '/data2/XLLi/Multiblock'
source = ROOT / '多块网格' / '2DRBOpenaccMultiblock.F90'
raw = source.read_bytes()
master_hash = hashlib.sha256(raw).hexdigest()
text = raw.decode('utf-8-sig')
assert 'EnableUseG' not in text
OUT.mkdir(exist_ok=True)
cases = []
for exponent, n in [(6, 256), (7, 512)]:
    name = f'Ra1e{exponent}/ratio2/sideheated_N{n}_20260916_v3'
    folder = OUT / name
    folder.mkdir(parents=True, exist_ok=False)
    for child in ['source', 'build', 'results']:
        (folder / child).mkdir()
    s = text
    for macro, enabled in [('steadyFlow', True), ('unsteadyFlow', False),
        ('RayleighBenardCell', False), ('HorizontalWallsConstT', False),
        ('VerticalWallsAdiabatic', False), ('SideHeatedCell', True),
        ('HorizontalWallsAdiabatic', True), ('VerticalWallsConstT', True)]:
        s, count = re.subn(r'^!?#define ' + macro + r'\s*$',
                          ('#define ' if enabled else '!#define ') + macro, s, flags=re.M)
        assert count == 1, (macro, count)
    changes = {'nx': str(n), 'ny': str(n), 'refineRatio': '2',
        'Rayleigh': f'1.0d{exponent}', 'Prandtl': '0.71d0',
        'fineLayerCellsLeft': str(n//8), 'fineLayerCellsBottom': str(n//8),
        'fineLayerCellsRight': str(n//8+1), 'fineLayerCellsTop': str(n//8+1)}
    for key, value in changes.items():
        s, count = re.subn(r'(\b' + key + r'\s*=\s*)[\d.]+(?:[dDeE][+-]?\d+)?',
                          lambda m: m[1]+value, s, count=1)
        assert count == 1, key
    s = s.replace('\r\n', '\n')
    data = s.encode('utf-8')
    (folder/'source/solver.F90').write_bytes(data)
    (folder/'source/local_master.F90').write_bytes(raw)
    case_hash = hashlib.sha256(data).hexdigest()
    meta = dict(case=name, local_master_sha256=master_hash, case_source_sha256=case_hash,
        ra=10**exponent, n=n, ratio=2, pr=0.71, mach=0.1, fine_indices=[n//8,n//8+1,n//8,n//8+1],
        epsU=1e-7, epsT=1e-7, max_steps=20000000, check_steps=2000,
        shared_gpu=True, source_changes=changes, numerical_kernels_changed=False,
        start='from_zero', node='node05', walltime='96:00:00',
        sampling_tff=10, checkpoint_tff=100, tecplot_tff=100,
        finite_check='owned macroscopic fields at Nu/Re sampling, every 10 tff',
        statistics='steady convergence; no unsteady half-window criterion')
    (folder/'manifest.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')
    run = '''#!/bin/bash
set -eu
CASE=__CASE__
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia
cd "$CASE/results"
sha256sum "$CASE/source/solver.F90" > runtime_source.sha256
test "$(cut -d ' ' -f1 runtime_source.sha256)" = '__HASH__'
nvfortran --version > compiler_version.txt 2>&1
printf '%s\n' '-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend' > build_flags.txt
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
printf 'job=%s\nhost=%s\n' "${PBS_JOBID:-unknown}" "$(hostname)" > job_identity.txt
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
'''.replace('__CASE__',BASE+'/'+name).replace('__HASH__',case_hash)
    (folder/'run.sh').write_text(run, encoding='utf-8', newline='\n')
    cases.append(meta)
pbs = '''#!/bin/bash
#PBS -N MBsideR2
#PBS -q batch
#PBS -l nodes=node05:ppn=1
#PBS -l walltime=96:00:00
#PBS -o /data2/XLLi/Multiblock/Ra1e6/ratio2/sideheated_N256_20260916_v3/results/pbs.stdout
#PBS -e /data2/XLLi/Multiblock/Ra1e6/ratio2/sideheated_N256_20260916_v3/results/pbs.stderr
set -u
cd /data2/XLLi/Multiblock/campaign_20260916_v3 || exit 90
date -Is > started
'''
for case in cases[:1]:
    pbs += f"bash '{BASE}/{case['case']}/run.sh'\nrc=$?\necho '{case['case']}' $rc >> case_exit_status.txt\nif [ \"$rc\" -ne 0 ]; then exit \"$rc\"; fi\n"
pbs += 'date -Is > finished\n'
(OUT/'campaign_20260916_v3').mkdir(exist_ok=True)
(OUT/'campaign_20260916_v3/run.pbs').write_text(pbs, encoding='utf-8', newline='\n')
(OUT/cases[0]['case']/'run.pbs').write_text(pbs, encoding='utf-8', newline='\n')
(OUT/'campaign_20260916_v3/manifest.json').write_text(json.dumps(cases, indent=2), encoding='utf-8')
with tarfile.open(OUT/'cases.tar', 'w') as archive:
    for name in ['Ra1e6','Ra1e7','campaign_20260916_v3']:
        archive.add(OUT/name, arcname=name)
print(json.dumps(cases, indent=2))
