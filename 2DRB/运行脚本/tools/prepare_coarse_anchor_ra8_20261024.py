"""Prepare protected-copy coarse-anchored and original-control RB GPU tests."""
from pathlib import Path
import hashlib, json, math, difflib, tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '运行脚本/run_cases/Multiblock_Ra1e8_coarse_anchor_N512_T30_20261024'
REMOTE = '/data2/XLLi/Multiblock/data2D/RB/diagnosis_coarse_anchor_N512_T30_20261024'
MASTER = ROOT / '多块网格/2DRBOpenaccMultiblock.F90'
raw = MASTER.read_bytes()
original = raw.decode('utf-8-sig').replace('\r\n', '\n')
assert not OUT.exists(), 'Never replace an existing campaign'
OUT.mkdir(parents=True)

def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding='utf-8', newline='\n')

def sha(data):
    return hashlib.sha256(data).hexdigest()

def once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)

baseline = once(original, 'nx=1024, ny=1024', 'nx=512, ny=512')
for side, old, new in [('Left',128,64), ('Right',129,65), ('Bottom',128,64), ('Top',129,65)]:
    baseline = once(baseline, f'fineLayerCells{side}={old}', f'fineLayerCells{side}={new}')
baseline = once(baseline, 'Rayleigh=1.0d7', 'Rayleigh=1.0d8')
baseline = once(baseline, 'unsteadyRunDuration=1000.0d0', 'unsteadyRunDuration=30.0d0')
assert '#define RayleighBenardCell\n' in baseline and '#define unsteadyFlow\n' in baseline
assert 'loadInitField=0' in baseline
assert (512-64-65+1)%2 == 0

anchor = once(baseline,
    'real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*sqrt(3.0d0*Prandtl/Rayleigh)',
    '! 独立诊断副本：先用全域等效粗网格长度计算粗格黏性松弛时间。\n'
    '        real(kind=8), parameter :: coarseLengthUnit=lengthUnit/dble(refineRatio)\n'
    '        real(kind=8), parameter :: taufCoarseAnchor=0.5d0+Mach*coarseLengthUnit*sqrt(3.0d0*Prandtl/Rayleigh)\n'
    '        ! 原始声学缩放得到细格参数；细格黏性、扩散、浮力和 TFF 标度与对照一致。\n'
    '        real(kind=8), parameter :: tauf=0.5d0+dble(refineRatio)*(taufCoarseAnchor-0.5d0)')
anchor = once(anchor, 'paraA=20.0d0*sqrt(3.0d0)*diffusivity-4.0d0',
    'paraA=20.0d0*sqrt(3.0d0)*(diffusivity/dble(refineRatio))-4.0d0')
start = anchor.index('        real(kind=8), parameter :: Snu=')
end = anchor.index('\n        !===============================================================================================', start)
anchor = anchor[:start] + '''        ! 本副本以粗格为基准：粗格流场满足原魔法参数，温度采用 Luo 固定参数。
        real(kind=8), parameter :: dxFine=1.0d0, dxCoarse=dble(refineRatio)
        real(kind=8), parameter :: SnuCoarse=1.0d0/taufCoarseAnchor
        real(kind=8), parameter :: SqCoarse=8.0d0*(2.0d0*taufCoarseAnchor-1.0d0)/(8.0d0*taufCoarseAnchor-1.0d0)
        real(kind=8), parameter :: QkCoarse=3.0d0-sqrt(3.0d0), QnuCoarse=4.0d0*sqrt(3.0d0)-6.0d0
        ! 所有细格松弛时间偏移量由粗格按原始粗细比缩放，不加正则化或联合反算。
        real(kind=8), parameter :: Snu=1.0d0/(0.5d0+dxCoarse*(1.0d0/SnuCoarse-0.5d0))
        real(kind=8), parameter :: Sq=1.0d0/(0.5d0+dxCoarse*(1.0d0/SqCoarse-0.5d0))
        real(kind=8), parameter :: Qk=1.0d0/(0.5d0+dxCoarse*(1.0d0/QkCoarse-0.5d0))
        real(kind=8), parameter :: Qnu=1.0d0/(0.5d0+dxCoarse*(1.0d0/QnuCoarse-0.5d0))
        ! paraA 按粗格扩散系数确定并在两级共用；保持两级物理扩散一致。
        real(kind=8), parameter :: thermalGeqCoeff=10.0d0/(4.0d0+paraA)''' + anchor[end:]

nu = .1*512*math.sqrt(3*.7/1e8)/3
alpha = nu/.7
tau = .5+3*nu
tau_c = .5+(tau-.5)/2
qk, qn = 3-math.sqrt(3), 4*math.sqrt(3)-6
sq = lambda t: 8*(2*t-1)/(8*t-1)
up = lambda s: 1/(.5+2*(1/s-.5))
down = lambda s: 1/(.5+(1/s-.5)/2)
parameters = {}
for name in ['original_control', 'coarse_anchor']:
    is_anchor = name == 'coarse_anchor'
    a = 20*math.sqrt(3)*alpha/(2 if is_anchor else 1)-4
    fine = [1/tau, up(sq(tau_c)) if is_anchor else sq(tau), up(qk) if is_anchor else qk, up(qn) if is_anchor else qn]
    coarse = [1/tau_c, sq(tau_c) if is_anchor else down(sq(tau)), qk if is_anchor else down(qk), qn if is_anchor else down(qn)]
    assert all(0<s<2 for s in fine+coarse)
    assert abs((4+a)/10*(1/fine[2]-.5)-alpha)<1e-15
    assert abs((4+a)/10*(1/coarse[2]-.5)-alpha/2)<1e-15
    assert abs((1/fine[0]-.5)/3-nu)<1e-15
    assert abs((1/coarse[0]-.5)/3-nu/2)<1e-15
    parameters[name] = dict(fine=fine, coarse=coarse, paraA=a, viscosity_fine=nu, diffusivity_fine=alpha,
        lambda_fine=(1/fine[0]-.5)*(1/fine[1]-.5), lambda_coarse=(1/coarse[0]-.5)*(1/coarse[1]-.5))

source_hashes = {}
for name, text in [('coarse_anchor',anchor), ('original_control',baseline)]:
    write(OUT/name/'source/solver.F90', text)
    (OUT/name/'results').mkdir()
    source_hashes[name] = sha(text.encode('utf-8'))
    write(OUT/name/'source/master_to_case.diff', ''.join(difflib.unified_diff(original.splitlines(True), text.splitlines(True), fromfile='protected_master', tofile=name)))
write(OUT/'coarse_anchor/source/control_to_anchor.diff', ''.join(difflib.unified_diff(baseline.splitlines(True), anchor.splitlines(True), fromfile='original_control', tofile='coarse_anchor')))
(OUT/'source_master').mkdir()
(OUT/'source_master/2DRBOpenaccMultiblock.F90').write_bytes(raw)

# Existing validator checks target, finite/contiguous history and final checkpoint alignment.
validator = (ROOT/'运行脚本/run_cases/Multiblock_Ra1e8_refined_N1024_T30_20261020/validate_completion.py').read_text(encoding='utf-8-sig')
write(OUT/'validate_completion.py', validator)
run = '''#!/bin/bash
set -euo pipefail
ROOT=REMOTE_ROOT
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia CUDA_VISIBLE_DEVICES=0
test "$(hostname -s)" = node07
test ! -e "$ROOT/campaign.started"
date -Is > "$ROOT/campaign.started"
pids=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits)
if [ -n "$pids" ]; then printf 'GPU occupied: %s\\n' "$pids" >&2; exit 75; fi
for name in coarse_anchor original_control; do
  CASE="$ROOT/$name"
  cd "$CASE/results"
  test ! -e run.started
  date -Is > run.started
  printf 'job=%s\\nhost=%s\\ncase=%s\\n' "${PBS_JOBID:-unknown}" "$(hostname)" "$CASE" > job_identity.txt
  sha256sum "$CASE/source/solver.F90" > runtime_source.sha256
  expected=$(python3 -c "import json; print(json.load(open('$ROOT/campaign.json'))['case_sha256']['$name'])")
  test "$(cut -d ' ' -f1 runtime_source.sha256)" = "$expected"
  nvfortran --version > compiler_version.txt 2>&1
  printf '%s\\n' '-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend' > build_flags.txt
  set +e
  nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend "$CASE/source/solver.F90" -o solver.exe > compile.log 2>&1
  rc=$?
  printf '%s\\n' "$rc" > compile.status
  if [ "$rc" -eq 0 ]; then
    /usr/bin/time -p -o timing.txt ./solver.exe > solver.stdout 2>&1
    rc=$?
    printf '%s\\n' "$rc" > solver.status
    python3 "$ROOT/assess.py" "$name" > assessment.log 2>&1
    arc=$?
    if [ "$rc" -eq 0 ] && [ "$arc" -ne 0 ]; then rc=$arc; fi
  fi
  printf '%s\\n' "$rc" > run.status
  date -Is > run.finished
  set -e
done
date -Is > "$ROOT/campaign.finished"
'''.replace('REMOTE_ROOT', REMOTE)
write(OUT/'run.sh',run)
write(OUT/'run.pbs',f'''#!/bin/bash
#PBS -N MB_Ra8_Canchor
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=02:00:00
#PBS -o {REMOTE}/pbs.stdout
#PBS -e {REMOTE}/pbs.stderr
set -euo pipefail
cd {REMOTE}
exec bash {REMOTE}/run.sh
''')
cfg = dict(remote_dir=REMOTE, master=str(MASTER), master_sha256=sha(raw), case_sha256=source_hashes,
    node='node07', order=['coarse_anchor','original_control'], Ra=1e8, Pr=.7, Ma=.1,
    fine_equivalent_N=512, coarse_equivalent_N=256, refineRatio=2,
    fine_layers=[64,65,64,65], center_fine_span=384, center_coarse_intervals=192,
    target_tff=30, history_snapshot_interval_tff=.5, loadInitField=0,
    initial_condition='unchanged conductive temperature, no added perturbation',
    filter='none', interface='unchanged original moment scaling/interpolation', parameters=parameters)
write(OUT/'campaign.json',json.dumps(cfg,indent=2)+'\n')
protected = [ROOT/p for p in ['多块网格/2DRBOpenaccMultiblock.F90','均匀网格/2DRBOpenmp.F90', '均匀网格/3DRBOpenmp.F90','均匀网格/3DRBOpenacc.F90','均匀网格/3DRBOpenaccMpi.F90']]
write(OUT/'protected_sources.json',json.dumps({str(f.relative_to(ROOT)):sha(f.read_bytes()) for f in protected},indent=2,ensure_ascii=False)+'\n')
write(OUT/'README.md', '''# Ra=1e8 粗网格基准原始缩放短测

用户确认：粗格满足原魔法参数和 Luo 温度参数，再通过原始声学缩放得到细格参数。只改副本。
Ra=1e8、Pr=0.7、Ma=0.1，细格全域等效512、粗格等效256，粗细比2；左/下64、右/上65，中心384细间隔=192粗间隔。
两个独立初场算例在 node07 串行执行，先 coarse_anchor，再 original_control，均从零算30 TFF；原初场不加扰动。
保留0.5 TFF快照与历史，原100 TFF检查点/Tecplot周期及结束补输出，接口和碰撞核不改，不加正则化。

粗格 tauf=.5+Ma*(Lfine/r)*sqrt(3Pr/Ra)，细格 tauf=.5+r*(tauf_coarse-.5)，与对照黏性相同。
coarse_anchor：粗格 Lambda=3/16，Qk=3-sqrt(3)，Qnu=4sqrt(3)-6；细格所有 1/s-.5 乘 r。
paraA=20sqrt(3)*(diffusivity_fine/r)-4，两级共用；细格有效扩散与对照相同，粗格为细格的1/r。
这会改变高阶矩和热平衡矩参数，不能称为仅改变运算顺序。接口继续使用其对应的各级实际松弛率。
30 TFF仅检测早期失稳，不提供最终统计Nu/Re、文献误差、长期稳定或网格独立性结论。
无效场检测沿用0.5 TFF统计采样：记录检测时刻及最后有限历史，不将其误称为严格首发时刻。
''')
assert MASTER.read_bytes() == raw
print(json.dumps({'local':str(OUT),'remote':REMOTE,'parameters':parameters},indent=2))
