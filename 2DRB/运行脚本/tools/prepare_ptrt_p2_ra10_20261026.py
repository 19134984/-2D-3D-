"""Prepare only the user-authorized P2 Ra=1e10 run; preserve all earlier jobs."""
from pathlib import Path
import difflib
import hashlib
import json
import tarfile
import time

ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'运行脚本/run_cases/LuoPTRT_Ra9_Ra10_20261025_r2'
OUT=ROOT/'运行脚本/run_cases/LuoPTRT_P2_Ra1e10_T1000_20261026'
REMOTE='/data2/XLLi/PTRT/P2_Ra1e10_T1000_20261026'
MASTER=ROOT/'PTRT/2DRBOpenaccLuoPTRT.F90'
NAME='P2_Ra1e10'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text,encoding='utf-8',newline='\n')


def change(text,old,new):
    assert text.count(old)==1,(old,text.count(old))
    return text.replace(old,new)


assert not OUT.exists(),'Refuse to replace an earlier attempt'
previous=json.loads((OLD/'campaign.json').read_text(encoding='utf-8'))
assert sha(MASTER)==previous['master_sha256'],'Revalidate a changed algorithm before submitting'
original=MASTER.read_text(encoding='utf-8')
source=change(original,'ptrtSqStrategy=1','ptrtSqStrategy=2')
source=change(source,'nx=1024, ny=1024','nx=2048, ny=2048')
source=change(source,'Rayleigh=1.0d9','Rayleigh=1.0d10')
assert 'loadInitField=0' in source and 'unsteadyRunDuration=1000.0d0' in source
assert hashlib.sha256(source.encode('utf-8')).hexdigest()==previous['cases'][NAME]['source_sha256']
write(OUT/NAME/'source/solver.F90',source)
write(OUT/NAME/'source/master_to_case.diff',''.join(difflib.unified_diff(
    original.splitlines(True),source.splitlines(True),fromfile='PTRT/2DRBOpenaccLuoPTRT.F90',
    tofile=NAME+'/source/solver.F90')))
(OUT/'source_master').mkdir()
(OUT/'source_master'/MASTER.name).write_bytes(MASTER.read_bytes())
(OUT/NAME/'results').mkdir()
cfg=dict(remote_dir=REMOTE,node='node07',order=[NAME],master=str(MASTER),master_sha256=sha(MASTER),
    cases={NAME:previous['cases'][NAME]},release_policy=previous['release_policy'],
    launch_authorization='User explicitly requested P2 Ra1e10 after finite 1000-tff Ra1e9; do not run P1.',
    launch_policy='Require completed finite P2 Ra1e9; its <1% statistical convergence is not a launch condition.',
    from_scratch=True,dependency=None,previous_results=previous['remote_dir'])
write(OUT/'campaign.json',json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
protected={p.relative_to(ROOT).as_posix():sha(p) for p in [MASTER,
    ROOT/'Xs/D2Q9TRT_D2Q5LuoTRT/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90',
    ROOT/'多块网格/2DRBOpenaccMultiblock.F90',*ROOT.glob('均匀网格/*.F90')]}
write(OUT/'protected_sources.json',json.dumps(protected,ensure_ascii=False,indent=2)+'\n')
for src,dst in [('运行脚本/tools/assess_luo_ptrt.py','assess_case.py'),('后处理/audit_flow_history.py','audit_flow_history.py')]:
    write(OUT/dst,(ROOT/src).read_text(encoding='utf-8-sig'))
write(OUT/'prior_validation.json',json.dumps(dict(source_unchanged=True,
    exact_source_syntax= json.loads((OLD/'full_source_host_validation.json').read_text())['syntax'][NAME],
    algebra_and_small_cases=json.loads((OLD/'local_validation.json').read_text())['pass'],
    GPU_validation_exit_code=0,checkpoint_check=json.loads((OLD/'result_check_20261026/result_check_20261026.json').read_text()),
    reference_dir=str(OLD)),indent=2)+'\n')
write(OUT/'preflight.py','''from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
cfg=json.loads((root/'campaign.json').read_text())
old=Path(cfg['previous_results'])
case=old/'P2_Ra1e9/results'
a=json.loads((case/'assessment.json').read_text())
assert int((case/'run.status').read_text())==0
assert a['structural_pass'] and a['target_reached'] and a['finite_histories']
assert abs(a['last_checkpoint_tff']-1000)<1e-3
gpu=json.loads((old/'GPU_validation/results/gpu_comparison.json').read_text())
assert gpu['P2']['pass']
source=root/'P2_Ra1e10/source/solver.F90'
assert hashlib.sha256(source.read_bytes()).hexdigest()==cfg['cases']['P2_Ra1e10']['source_sha256']
print(json.dumps(dict(authorization=cfg['launch_authorization'],Ra9_completed=True,
    Ra9_statistics_converged=a['statistics_converged'],statistical_gate_required=False,
    target_from_scratch=True,P1_submitted=False),indent=2))
''')
runner='''#!/bin/bash
set -euo pipefail
ROOT=__REMOTE__
NAME=P2_Ra1e10
RESULT="$ROOT/$NAME/results"
test "$(hostname -s)" = node07
cd "$RESULT"
test ! -e run.started
date -Is > run.started
trap 'rc=$?; printf "%s\\n" "$rc" > "$RESULT/run.status"; date -Is > "$RESULT/run.finished"' EXIT
printf 'job=%s\\nhost=%s\\ncase=%s\\n' "${PBS_JOBID:-unknown}" "$(hostname)" "$NAME" > job_identity.txt
(cd "$ROOT" && sha256sum -c upload.sha256) > input_hash_verification.log
/usr/bin/python3 "$ROOT/preflight.py" > launch_authorization_check.json
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia
export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}
nvfortran --version > compiler_version.txt 2>&1
nvidia-smi > nvidia_smi_before.log 2>&1 || true
for attempt in $(seq 1 120); do
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits | head -1)
    test "$util" -lt 10 && break
    if test "$attempt" -eq 120; then echo 'GPU remained busy; solver not launched' >&2; exit 75; fi
    sleep 30
done
sha256sum "$ROOT/$NAME/source/solver.F90" > runtime_source.sha256
FLAGS=(-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend)
printf '%s\\n' "${FLAGS[*]}" > build_flags.txt
set +e
nvfortran "${FLAGS[@]}" "$ROOT/$NAME/source/solver.F90" -o solver.exe > compile.log 2>&1
compile_rc=$?
printf '%s\\n' "$compile_rc" > compile.status
set -e
test "$compile_rc" -eq 0 || exit "$compile_rc"
set +e
/usr/bin/time -p -o timing.txt ./solver.exe > solver.stdout 2>&1
solver_rc=$?
printf '%s\\n' "$solver_rc" > solver.status
/usr/bin/python3 "$ROOT/assess_case.py" "$NAME" > assessment.log 2>&1
assessment_rc=$?
set -e
nvidia-smi > nvidia_smi_after.log 2>&1 || true
test "$solver_rc" -eq 0 || exit "$solver_rc"
exit "$assessment_rc"
'''.replace('__REMOTE__',REMOTE)
write(OUT/'run_case.sh',runner)
write(OUT/'run.pbs',f'''#!/bin/bash
#PBS -N PTRT_P2_Ra10
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=96:00:00
#PBS -o {REMOTE}/{NAME}/results/pbs.stdout
#PBS -e {REMOTE}/{NAME}/results/pbs.stderr
set -euo pipefail
exec bash {REMOTE}/run_case.sh
''')
write(OUT/'submit.sh',f'''#!/bin/bash
set -euo pipefail
cd {REMOTE}
test ! -s submission.txt
sha256sum -c upload.sha256 > upload_verification.log
/usr/bin/python3 preflight.py > submission_preflight.json
qsub run.pbs > submission.txt
date -Is > submission_time.txt
cat submission.txt
''')
write(OUT/'provenance.md','''# P2 Ra=10^10 独立测试

用户于 2026-10-26 明确要求单独计算 P2 的 Ra=10^10，不再提交 P1。
P2 的 Ra=10^9 已完成 1000 t_ff 且最终 f/g 检查点有限；统计窗口偏差仍大于 1%。
此次明确取消“Ra9 统计偏差必须先小于 1%”这一启动条件，保留所有数值有限性检查。
稳定运行时长、统计收敛与文献符合程度分别报告。

从提交时最新本地 PTRT/2DRBOpenaccLuoPTRT.F90 复制，仅参数化为 P2、2048²、Ra=1e10。
Pr=.7、Ma=.1、chi_nu=chi_b=0、Sq=1.6、原 Luo D2Q5 温度与边界不变。
从零启动，目标 1000 t_ff；此任务不读取 Ra9 检查点，不依赖 P1。
数值源与此前已通过语法检查的 P2 Ra10 源逐字一致；PBS 使用 NVHPC 24.3 重新编译。
node07 实时 GPU 显存为空且无 PBS 作业；node05 另有进程使用约 448 MiB。
旧算例、旧提交状态及其源快照保留不动。本次只提交一个正式作业。
''')
paths=sorted(p for p in OUT.rglob('*') if p.is_file())
write(OUT/'upload.sha256',''.join(sha(p)+'  '+p.relative_to(OUT).as_posix()+'\n' for p in paths))
with tarfile.open(OUT/'campaign.tar.gz','w:gz') as tar:
    for p in paths+[OUT/'upload.sha256']: tar.add(p,arcname=p.relative_to(OUT).as_posix(),recursive=False)
    info=tarfile.TarInfo(NAME+'/results');info.type=tarfile.DIRTYPE;info.mode=0o755;info.mtime=int(time.time());tar.addfile(info)
assert all(sha(ROOT/p)==h for p,h in protected.items())
print(json.dumps(dict(local=str(OUT),remote=REMOTE,node='node07',master_sha256=sha(MASTER),
    source_sha256=sha(OUT/NAME/'source/solver.F90'),archive_sha256=sha(OUT/'campaign.tar.gz'),
    exact_previously_checked_source=True,protected_sources_unchanged=True),ensure_ascii=False,indent=2))
