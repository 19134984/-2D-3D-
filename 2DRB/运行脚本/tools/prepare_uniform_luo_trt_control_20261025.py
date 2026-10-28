"""Restore TRT in a copy of the failed MRT run for a controlled Ra9 test."""
from pathlib import Path
import difflib
import hashlib
import json
import tarfile

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT/'运行脚本/run_cases/Uniform_LuoMRT_Ra9_Ra10_T1000_20261025'
OUT = ROOT/'运行脚本/run_cases/Uniform_LuoTRT_Ra9_control_T1000_20261025'
REMOTE = '/data2/XLLi/LBMCDE/Uniform_LuoTRT_Ra9_control_T1000_20261025'
assert not OUT.exists(), 'Refuse to overwrite an existing control'
OUT.mkdir(parents=True)
source_path = PARENT/'Ra1e9/source/solver.F90'
raw = source_path.read_bytes()
assert hashlib.sha256(raw).hexdigest() == '029a1134e266c9a27f3e979f490074a42ab0d0aec1a468669f936f84924cf7ca'
original = raw.decode('utf-8')


def once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old,new)


def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text,encoding='utf-8',newline='\n')


text = once(original,'!!!    均匀网格独立测试：流场 Luo D2Q9-MRT；温度场 Luo D2Q5-TRT',
    '!!!    同源对照：流场恢复原始 D2Q9-TRT；温度场仍为 Luo D2Q5-TRT')
start = text.index('        ! Contrino/Luo 2014 §5.3：')
end = text.index('\n',text.index('        real(kind=8), parameter :: bulkViscosity=',start))
text = text[:start]+('        ! 同源 TRT 对照：能量矩、四阶矩恢复为 Snu；Sq 保留原始魔法参数关系。\n'
    '        real(kind=8), parameter :: bulkViscosity=(tauf-0.5d0)*(1.0d0-chi_b)*cs2')+text[end:]
text = once(text,'            s(1) = flowEnergyRate','            s(1) = Snu')
text = once(text,'            s(2) = flowEpsilonRate','            s(2) = Snu')
text = once(text,'    ! 本诊断副本只允许无 CDE 修正的论文 MRT 配置。',
    '    ! 本 TRT 对照与失败 MRT 一样，只允许 chi_nu=chi_b=0。')
text = once(text,'error stop "MRT test requires chi=0"','error stop "TRT control requires chi=0"')
text = once(text,'Flow collision = Luo MRT; sq follows paper Eq. (7b)',
    'Flow collision = Luo TRT control; se=sepsilon=snu; original magic')
text = once(text,"'flow s_e / s_epsilon =',flowEnergyRate,flowEpsilonRate",
    "'flow s_e / s_epsilon =',Snu,Snu")
text = once(text,'    ! 保留 3/16 的参数关系不等于证明 MRT 具有全部 TRT 边界性质。',
    '    ! 原 TRT：所有非守恒偶矩率为 Snu，非守恒奇矩率为 Sq。')
assert 'flowEnergyRate' not in text and 'flowEpsilonRate' not in text
assert 'nx=1025, ny=1025' in text and 'Rayleigh=1.0d9' in text
assert 'unsteadyRunDuration=1000.0d0' in text
write(OUT/'Ra1e9/source/solver.F90',text)
write(OUT/'Ra1e9/source/MRT_to_TRT.diff',''.join(difflib.unified_diff(
    original.splitlines(True),text.splitlines(True),fromfile='failed_MRT',tofile='same_source_TRT')))
master_raw=(PARENT/'source_master/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90').read_bytes()
write(OUT/'Ra1e9/source/master_to_control.diff',''.join(difflib.unified_diff(
    master_raw.decode('utf-8-sig').replace('\r\n','\n').splitlines(True),text.splitlines(True),
    fromfile='unchanged_Xs_master',tofile='TRT_control')))
(OUT/'Ra1e9/results').mkdir()
(OUT/'source_parent').mkdir()
(OUT/'source_parent/failed_MRT_solver.F90').write_bytes(raw)
cfg=json.loads((PARENT/'campaign.json').read_text(encoding='utf-8'))
case=cfg['cases']['Ra1e9']
case.update(s_e=case['s_nu'],s_epsilon=case['s_nu'],bulk_viscosity=case['viscosity'],
    source_sha256=hashlib.sha256(text.encode('utf-8')).hexdigest(),flow='original TRT')
cfg.update(remote_dir=REMOTE,order=['Ra1e9'],cases={'Ra1e9':case},
    parent_failed_MRT_source_sha256=hashlib.sha256(raw).hexdigest(),
    parent_frozen_source=str(source_path),
    comparison='Only restore se=sepsilon=snu and its bulk diagnostic; same grid, initial state, thermal solver, forcing and output schedule.')
write(OUT/'campaign.json',json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
for name in ['audit_flow_history.py','protected_sources.json']:
    write(OUT/name,(PARENT/name).read_text(encoding='utf-8'))
validator=(PARENT/'validate_case.py').read_text(encoding='utf-8')
validator=once(validator,"if 'Luo MRT' not in settings", "if 'Luo TRT control' not in settings")
validator=once(validator,"'MRT identity missing from runtime settings'","'TRT identity missing from runtime settings'")
write(OUT/'validate_case.py',validator)
runner=(PARENT/'run_case.sh').read_text(encoding='utf-8')
runner=once(runner,'ROOT=/data2/XLLi/LBMCDE/Uniform_LuoMRT_Ra9_Ra10_T1000_20261025',
            'ROOT='+REMOTE)
runner=once(runner,'case "$name" in Ra1e9|Ra1e10) ;; *) exit 2;; esac',
    'case "$name" in Ra1e9) ;; *) exit 2;; esac')
start=runner.index('if [ "$name" = Ra1e9 ]; then')
end=runner.index('\nset +e\n/usr/bin/time',start)
runner=runner[:start]+runner[end:]
write(OUT/'run_case.sh',runner)
write(OUT/'Ra1e9/run.pbs',f'''#!/bin/bash
#PBS -N LuoTRT_Ra9_Ctl
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=96:00:00
#PBS -o {REMOTE}/Ra1e9/results/pbs.stdout
#PBS -e {REMOTE}/Ra1e9/results/pbs.stderr
set -euo pipefail
exec bash {REMOTE}/run_case.sh Ra1e9
''')
write(OUT/'submit.sh',f'''#!/bin/bash
set -euo pipefail
cd {REMOTE}
test "$(hostname -s)" = master
sha256sum -c upload.sha256
bash -n run_case.sh
bash -n Ra1e9/run.pbs
mkdir .submission.lock
trap 'rmdir .submission.lock' EXIT
test ! -s Ra1e9/job.id
job=$(qsub {REMOTE}/Ra1e9/run.pbs)
printf '%s\n' "$job" | tee Ra1e9/job.id
qstat -f "$job" | tee submission_qstat.txt
''')
write(OUT/'README.md','''# 同源 TRT 九次方对照

用户要求：怀疑测试代码有问题，将失败MRT副本中的1.2、1.3恢复原TRT关系，再测Ra=10^9。
使用冻结的失败Ra9源码作唯一母本；恢复s(1)=s(2)=Snu，体黏性诊断同步恢复原公式。
不另换求解器，不改变初始化、网格1025²、Pr=.7、Ma=.1、chi_nu=chi_b=0、D2Q5温度场、力项和采样。
Sq仍按原始魔法参数，(1/Snu-.5)(1/Sq-.5)=3/16。没有读取旧场，从零计算。
目标仍1000 TFF；先重点检查能否越过MRT的15 TFF失败点，通过短时不代表长期稳定。
只有一组PBS作业，node07，96小时；非数保护、快照和100 TFF重启输出沿用母本。
失败MRT记录、多块源码及正式Xs/均匀母版均不修改。
源码差异在Ra1e9/source/MRT_to_TRT.diff；来源SHA256见campaign.json。
运行前校验源哈希，结束后检查完整历史、Nu/Re RMS统计、耗散/剖面与检查点账本。
''')
for rel,expected in json.loads((OUT/'protected_sources.json').read_text(encoding='utf-8')).items():
    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==expected,rel
files=[p for p in OUT.rglob('*') if p.is_file()]
manifest=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(OUT).as_posix()+'\n' for p in sorted(files))
write(OUT/'upload.sha256',manifest)
with tarfile.open(OUT/'campaign.tar.gz','w:gz') as archive:
    for p in files+[OUT/'upload.sha256']:
        archive.add(p,arcname=p.relative_to(OUT).as_posix(),recursive=False)
    info=tarfile.TarInfo('Ra1e9/results')
    info.type=tarfile.DIRTYPE
    info.mode=0o755
    archive.addfile(info)
print(json.dumps(dict(remote=REMOTE,source_sha256=case['source_sha256'],
    archive_sha256=hashlib.sha256((OUT/'campaign.tar.gz').read_bytes()).hexdigest()),indent=2))
