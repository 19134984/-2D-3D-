"""Prepare the requested FPC Ra=1e10, 2048-square RB case from the latest master."""
from pathlib import Path
import difflib
import hashlib
import json
import math
import re
import shutil
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'运行脚本/run_cases/FPC_RB_Ra1e10_N2048_T1000_20261101'
REMOTE='/data2/XLLi/FPC/FPC_RB_Ra1e10_N2048_T1000_20261101'
NAME='RB_Ra1e10_N2048'
NODE='node07'
OLD=ROOT/'运行脚本/run_cases/FPC_P100_node07_20261029'
MASTER=ROOT/'FPC/2DRBOpenaccFPC.F90'
raw=MASTER.read_bytes()
baseline=raw.decode('utf-8-sig').replace('\r\n','\n')
assert not (OUT/'campaign.json').exists(), 'Never overwrite a prepared campaign'
OUT.mkdir(parents=True,exist_ok=True)
case_dir=OUT/NAME
case_dir.mkdir()
(case_dir/'results').mkdir()

def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text,encoding='utf-8',newline='\n')

def sha(data):
    return hashlib.sha256(data).hexdigest()

def once(text,a,b):
    assert text.count(a)==1,a
    return text.replace(a,b)

text=baseline
for a,b in [('\n#define steadyFlow\n','\n!#define steadyFlow\n'),
            ('\n!#define unsteadyFlow\n','\n#define unsteadyFlow\n'),
            ('\n#define SideHeatedCell\n','\n!#define SideHeatedCell\n'),
            ('\n!#define RayleighBenardCell\n','\n#define RayleighBenardCell\n'),
            ('nx=256, ny=256','nx=2048, ny=2048'),
            ('Rayleigh=1.0d6','Rayleigh=1.0d10'),
            ('Prandtl=0.71d0','Prandtl=0.7d0')]:
    text=once(text,a,b)
write(case_dir/'source/solver.F90',text)
write(case_dir/'source/master_to_case.diff',''.join(difflib.unified_diff(
    baseline.splitlines(True),text.splitlines(True),fromfile='latest_FPC_master',tofile=NAME)))
old_smoke=(OLD/'RB_Ra1e9_N1024/source/smoke.F90').read_text(encoding='utf-8')
driver=re.search(r'^    program main\b.*?^    end program main\b',old_smoke,re.M|re.S).group(0)
smoke,count=re.subn(r'^    program main\b.*?^    end program main\b',lambda m:driver,text,flags=re.M|re.S)
assert count==1
write(case_dir/'source/smoke.F90',smoke)

# Prove that the previous validated RB solver differs only in Ra/grid and its inactive steady cap.
previous=(OLD/'RB_Ra1e9_N1024/source/solver.F90').read_text(encoding='utf-8')
expected=previous.replace('nx=1024, ny=1024','nx=2048, ny=2048').replace('Rayleigh=1.0d9','Rayleigh=1.0d10')
expected=once(expected,'itc_max=nint(1000.0d0*timeUnit)','itc_max=20000000')
assert text==expected, 'Unexpected change to FPC physics/infrastructure'
write(case_dir/'source/previous_RB9_to_RB10.diff',''.join(difflib.unified_diff(
    previous.splitlines(True),text.splitlines(True),fromfile='validated_RB_Ra1e9_N1024',tofile=NAME)))

nu=.01*2048*math.sqrt(.7/1e10)
case=dict(name=NAME,case_type='rb',nx=2048,ny=2048,Ra=1e10,Pr=.7,
    Mach=math.sqrt(3)*.01,Uff=.01,Thot=2,Tcold=1,viscosity=nu,diffusivity=nu/.7,
    tauf=.5+3*nu,taug=.5+3*nu/.7,target_tff=1000,step_cap=204800000,
    restart=False,node=NODE,epsU=1e-7,epsT=1e-7,finite_check_interval_steps=2000,
    sample_interval_tff=.5,checkpoint_interval_tff=100,snapshot_enabled=True,
    tecplot_enabled=True,source_sha256=sha(text.encode()),smoke_sha256=sha(smoke.encode()),
    remote_dir=REMOTE+'/'+NAME,reference_Nu=93.88,reference_Re=61797,
    initial_condition='zero physical velocity; conductive profile plus deterministic 1e-3 relative temperature perturbation',
    walltime='720:00:00')
write(case_dir/'case.json',json.dumps(case,indent=2)+'\n')
old_case=json.loads((OLD/'RB_Ra1e9_N1024/case.json').read_text())
pbs=(OLD/'RB_Ra1e9_N1024/run.pbs').read_text(encoding='utf-8')
pbs=pbs.replace(old_case['remote_dir'],case['remote_dir'])
pbs=pbs.replace('/data2/XLLi/FPC/FPC_P100_node07_20261029',REMOTE)
pbs=once(pbs,'#PBS -N FPC_R9','#PBS -N FPC_R10')
pbs=once(pbs,'#PBS -l walltime=240:00:00','#PBS -l walltime=720:00:00')
pbs=once(pbs,old_case['source_sha256'],case['source_sha256'])
write(case_dir/'run.pbs',pbs)
shutil.copyfile(OLD/'assess.py',OUT/'assess.py')
gate=(OLD/'gpu_gate.py').read_text(encoding='utf-8')
write(OUT/'gpu_gate.py',once(gate,'free >= 2048','free >= 4096'))
(OUT/'source_master').mkdir()
(OUT/'source_master/2DRBOpenaccFPC.original').write_bytes(raw)
write(OUT/'source_master/2DRBOpenaccFPC.F90',baseline)
protected=['FPC/2DRBOpenaccFPC.F90','均匀网格/2DRBOpenmp.F90','均匀网格/2DRBOpenacc.F90']
cfg=dict(created='2026-11-01',remote_dir=REMOTE,node=NODE,latest_master=str(MASTER),
    master_sha256=sha(raw),case=case,cases=[case],
    protected_sources={p:sha((ROOT/p).read_bytes()) for p in protected},
    node_selection='2026-11-01 live check: node07 has no PBS jobs or GPU compute processes; node05 is occupied',
    intentional_changes='Only case macros, Ra, nx/ny and Pr; no change to collision, forcing, boundaries, initialization or post-processing',
    parent_campaign='FPC_P100_node07_20261029',
    statistical_windows_tff=[[500,750],[750,1000]],final_statistics_window_tff=[750,1000],
    Re_statistic='sqrt(mean(Re(t)^2)); independent assess.py overrides the built-in arithmetic-Re summary',
    walltime='720:00:00',estimated_lattice_updates=204800000*2048**2,
    relative_updates_vs_Ra9_N1024=8,
    reference='Zhang, Zhou and Sun (2017), JFM 814, Table 1; Pr=0.7, Ra=1e10; user requests 2048 rather than 2049')
write(OUT/'campaign.json',json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
write(OUT/'README.md','''# FPC，Ra=10^10，2048² RB 非稳态

2026-11-01 用户要求新增算例；用提交时最新本地 FPC 主程序生成独立参数化副本，不覆盖原七组结果。
Ra=1e10，nx=ny=2048，Pr=0.7，Uff=0.01，Ma=sqrt(3)*0.01，Thot/Tcold=2/1。
下热上冷，左右绝热，四壁无滑移；从零物理速度和原线性温度加确定性1e-3扰动启动。
不加载Ra=1e9检查点；FPC碰撞、力项、边界、后处理和输出设置保持原样。
目标1000TFF，Nu/Re及宏观快照0.5TFF；检查点与Tecplot100TFF；每2000步检查场有限性。
统计比较500–750与750–1000TFF，最终取后一窗；Re使用sqrt(mean(Re(t)^2))。
参考JFM 2017表1的Pr=0.7行：Nu=93.88，Re=61797。原文2049²按用户要求使用2048²。
实时资源检查后选node07；720h是PBS运行时限，1000TFF仍是求解目标。
每个算例独立保存源码、差异、主源/参数化源/编译源哈希；先实际2048²的40步CPU/GPU对照，再进入正式求解。
该网格和同样1000TFF的更新量为前一1024²工况的8倍，耗时需以实际运行测量。
编译、短程对照、完成目标、有限性、半窗检查及文献精度分别验收，不以短程通过代替长期结果。
''')

build=Path('C:/Users/User/.codex/visualizations/2026/10/28/01a180ca-1d4f-7000-8edc-cffeb3392d7c/fpc_ra10_2048_syntax')
build.mkdir(parents=True,exist_ok=True)
for name in ['solver','smoke']:
    target=build/(name+'.F90')
    shutil.copyfile(case_dir/('source/'+name+'.F90'),target)
    result=subprocess.run(['C:/msys64/ucrt64/bin/gfortran.exe','-cpp','-fopenacc',
        '-ffree-line-length-none','-fsyntax-only',str(target)],cwd=build,capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    write(OUT/('evidence/'+name+'_syntax.log'),result.stdout+result.stderr)
for path,expected_sha in cfg['protected_sources'].items():
    assert sha((ROOT/path).read_bytes())==expected_sha
write(OUT/'evidence/local_preflight.json',json.dumps(dict(passed=True,
    production_and_smoke_syntax=True,master_unchanged=True,
    equivalent_to_previous_RB_except_Ra_grid=True,source_sha256=case['source_sha256']),indent=2)+'\n')
with tarfile.open(OUT/'campaign.tar','w') as archive:
    for path in sorted(OUT.rglob('*')):
        relative=path.relative_to(OUT)
        if path.is_file() and relative.parts[0]!='evidence' and path.name!='campaign.tar':
            archive.add(path,arcname=OUT.name+'/'+relative.as_posix())
print(json.dumps(dict(out=str(OUT),remote=REMOTE,case=case,preflight_passed=True),indent=2))
