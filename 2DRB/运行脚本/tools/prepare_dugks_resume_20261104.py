"""Prepare unchanged RB numerics: Ra=1e8 checkpoint resume, then fresh Ra=1e9."""
from pathlib import Path
import ast,copy,difflib,hashlib,json,re,tarfile
ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'运行脚本/run_cases/DUGKS_P100_node05_20261030'
OUT=ROOT/'运行脚本/run_cases/DUGKS_RB_resume_P100_node05_20261104'
REMOTE='/data2/XLLi/DUGKS/DUGKS_RB_resume_P100_node05_20261104'
OLDREMOTE='/data2/XLLi/DUGKS/DUGKS_P100_node05_20261030'
MASTER=ROOT/'DUGKS/2DRBOpenaccDUGKS.F90'
raw=MASTER.read_bytes(); base=raw.decode('utf-8-sig').replace('\r\n','\n')
sha=lambda b:hashlib.sha256(b).hexdigest()
assert not OUT.exists(),'Never overwrite a prepared campaign'
def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text,encoding='utf-8',newline='\n')
def once(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new)
tree=ast.parse((ROOT/'运行脚本/tools/prepare_dugks_p100_20261030.py').read_text(encoding='utf-8'))
driver=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='DRIVER' for t in n.targets))
driver=driver.replace('integer :: k, unit','integer :: k, unit, startStep')
driver=driver.replace('    do k=1,5','    startStep=itc\n    do k=1,5').replace('        itc=k','        itc=startStep+k')
cfg=json.loads((OLD/'campaign.json').read_text(encoding='utf-8'))
cases=[]
for exponent,n,resume in [(8,512,True),(9,1024,False)]:
    oldcase=next(c for c in cfg['cases'] if c['case_type']=='rb' and c['Ra']==10**exponent)
    name='RB_Ra1e%d_N%d_%sT1000'%(exponent,n,'restart100_' if resume else '')
    text=once(base,'nx=256, ny=256','nx=%d, ny=%d'%(n,n))
    for old,new in [('Rayleigh=1.0d6','Rayleigh=1.0d%d'%exponent),
                    ('flowCode=1, itc_max=20000000','flowCode=1, itc_max=ceiling(1000.0d0*timeUnit/dt)'),
                    ('\n#define steadyFlow\n','\n!#define steadyFlow\n'),
                    ('\n!#define unsteadyFlow\n','\n#define unsteadyFlow\n'),
                    ('\n#define SideHeatedCell\n','\n!#define SideHeatedCell\n'),
                    ('\n!#define RayleighBenardCell\n','\n#define RayleighBenardCell\n'),
                    ('Prandtl=0.71d0','Prandtl=0.7d0')]: text=once(text,old,new)
    oldtext=(OLD/oldcase['name']/'source/solver.F90').read_text(encoding='utf-8')
    assert text==oldtext,'Latest master numerics differ from cancelled case; checkpoint compatibility needs review'
    if resume: text=once(text,'parameter :: loadInitField=0','parameter :: loadInitField=1')
    smoke,count=re.subn(r'(?ms)^program main\s*$.*?^end program main\s*$',driver,text)
    assert count==1
    case=copy.deepcopy(oldcase)
    case.update(name=name,restart=resume,remote_dir=REMOTE+'/'+name,source_sha256=sha(text.encode()),
                smoke_sha256=sha(smoke.encode()),previous_case=oldcase['remote_dir'],
                previous_cancelled_job='7469.master' if resume else '7470.master')
    if resume: case.update(restart_tff=100.0000001351535,restart_step=25082775,
                          expected_seed_sha256='7f1cafbab071ced0a741ed12f29cd854a609271914c65385906274dcd615d249')
    folder=OUT/name
    write(folder/'source/solver.F90',text); write(folder/'source/smoke.F90',smoke)
    write(folder/'source/master_to_case.diff',''.join(difflib.unified_diff(base.splitlines(True),text.splitlines(True),fromfile='latest_local_DUGKS_master',tofile=name)))
    write(folder/'source/old_case_to_new.diff',''.join(difflib.unified_diff(oldtext.splitlines(True),text.splitlines(True),fromfile=oldcase['name'],tofile=name)))
    write(folder/'case.json',json.dumps(case,indent=2)+'\n')
    (folder/'results').mkdir()
    pbs=(OLD/oldcase['name']/'run.pbs').read_text(encoding='utf-8')
    pbs=pbs.replace(OLDREMOTE,REMOTE).replace(oldcase['name'],name).replace(oldcase['source_sha256'],case['source_sha256'])
    pbs=pbs.replace('export PYTHONHOME="$CAMPAIGN/runtime/python36"','export PYTHONHOME="'+OLDREMOTE+'/runtime/python36"')
    pbs=pbs.replace('#PBS -N DG_R'+str(exponent),'#PBS -N DG_R'+str(exponent)+('_RS100' if resume else '_NEW'))
    if resume:
        pbs=pbs.replace('nvfortran --version > compiler_version.txt 2>&1','sha256sum -c input_restart.sha256\nnvfortran --version > compiler_version.txt 2>&1')
    write(folder/'run.pbs',pbs)
    cases.append(case)
newcfg=dict(date='2026-11-04',node='node05',remote_dir=REMOTE,master=str(MASTER),master_sha256=sha(raw),
    normalized_master_sha256=sha(base.encode()),cases=cases,
    authorization='Resume the two unfinished DUGKS RB cases on node05. Do not rerun the five completed cases.',
    target_tff=1000,statistics_windows=[[500,750],[750,1000]],formal_dependencies='afterany: two independent Ra cases serial',
    previous_campaign=OLDREMOTE,numerics_unchanged=True)
write(OUT/'campaign.json',json.dumps(newcfg,indent=2,ensure_ascii=False)+'\n')
write(OUT/'source_master/2DRBOpenaccDUGKS.F90',base)
(OUT/'source_master/2DRBOpenaccDUGKS.original').write_bytes(raw)
write(OUT/'assess.py',(ROOT/'运行脚本/tools/dugks_campaign_assess.py').read_text(encoding='utf-8'))
preflight=(OLD/'preflight.py').read_text(encoding='utf-8')
preflight=preflight.replace("import hashlib,json,os,subprocess,sys","import hashlib,json,os,shutil,subprocess,sys")
preflight=preflight.replace("        env=os.environ.copy()", "        if case.get('restart'):\n            seed=json.loads((directory/'results/restart_preparation.json').read_text())\n            for name in [seed['checkpoint_name'],'reloadFile2DOpenaccDUGKS-latest.meta','Convergence_2DOpenaccDUGKS.dat']:\n                shutil.copy2(str(directory/'results'/name),str(work/name))\n        env=os.environ.copy()")
write(OUT/'preflight.py',preflight)
write(OUT/'gpu_gate.py',(OLD/'gpu_gate.py').read_text(encoding='utf-8'))
pbs=(OLD/'preflight/run.pbs').read_text(encoding='utf-8').replace(OLDREMOTE,REMOTE)
pbs=pbs.replace('#PBS -N DUGKS_GATE','#PBS -N DG_RESTART_GATE')
pbs=pbs.replace('export NVHPC_ROOT=', 'export PYTHONHOME="'+OLDREMOTE+'/runtime/python36"\nexport NVHPC_ROOT=')
pbs=pbs.replace('export PATH="$NVHPC_ROOT/compilers/bin:$PATH"','export PATH="$PYTHONHOME/bin:$NVHPC_ROOT/compilers/bin:$PATH"')
pbs=pbs.replace('export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:','export LD_LIBRARY_PATH="$PYTHONHOME/lib64:$NVHPC_ROOT/compilers/lib:')
write(OUT/'preflight/run.pbs',pbs)
write(OUT/'initial_node_check.json',(ROOT/'tmp/dugks_resume_nodecheck_20261104.json').read_text(encoding='utf-8-sig'))
write(OUT/'README.md','''# DUGKS RB 续算与重新提交（2026-11-04）

node05串行：Ra=1e8，512²，从100tff检查点续算至总计1000tff；Ra=1e9，1024²，从零开始至1000tff。
已完成的三个侧壁工况和两个低Ra的RB工况不重算。
所有数值参数与旧工况一致；唯一源差异为Ra=1e8的loadInitField=1。母版不改。
原始目录、取消记录与备份保持；续算使用新目录。Nu/Re历史由检查点缓存重建，原107tff未提交尾部不进入新轨迹。
重启前验证网格、宏、20项物理/输出参数、完整二进制状态、历史边界、签名与SHA-256。
早期场快照和Tecplot继续留在原目录，新轨迹恢复编号后接续输出，不伪造或覆盖旧场文件。
总时长1000tff不是额外1000tff；统计窗口仍为500–750、750–1000tff。
先在PBS/P100中编译两个精确工况并执行完整45步主机/GPU对照，再提交正式任务。
编译、短程一致性、完成时长、统计收敛、文献误差和网格无关性分别判断。
''')
archive=ROOT/'tmp/dugks_resume_20261104.tar.gz'
with tarfile.open(archive,'w:gz') as tar:
    for path in OUT.rglob('*'): tar.add(path,arcname=str(path.relative_to(OUT)),recursive=False)
assert MASTER.read_bytes()==raw
print(json.dumps(dict(directory=str(OUT),remote=REMOTE,archive=str(archive),cases=[dict(name=c['name'],restart=c['restart'],source_sha256=c['source_sha256']) for c in cases]),indent=2))
