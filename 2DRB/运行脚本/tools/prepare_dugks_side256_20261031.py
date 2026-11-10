"""Prepare three replacement steady cases; preserve the original seven-case campaign."""
from pathlib import Path
import ast
import copy
import difflib
import hashlib
import json
import math
import re
import tarfile

ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'运行脚本/run_cases/DUGKS_P100_node05_20261030'
OUT=ROOT/'运行脚本/run_cases/DUGKS_Side256_P100_node05_20261031'
REMOTE='/data2/XLLi/DUGKS/DUGKS_Side256_P100_node05_20261031'
OLDREMOTE='/data2/XLLi/DUGKS/DUGKS_P100_node05_20261030'
MASTER=ROOT/'DUGKS/2DRBOpenaccDUGKS.F90'
raw=MASTER.read_bytes()
base=raw.decode('utf-8-sig').replace('\r\n','\n')
sha=lambda b: hashlib.sha256(b).hexdigest()
assert not OUT.exists(), 'Never overwrite a prepared campaign'
def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text,encoding='utf-8',newline='\n')
def once(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new)
tree=ast.parse((ROOT/'运行脚本/tools/prepare_dugks_p100_20261030.py').read_text(encoding='utf-8'))
driver=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='DRIVER' for t in n.targets))
cfg=json.loads((OLD/'campaign.json').read_text(encoding='utf-8'))
cases=[]
for oldcase in cfg['cases'][:3]:
    case=copy.deepcopy(oldcase)
    name=oldcase['name'].replace('_N512','_N256')
    exponent=int(math.log10(oldcase['Ra']))
    text=once(base,'Rayleigh=1.0d6','Rayleigh=1.0d%d'%exponent)
    text=once(text,'flowCode=1, itc_max=20000000','flowCode=1, itc_max=ceiling(1000.0d0*timeUnit/dt)')
    oldtext=(OLD/oldcase['name']/'source/solver.F90').read_text(encoding='utf-8')
    assert text==once(oldtext,'nx=512, ny=512','nx=256, ny=256'), 'Unexpected parameter or method change'
    smoke,count=re.subn(r'(?ms)^program main\s*$.*?^end program main\s*$',driver,text)
    assert count==1
    case.update(name=name,nx=256,ny=256,remote_dir=REMOTE+'/'+name,
                source_sha256=sha(text.encode()),smoke_sha256=sha(smoke.encode()),
                step_cap=math.ceil(1000*256/case['Uff']/case['dt']),
                estimated_device_array_MiB=math.ceil((112*256*256+36*256)*8/2**20),
                walltime_hours=264,walltime='264:00:00',replaces_job=str(7464+len(cases))+'.master',
                replaced_case=oldcase['name'],restart=False)
    case.pop('preflight_estimate_1000_tff_hours',None)
    folder=OUT/name
    write(folder/'source/solver.F90',text)
    write(folder/'source/smoke.F90',smoke)
    write(folder/'source/master_to_case.diff',''.join(difflib.unified_diff(base.splitlines(True),text.splitlines(True),fromfile='latest_local_master',tofile=name)))
    write(folder/'source/old512_to_new256.diff',''.join(difflib.unified_diff(oldtext.splitlines(True),text.splitlines(True),fromfile=oldcase['name'],tofile=name)))
    write(folder/'case.json',json.dumps(case,indent=2)+'\n')
    (folder/'results').mkdir()
    pbs=(OLD/oldcase['name']/'run.pbs').read_text(encoding='utf-8')
    pbs=pbs.replace(OLDREMOTE,REMOTE).replace(oldcase['name'],name).replace(oldcase['source_sha256'],case['source_sha256'])
    pbs=pbs.replace('export PYTHONHOME="$CAMPAIGN/runtime/python36"','export PYTHONHOME="'+OLDREMOTE+'/runtime/python36"')
    pbs=pbs.replace('#PBS -N DG_S'+str(exponent),'#PBS -N DG256_S'+str(exponent))
    # Preserve the established scheduler ceiling; this is a cap, not a timing prediction.
    write(folder/'run.pbs',pbs)
    cases.append(case)
newcfg=dict(date='2026-10-31',node='node05',remote_dir=REMOTE,master=str(MASTER),
    master_sha256=sha(raw),normalized_master_sha256=sha(base.encode()),cases=cases,
    user_authorized_gpu_sharing=True,original_campaign=OLDREMOTE,
    authorization='Replace only three steady 512-square cases with 256-square cases. Keep all four existing RB jobs and their numerical settings.',
    changes=['nx=ny=256; derived physical scaling and step cap recomputed'],
    runtime_python=OLDREMOTE+'/runtime/python36',formal_dependencies='three steady cases serial, then the four unchanged existing RB cases')
write(OUT/'campaign.json',json.dumps(newcfg,indent=2,ensure_ascii=False)+'\n')
write(OUT/'source_master/2DRBOpenaccDUGKS.F90',base)
(OUT/'source_master/2DRBOpenaccDUGKS.original').write_bytes(raw)
write(OUT/'assess.py',(ROOT/'运行脚本/tools/dugks_campaign_assess.py').read_text(encoding='utf-8'))
write(OUT/'preflight.py',(OLD/'preflight.py').read_text(encoding='utf-8'))
write(OUT/'gpu_gate.py', '''"""Record current shared-GPU resources; sharing explicitly authorized."""
import datetime,json,subprocess
gpu=subprocess.check_output(['nvidia-smi','--query-gpu=memory.free,utilization.gpu','--format=csv,noheader,nounits'],universal_newlines=True).strip()
processes=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv'],universal_newlines=True)
free,util=map(int,gpu.splitlines()[0].split(','))
result=dict(time=datetime.datetime.now().isoformat(),sharing_authorized=True,free_MiB=free,utilization=util,processes=processes,passed=free>=1024)
print(json.dumps(result,indent=2))
if not result['passed']: raise SystemExit('Insufficient free memory for the 256-square case')
''')
pbs=(OLD/'preflight/run.pbs').read_text(encoding='utf-8').replace(OLDREMOTE,REMOTE)
pbs=pbs.replace('#PBS -N DUGKS_GATE','#PBS -N DG256_GATE')
pbs=pbs.replace('export NVHPC_ROOT=', 'export PYTHONHOME="'+OLDREMOTE+'/runtime/python36"\nexport NVHPC_ROOT=')
pbs=pbs.replace('export PATH="$NVHPC_ROOT/compilers/bin:$PATH"','export PATH="$PYTHONHOME/bin:$NVHPC_ROOT/compilers/bin:$PATH"')
pbs=pbs.replace('export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:', 'export LD_LIBRARY_PATH="$PYTHONHOME/lib64:$NVHPC_ROOT/compilers/lib:')
write(OUT/'preflight/run.pbs',pbs)
write(OUT/'README.md','''# DUGKS 侧壁稳态 256² 替代任务

Ra=1e6/1e7/1e8，Pr=0.71，单 GPU OpenACC，node05 共享运行，三个独立工况串行。
仅网格从512²改为256²；Mach=.01、CFL=.5、冷热壁±.5、收敛标准和输出间隔均保持。
从初始场重新计算，最多1000自由落体时间；连续三次errorU/errorT<1e-7停止。
保留旧512²目录；旧作业取消记录单列。既有四个RB作业源码、网格、输出路径不变，仅首个RB前驱改为新的侧壁末作业。
264小时是调度上限，不是预期耗时。短程编译/GPU主机一致性不代表已收敛或已通过文献验证。
''')
write(OUT/'initial_node_check.json',(ROOT/'tmp/dugks_replace256_initial.txt').read_text(encoding='utf-8-sig'))
archive=ROOT/'tmp/dugks_side256_20261031.tar.gz'
with tarfile.open(archive,'w:gz') as tar:
    for path in OUT.rglob('*'):
        tar.add(path,arcname=str(path.relative_to(OUT)),recursive=False)
assert MASTER.read_bytes()==raw
print(json.dumps(dict(local=str(OUT),remote=REMOTE,archive=str(archive),master_sha256=sha(raw),cases=[dict(name=c['name'],source_sha256=c['source_sha256'],step_cap=c['step_cap']) for c in cases]),indent=2))
