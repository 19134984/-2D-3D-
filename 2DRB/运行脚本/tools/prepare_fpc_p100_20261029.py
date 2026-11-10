"""Prepare the user's seven FPC cases without changing any solver master."""
from pathlib import Path
import difflib
import hashlib
import json
import math
import re
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '运行脚本/run_cases/FPC_P100_node07_20261029'
REMOTE = '/data2/XLLi/FPC/FPC_P100_node07_20261029'
MASTER = ROOT / 'FPC/2DRBOpenaccFPC.F90'
raw = MASTER.read_bytes()
baseline = raw.decode('utf-8-sig').replace('\r\n', '\n')
assert not (OUT / 'campaign.json').exists(), 'Do not replace an existing campaign'
OUT.mkdir(parents=True, exist_ok=True)

def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8', newline='\n')

def digest(data):
    return hashlib.sha256(data).hexdigest()

def once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)

driver = '''    program main
    use openacc
    use commondata
    implicit none
    integer(kind=8) :: clock0, clock1, rate
    integer :: k
    call initial()
    call enter_data_2d_openacc()
    call system_clock(clock0,rate)
    do k=1,40
        itc=k
        call collision()
        call collisionT()
        call streaming()
        call bounceback()
        call streamingT()
        call bouncebackT()
        call macroT()
        call macro()
    enddo
    !$acc wait
    call system_clock(clock1)
    call check_fields_fpc()
    call update_host_snapshot_2d_openacc()
    call update_host_reload_2d_openacc()
    open(unit=77,file='smoke.bin',access='sequential',form='unformatted',status='replace')
    write(77) nx,ny
    write(77) f
    write(77) g
    write(77) u
    write(77) v
    write(77) T
    write(77) rho
    close(77)
    write(*,*) 'SMOKE_ELAPSED_SECONDS',dble(clock1-clock0)/dble(rate)
    call exit_data_2d_openacc()
    end program main
'''

pbs_template = '''#!/usr/bin/env bash
#PBS -N {job_name}
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=240:00:00
#PBS -o {case_root}/results/pbs.stdout
#PBS -e {case_root}/results/pbs.stderr
set -euo pipefail
CASE_ROOT={case_root}
CAMPAIGN={remote}
RESULT="$CASE_ROOT/results"
BUILD="$CASE_ROOT/build"
mkdir -p "$BUILD" "$RESULT" "$CASE_ROOT/smoke/host" "$CASE_ROOT/smoke/gpu"
cd "$RESULT"
test "$(hostname -s)" = node07
test ! -e run.started
date -Is > run.started
finish() {{ rc=$?; printf '%s\n' "$rc" > "$RESULT/job.exit.status"; date -Is > "$RESULT/run.finished"; }}
trap finish EXIT
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${{LD_LIBRARY_PATH:-}}"
export OMP_NUM_THREADS=1
export ACC_DEVICE_TYPE=nvidia
export CUDA_VISIBLE_DEVICES="${{CUDA_VISIBLE_DEVICES:-0}}"
{{ date -Is; hostname; printf 'PBS_JOBID=%s\nPBS_NODEFILE=%s\n' "$PBS_JOBID" "$PBS_NODEFILE"; cat "$PBS_NODEFILE"; }} > job_identity.txt
nvidia-smi > gpu.start.txt 2>&1 || true
nvfortran --version > compiler_version.txt 2>&1
printf '%s\n' '-cpp -O3 -acc -gpu=cc60 -Mextend -Minfo=accel' > build_flags.txt
sha256sum "$CASE_ROOT/source/solver.F90" > runtime_source.sha256
printf '%s  %s\n' '{case_sha}' "$CASE_ROOT/source/solver.F90" | sha256sum -c -
python3 "$CAMPAIGN/gpu_gate.py" > gpu_gate.json
cd "$BUILD"
set +e
nvfortran -cpp -O3 -acc -gpu=cc60 -Mextend -Minfo=accel "$CASE_ROOT/source/solver.F90" -o solver.exe > "$RESULT/compile.log" 2>&1
rc=$?
set -e
printf '%s\n' "$rc" > "$RESULT/compile.status"
test "$rc" -eq 0
nvfortran -cpp -O2 -acc=host -Mextend "$CASE_ROOT/source/smoke.F90" -o smoke_host.exe > "$RESULT/smoke_host_compile.log" 2>&1
nvfortran -cpp -O3 -acc -gpu=cc60 -Mextend -Minfo=accel "$CASE_ROOT/source/smoke.F90" -o smoke_gpu.exe > "$RESULT/smoke_gpu_compile.log" 2>&1
( cd "$CASE_ROOT/smoke/host"; ACC_DEVICE_TYPE=host "$BUILD/smoke_host.exe" > solver.stdout 2>&1 )
( cd "$CASE_ROOT/smoke/gpu"; ACC_DEVICE_TYPE=nvidia "$BUILD/smoke_gpu.exe" > solver.stdout 2>&1 )
python3 "$CAMPAIGN/assess.py" --smoke "$CASE_ROOT" > "$RESULT/smoke_validation.json"
cd "$RESULT"
date -Is > solver.started
set +e
/usr/bin/time -v -o timing.txt "$BUILD/solver.exe" > solver.stdout 2> solver.stderr
rc=$?
set -e
printf '%s\n' "$rc" > run.status
date -Is > solver.finished
python3 "$CAMPAIGN/assess.py" --case "$CASE_ROOT" > assessment.json || {{ printf '%s\n' 87 > validation.status; exit 87; }}
printf '%s\n' 0 > validation.status
exit "$rc"
'''

cases = []
for mode, exponent, n in [('side',6,512),('side',7,512),('side',8,512),
                          ('rb',6,128),('rb',7,256),('rb',8,512),('rb',9,1024)]:
    name = '{}_Ra1e{}_N{}'.format('SideHeated' if mode=='side' else 'RB', exponent, n)
    folder = OUT / name
    folder.mkdir()
    (folder/'results').mkdir()
    text = baseline
    text = once(text,'nx=256, ny=256','nx={}, ny={}'.format(n,n))
    text = once(text,'Rayleigh=1.0d6','Rayleigh=1.0d{}'.format(exponent))
    if mode=='rb':
        text=once(text,'\n#define steadyFlow\n','\n!#define steadyFlow\n')
        text=once(text,'\n!#define unsteadyFlow\n','\n#define unsteadyFlow\n')
        text=once(text,'\n#define SideHeatedCell\n','\n!#define SideHeatedCell\n')
        text=once(text,'\n!#define RayleighBenardCell\n','\n#define RayleighBenardCell\n')
        text=once(text,'Prandtl=0.71d0','Prandtl=0.7d0')
    text=once(text,'itc_max=20000000','itc_max=nint(1000.0d0*timeUnit)')
    # No changes to collisions, forcing, boundaries, or post-processing.
    write(folder/'source/solver.F90',text)
    smoke,count=re.subn(r'^    program main\b.*?^    end program main\b',driver,text,flags=re.M|re.S)
    assert count==1
    write(folder/'source/smoke.F90',smoke)
    write(folder/'source/master_to_case.diff',''.join(difflib.unified_diff(
        baseline.splitlines(True),text.splitlines(True),fromfile='FPC_master',tofile=name)))
    case_sha=digest(text.encode('utf-8'))
    pr=.71 if mode=='side' else .7
    nu=.01*n*math.sqrt(pr/10**exponent)
    case=dict(name=name,case_type=mode,nx=n,ny=n,Ra=10**exponent,Pr=pr,
              Mach=math.sqrt(3)*.01,Uff=.01,Thot=2,Tcold=1,
              viscosity=nu,diffusivity=nu/pr,tauf=.5+3*nu,
              taug=.5+3*nu/pr,target_tff=1000,step_cap=100000*n,
              restart=False,node='node07',epsU=1e-7,epsT=1e-7,
              finite_check_interval_steps=2000,sample_interval_tff=10 if mode=='side' else .5,
              checkpoint_interval_tff=100,snapshot_enabled=True,tecplot_enabled=True,
              source_sha256=case_sha,smoke_sha256=digest(smoke.encode('utf-8')),
              remote_dir=REMOTE+'/'+name)
    if mode=='rb':
        case['reference_Nu'],case['reference_Re']={6:(6.30,279),7:(11.37,968),8:(25.25,3662),9:(53.51,15101)}[exponent]
    write(folder/'case.json',json.dumps(case,indent=2)+'\n')
    write(folder/'run.pbs',pbs_template.format(job_name='FPC_{}{}'.format('S' if mode=='side' else 'R',exponent),
            case_root=case['remote_dir'],remote=REMOTE,case_sha=case_sha))
    cases.append(case)

write(OUT/'source_master/2DRBOpenaccFPC.F90',baseline)
(OUT/'source_master/2DRBOpenaccFPC.original').write_bytes(raw)
protected=['FPC/2DRBOpenaccFPC.F90','均匀网格/2DRBOpenmp.F90','均匀网格/2DRBOpenacc.F90']
cfg=dict(created='2026-10-29',remote_dir=REMOTE,node='node07',walltime_each='240:00:00',
    latest_master=str(MASTER),master_sha256=digest(raw),cases=cases,
    reference_pdf=str(ROOT/'pdf/非稳态 RB 数据/Statistics of kinetic and thermal energy dissipation rates in two-dimensional turbulent Rayleigh–Bénard convection.pdf'),
    reference_pdf_page=5,reference_printed_page=169,
    reference_pdf_sha256=digest((ROOT/'pdf/非稳态 RB 数据/Statistics of kinetic and thermal energy dissipation rates in two-dimensional turbulent Rayleigh–Bénard convection.pdf').read_bytes()),
    dependencies='afterany between independent cases, in listed order',
    final_statistics='independent assessment uses Nu mean and sqrt(mean(Re(t)^2)); windows [500,750), [750,1000]',
    protected_sources={p:digest((ROOT/p).read_bytes()) for p in protected})
write(OUT/'campaign.json',json.dumps(cfg,indent=2,ensure_ascii=False)+'\n')
write(OUT/'README.md','''# FPC：node07 七工况测试（2026-10-29）

侧壁稳态：Ra=1e6/1e7/1e8，全部512²，Pr=0.71。RB非稳态：Ra=1e6/1e7/1e8/1e9，分别128²/256²/512²/1024²，Pr=0.7。
RB网格按用户要求，由Zhang、Zhou、Sun（2017）JFM原文表1对应的129/257/513/1025改用偶数尺寸。
均使用最新本地FPC源，Uff=0.01、Ma=sqrt(3)*0.01、Thot/Tcold=2/1、全部封闭无滑移；RB侧壁绝热、初始扰动1e-3。
侧壁epsU=epsT=1e-7，每2000步检测，收敛提前结束；1000TFF上限。RB目标1000TFF。
保留原输出设置：RB每0.5TFF记录Nu/Re与宏观快照，侧壁快照10TFF；检查点和Tecplot每100TFF。
独立算例从零启动，按afterany串行排在node07；前一工况失败不阻断后续独立工况。
每个作业先编译正式源和40步CPU/GPU对照副本，对照通过后才运行正式源。GPU上有其他进程时等待，避免并发挤占。
原算法与主源保持不变，参数化差异在各source/master_to_case.diff；三层源码SHA-256均留存。
原程序的非稳态汇总对瞬时Re取算术平均，本轮以独立assessment.json中的sqrt(mean(Re²))与JFM比较。
统计窗口分为500–750与750–1000TFF，半窗差异与文献误差分开报告；达到1000TFF不自动视为统计收敛。
实际编译、40步GPU对照、正式运行、有限性、达到时限、稳态/统计收敛分别判断。
所有正式数据位于/data2/XLLi/FPC/FPC_P100_node07_20261029；尚未提交时不写作业ID或运行成功结论。
''')
assert MASTER.read_bytes()==raw
print(json.dumps({'out':str(OUT),'cases':len(cases),'master_sha256':digest(raw)},indent=2))
