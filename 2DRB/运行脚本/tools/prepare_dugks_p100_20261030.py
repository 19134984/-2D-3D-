"""Prepare seven authorized DUGKS cases and a PBS-only GPU preflight."""
from pathlib import Path
import difflib
import hashlib
import json
import math
import re
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "运行脚本/run_cases/DUGKS_P100_node05_20261030"
REMOTE = "/data2/XLLi/DUGKS/DUGKS_P100_node05_20261030"
MASTER = ROOT / "DUGKS/2DRBOpenaccDUGKS.F90"
RAW = MASTER.read_bytes()
BASE = RAW.decode("utf-8-sig").replace("\r\n", "\n")
assert not (OUT / "campaign.json").exists(), "Do not overwrite an existing campaign."


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)


DRIVER = """program main
    use openacc
    use commondata
    implicit none
    integer :: k, unit
    integer(kind=8) :: c0,c1,rate
    call acc_init(acc_device_default)
    call initial()
    call enter_data_2d_openacc()
    do k=1,5
        itc=k
        call collision()
        call reconstruction_slopes()
        call face_flux_x()
        call face_flux_y()
        call streaming()
        call macro()
    enddo
    !$acc wait(1)
    call system_clock(c0,rate)
    do k=6,45
        itc=k
        call collision()
        call reconstruction_slopes()
        call face_flux_x()
        call face_flux_y()
        call streaming()
        call macro()
    enddo
    !$acc wait(1)
    call system_clock(c1)
    call check_fields()
    !$acc update self(f,g)
    open(newunit=unit,file='smoke.bin',form='unformatted',status='replace')
    write(unit) nx,ny
    write(unit) f
    write(unit) g
    write(unit) u
    write(unit) v
    write(unit) T
    write(unit) rho
    close(unit)
    write(*,*) 'SMOKE_ELAPSED_SECONDS',dble(c1-c0)/dble(rate)
    write(*,*) 'SMOKE_MEASURED_STEPS',40
    write(*,*) 'SMOKE_TOTAL_STEPS',45
    call exit_data_2d_openacc()
end program main
"""

ENV = """export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:$(printenv LD_LIBRARY_PATH || true)"
export OMP_NUM_THREADS=1
export ACC_DEVICE_TYPE=nvidia
export CUDA_VISIBLE_DEVICES="$(printenv CUDA_VISIBLE_DEVICES || printf 0)"
"""

cases = []
for kind, exponent, n in [
    ("side",6,512), ("side",7,512), ("side",8,512),
    ("rb",6,128), ("rb",7,256), ("rb",8,512), ("rb",9,1024)
]:
    name = ("SideHeated" if kind == "side" else "RB") + "_Ra1e%d_N%d" % (exponent,n)
    directory = OUT / name
    text = once(BASE,"nx=256, ny=256","nx=%d, ny=%d" % (n,n))
    text = once(text,"Rayleigh=1.0d6","Rayleigh=1.0d%d" % exponent)
    text = once(text,"flowCode=1, itc_max=20000000",
                "flowCode=1, itc_max=ceiling(1000.0d0*timeUnit/dt)")
    if kind == "rb":
        for old,new in [
            ("\n#define steadyFlow\n","\n!#define steadyFlow\n"),
            ("\n!#define unsteadyFlow\n","\n#define unsteadyFlow\n"),
            ("\n#define SideHeatedCell\n","\n!#define SideHeatedCell\n"),
            ("\n!#define RayleighBenardCell\n","\n#define RayleighBenardCell\n"),
            ("Prandtl=0.71d0","Prandtl=0.7d0")
        ]:
            text = once(text,old,new)
    smoke, count = re.subn(r"(?ms)^program main\s*$.*?^end program main\s*$",DRIVER,text)
    assert count == 1
    case = dict(name=name,case_type=kind,Ra=10**exponent,nx=n,ny=n,
                Pr=.71 if kind=="side" else .7,Mach=.01,Uff=.01/math.sqrt(3),
                Thot=.5,Tcold=-.5,Courant=.5,dt=.5/math.sqrt(2),
                target_tff=1000,average_windows=[[500,750],[750,1000]],
                epsU=1e-7,epsT=1e-7,consecutive_checks=3,check_interval_steps=2000,
                sample_interval_tff=10 if kind=="side" else .5,
                checkpoint_interval_tff=100,tecplot_interval_tff=100,
                snapshot_enabled=True,restart=False,node="node05",
                boundary="closed no-slip; adiabatic nonheated walls",
                remote_dir=REMOTE+"/"+name,source_sha256=sha(text.encode()),
                smoke_sha256=sha(smoke.encode()))
    case["step_cap"] = math.ceil(1000*n/case["Uff"]/case["dt"])
    case["estimated_device_array_MiB"] = math.ceil((108*n*n+36*n)*8/2**20+4*n*n*8/2**20)
    if kind == "rb":
        case["reference_Nu"],case["reference_Re"] = {
            6:(6.30,279),7:(11.37,968),8:(25.25,3662),9:(53.51,15101)}[exponent]
    write(directory/"source/solver.F90",text)
    write(directory/"source/smoke.F90",smoke)
    write(directory/"source/master_to_case.diff","".join(difflib.unified_diff(
        BASE.splitlines(True),text.splitlines(True),fromfile="latest_DUGKS_master",tofile=name)))
    write(directory/"case.json",json.dumps(case,indent=2)+"\n")
    (directory/"results").mkdir(exist_ok=True)
    cases.append(case)

write(OUT/"source_master/2DRBOpenaccDUGKS.F90",BASE)
(OUT/"source_master/2DRBOpenaccDUGKS.original").write_bytes(RAW)
paper = ROOT/"pdf/非稳态 RB 数据/Statistics of kinetic and thermal energy dissipation rates in two-dimensional turbulent Rayleigh–Bénard convection.pdf"
manifest = dict(date="2026-10-30",node="node05",remote_dir=REMOTE,
    master=str(MASTER),master_sha256=sha(RAW),normalized_master_sha256=sha(BASE.encode()),
    reference_pdf=str(paper),reference_pdf_sha256=sha(paper.read_bytes()),
    reference_table="Table 1, printed page 169 / PDF page 5, Pr=0.7",
    cases=cases,user_authorized_gpu_sharing=True,
    unchanged_solver_parameters=["Mach=.01","Courant=.5","Thot=.5","Tcold=-.5",
        "closed RB","output intervals","collision/forcing/face reconstruction/boundaries"],
    formal_dependencies="afterany, seven independent cases in listed order",
    validation="Separate compilation, host/GPU agreement, finite fields, completion, convergence and benchmark agreement.")
write(OUT/"campaign.json",json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")

preflight = f"""#!/usr/bin/env bash
#PBS -N DUGKS_GATE
#PBS -q batch
#PBS -l nodes=node05:ppn=1
#PBS -l walltime=01:00:00
#PBS -o {REMOTE}/preflight/pbs.stdout
#PBS -e {REMOTE}/preflight/pbs.stderr
set -euo pipefail
CAMPAIGN={REMOTE}
cd "$CAMPAIGN/preflight"
test "$(hostname -s)" = node05
test ! -e started
date -Is > started
finish() {{ rc=$?; printf '%s\\n' "$rc" > job.exit.status; date -Is > finished; }}
trap finish EXIT
{ENV}
{{ date -Is; hostname; echo "$PBS_JOBID"; cat "$PBS_NODEFILE"; }} > job_identity.txt
nvfortran --version > compiler_version.txt 2>&1
python3 "$CAMPAIGN/gpu_gate.py" > gpu_gate.json
python3 "$CAMPAIGN/preflight.py" "$CAMPAIGN" > summary.stdout 2>&1
"""
write(OUT/"preflight/run.pbs",preflight)
preflight_script = '''"""Compile every exact case and compare actual full-grid host/GPU steps."""
from pathlib import Path
import hashlib,json,os,subprocess,sys
from assess import smoke

root=Path(sys.argv[1])
cfg=json.loads((root/'campaign.json').read_text())
compiler='/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/bin/nvfortran'
reports=[]
for case in cfg['cases']:
    directory=root/case['name']
    build=directory/'build'
    build.mkdir(exist_ok=True)
    (directory/'results').mkdir(exist_ok=True)
    source=directory/'source/solver.F90'
    assert hashlib.sha256(source.read_bytes()).hexdigest()==case['source_sha256']
    (directory/'results/runtime_source.sha256').write_text(case['source_sha256']+'  '+str(source)+'\\n')
    flags=['-cpp','-O3','-acc','-gpu=cc60','-Mextend','-Minfo=accel']
    for name,src,extra in [
        ('solver.exe',source,flags),
        ('smoke_gpu.exe',directory/'source/smoke.F90',flags),
        ('smoke_host.exe',directory/'source/smoke.F90',['-cpp','-O2','-acc=host','-Mextend'])
    ]:
        logfile=directory/'results'/('compile.log' if name=='solver.exe' else name+'.compile.log')
        with logfile.open('w') as log:
            result=subprocess.run([compiler]+extra+[str(src),'-o',str(build/name)],cwd=str(build),stdout=log,stderr=subprocess.STDOUT)
        if name=='solver.exe':
            (directory/'results/compile.status').write_text(str(result.returncode)+'\\n')
        if result.returncode:
            raise RuntimeError('Compile failed: '+str(logfile))
    for device in ['host','gpu']:
        work=directory/'smoke'/device
        work.mkdir(parents=True,exist_ok=True)
        env=os.environ.copy()
        env['ACC_DEVICE_TYPE']='host' if device=='host' else 'nvidia'
        with (work/'solver.stdout').open('w') as log:
            subprocess.check_call([str(build/('smoke_'+device+'.exe'))],cwd=str(work),env=env,stdout=log,stderr=subprocess.STDOUT)
    report=smoke(directory)
    (directory/'results/smoke_validation.json').write_text(json.dumps(report,indent=2)+'\\n')
    reports.append(report)
    assert report['passed'],case['name']
    print(json.dumps(report),flush=True)
(root/'preflight/assessment.json').write_text(json.dumps(dict(passed=True,cases=reports),indent=2)+'\\n')
'''
write(OUT/"preflight.py",preflight_script)
(OUT/"assess.py").write_bytes((ROOT/"运行脚本/tools/dugks_campaign_assess.py").read_bytes())
write(OUT/"gpu_gate.py",'''"""User-authorized GPU sharing: verify free memory and record competing processes."""
import datetime,json,subprocess,time
checks=[]
consecutive=0
for attempt in range(61):
    gpu=subprocess.check_output(['nvidia-smi','--query-gpu=memory.free,utilization.gpu','--format=csv,noheader,nounits'],universal_newlines=True).strip()
    proc=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],universal_newlines=True).strip()
    free,util=[int(x.strip()) for x in gpu.splitlines()[0].split(',')]
    checks.append(dict(time=datetime.datetime.now().isoformat(),free_MiB=free,utilization=util,other_processes=proc))
    # Sharing is explicitly authorized; a resident process does not imply an idle GPU.
    consecutive=consecutive+1 if free>=4096 and util<15 else 0
    if consecutive>=3:
        print(json.dumps(dict(passed=True,sharing_authorized=True,checks=checks),indent=2))
        break
    time.sleep(5 if consecutive else 30)
else:
    print(json.dumps(dict(passed=False,checks=checks),indent=2))
    raise SystemExit('Insufficient free GPU resources; no other jobs were changed.')
''')
write(OUT/"README.md","""# DUGKS：node05 七工况测试（2026-10-30）

侧壁稳态：Ra=1e6/1e7/1e8，全部512²，Pr=0.71。
RB非稳态：Ra=1e6/1e7/1e8/1e9，128²/256²/512²/1024²，Pr=0.7。
网格由2017 JFM表1的129/257/513/1025按用户要求改为偶数。
RB采用原项目封闭方腔，四壁无滑移、左右绝热、上下恒温。
所有算例从零开始，目标1000TFF；侧壁达到连续三次epsU/epsT<1e-7可提前停止。
RB窗口500–750与750–1000TFF；Re统计采用sqrt(time mean(Re(t)^2))。
保持DUGKS母版Ma=.01、CFL=.5、冷热壁温度±.5及输出设置。
用户已明确允许在node05与低占用常驻GPU进程共享，七个正式任务用afterany串行。
先通过PBS进行NVHPC/P100全尺寸45步主机/GPU对照与计时，再设置正式任务walltime。
母版、规范化算例及PBS编译源的SHA-256分别记录；本轮不改求解器母版。
本目录不包含连接密码。编译/短程一致性/有限性/完成时长/统计收敛/文献一致性分别判断。
""")
initial=ROOT/"tmp/dugks_p100_20261030/node_check_initial.txt"
if initial.exists():
    (OUT/"pre_submit_node_check_initial.txt").write_bytes(initial.read_bytes())
assert MASTER.read_bytes()==RAW
print(json.dumps(dict(directory=str(OUT),master_sha256=sha(RAW),cases=len(cases)),indent=2))
