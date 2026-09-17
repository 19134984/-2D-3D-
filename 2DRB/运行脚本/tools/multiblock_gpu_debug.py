"""P100 short regression jobs. Credentials come only from SSHPASS."""
from pathlib import Path
import subprocess, sys, tarfile

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/'运行脚本/Multiblock_gpu_debug'
REMOTE='/data2/XLLi/Multiblock/gpu_debug'
SSH=['C:/msys64/usr/bin/sshpass.exe','-e','C:/msys64/usr/bin/ssh.exe','-o','ConnectTimeout=15','xlli@10.169.174.118']
action=sys.argv[1]
name=sys.argv[2]
folder=LOCAL/name
remote=REMOTE+'/'+name
if action=='prepare':
    import re, hashlib
    src=Path(sys.argv[3])
    s=src.read_text(encoding='utf-8-sig')
    s=re.sub(r'(parameter :: itc_max = )20000000',r'\g<1>2',s)
    s=re.sub(r'(parameter :: output(?:Snapshot|Plt|Reload)File = )1',r'\g<1>0',s)
    folder.mkdir(parents=True,exist_ok=False)
    (folder/'solver.F90').write_text(s,encoding='utf-8',newline='\n')
    (folder/'source.sha256').write_text(hashlib.sha256((folder/'solver.F90').read_bytes()).hexdigest()+'  solver.F90\n',newline='\n')
    pbs='''#!/bin/bash
#PBS -N MBdebug
#PBS -q batch
#PBS -l nodes=node05:ppn=1
#PBS -l walltime=00:20:00
#PBS -o __REMOTE__/pbs.stdout
#PBS -e __REMOTE__/pbs.stderr
set -eu
cd __REMOTE__
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/bin:$PATH
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=0 ACC_DEVICE_TYPE=nvidia
sha256sum -c source.sha256
nvfortran -cpp -O3 -acc -gpu=cc60,lineinfo -Minfo=accel -Mextend solver.F90 -o solver.exe > compile.log 2>&1
set +e
./solver.exe > solver.stdout 2>&1
rc=$?
echo "$rc" > run.status
exit "$rc"
'''.replace('__REMOTE__',remote)
    (folder/'run.pbs').write_text(pbs,encoding='utf-8',newline='\n')
elif action=='upload':
    archive=LOCAL/(name+'.tar')
    with tarfile.open(archive,'w') as t:
        t.add(folder,arcname=name)
    subprocess.run(SSH+[f"mkdir -p {REMOTE}; test ! -e {remote}/run.pbs && tar -xf - -C {REMOTE}"],input=archive.read_bytes(),check=True)
elif action=='submit':
    r=subprocess.run(SSH+[f"cd {remote} && bash -n run.pbs && test ! -e job.id && qsub run.pbs"],capture_output=True,check=True)
    job=r.stdout.decode().strip()
    (folder/'job.id').write_text(job)
    subprocess.run(SSH+[f"printf '%s\\n' '{job}' > {remote}/job.id"],check=True)
    print(job)
elif action=='collect':
    job=(folder/'job.id').read_text().strip()
    q=subprocess.run(SSH+['qstat -f '+job],capture_output=True)
    (folder/'qstat.txt').write_bytes(q.stdout)
    print(q.stdout.decode())
    for f in ['run.status','solver.stdout','compile.log','pbs.stdout','pbs.stderr','sync.stdout','sync.status','sanitizer.stdout','sanitizer.status','SimulationSettings2DOpenaccMultiblock.txt']:
        r=subprocess.run(SSH+[f'cat {remote}/{f}'],capture_output=True)
        if r.returncode==0:
            (folder/f).write_bytes(r.stdout)
            if f.endswith('.status') or f in ['solver.stdout','sync.stdout','sanitizer.stdout']:
                print(f, r.stdout.decode(errors='replace')[-12000:])
elif action=='compare':
    import numpy as np, json
    r=subprocess.run(SSH+[f'cat {remote}/state.bin'],capture_output=True,check=True)
    gpu=np.frombuffer(r.stdout,dtype='<f8')
    host=np.fromfile(LOCAL/'host_state.bin',dtype='<f8')
    assert gpu.shape==host.shape
    assert np.all(np.isfinite(gpu)) and np.all(np.isfinite(host))
    report=dict(values=len(gpu),max_absolute_error=float(np.max(abs(gpu-host))),
                tolerance=1e-11,passed=bool(np.allclose(gpu,host,rtol=1e-11,atol=1e-11)))
    (folder/'gpu_host_comparison.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report))
    assert report['passed']
else: raise SystemExit('unknown action')
