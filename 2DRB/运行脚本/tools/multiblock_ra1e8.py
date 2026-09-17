"""Submit and inspect the authorized Ra=1e8 multiblock case; credentials via SSHPASS."""
from pathlib import Path
import subprocess,sys,json,hashlib,tarfile,io,re

ROOT=Path(__file__).resolve().parents[2]
REL='Ra1e8/ratio2/sideheated_N512_20260916_v3'
LOCAL=ROOT/'运行脚本/Multiblock_sideheated_20260916_v3'/REL
REMOTE='/data2/XLLi/Multiblock/'+REL
SSH=['C:/msys64/usr/bin/sshpass.exe','-e','C:/msys64/usr/bin/ssh.exe','-o','ConnectTimeout=15','xlli@10.169.174.118']
action=sys.argv[1]
if action=='prepare':
    meta=json.loads((LOCAL/'manifest.json').read_text())
    current=(ROOT/'多块网格/2DRBOpenaccMultiblock.F90').read_bytes()
    assert hashlib.sha256(current).hexdigest()==meta['local_master_sha256'],'Master changed: regenerate case first'
    assert current==(LOCAL/'source/local_master.F90').read_bytes()
    assert hashlib.sha256((LOCAL/'source/solver.F90').read_bytes()).hexdigest()==meta['case_source_sha256']
    s=(LOCAL/'source/solver.F90').read_text(encoding='utf-8')
    enabled=set(re.findall(r'^#define\s+(\w+)',s,re.M))
    assert {'steadyFlow','SideHeatedCell','HorizontalWallsAdiabatic','VerticalWallsConstT'}<=enabled
    assert not {'unsteadyFlow','RayleighBenardCell','EnableUseG'}&enabled
    assert re.search(r'nx\s*=\s*512,\s*ny\s*=\s*512',s)
    assert re.search(r'loadInitField\s*=\s*0',s)
    pbs='''#!/bin/bash
#PBS -N MB8N512R2
#PBS -q batch
#PBS -l nodes=node05:ppn=1
#PBS -l walltime=96:00:00
#PBS -o __CASE__/results/pbs.stdout
#PBS -e __CASE__/results/pbs.stderr
set -eu
cd __CASE__/results
exec bash __CASE__/run.sh
'''.replace('__CASE__',REMOTE)
    (LOCAL/'run.pbs').write_text(pbs,encoding='utf-8',newline='\n')
    assert b'\r' not in (LOCAL/'run.sh').read_bytes()
    print(meta['case_source_sha256'])
elif action=='upload':
    buf=io.BytesIO()
    with tarfile.open(fileobj=buf,mode='w') as t:
        for p in ['run.pbs','run.sh','manifest.json','source']:
            t.add(LOCAL/p,arcname=p)
    command=f'mkdir -p {REMOTE}/results && test ! -e {REMOTE}/job.id && test ! -e {REMOTE}/results/run.started && tar -xf - -C {REMOTE} && bash -n {REMOTE}/run.pbs && bash -n {REMOTE}/run.sh && sha256sum {REMOTE}/source/solver.F90'
    subprocess.run(SSH+[command],input=buf.getvalue(),check=True)
elif action=='submit':
    r=subprocess.run(SSH+[f'cd {REMOTE} && test ! -e job.id && test ! -e results/run.started && qsub run.pbs'],capture_output=True,check=True)
    job=r.stdout.decode().strip()
    assert re.fullmatch(r'\d+\.master',job),job
    (LOCAL/'job.id').write_text(job)
    subprocess.run(SSH+[f"printf '%s\\n' '{job}' > {REMOTE}/job.id"],check=True)
    print(job)
elif action=='inspect':
    job=(LOCAL/'job.id').read_text().strip()
    r=subprocess.run(SSH+['qstat -f '+job],capture_output=True)
    (LOCAL/'results/qstat.txt').write_bytes(r.stdout)
    print(r.stdout.decode(errors='replace'))
    for name in ['compile.status','runtime_source.sha256','job_identity.txt','solver.stdout',
                 'run.status','SimulationSettings2DOpenaccMultiblock.txt','compile.log',
                 'Convergence_2DOpenaccMultiblock.dat','NuRe_2DOpenaccMultiblock.dat']:
        r=subprocess.run(SSH+[f'cat {REMOTE}/results/{name}'],capture_output=True)
        if r.returncode==0:
            (LOCAL/'results'/name).write_bytes(r.stdout)
            if name in ['compile.status','runtime_source.sha256','job_identity.txt','solver.stdout','run.status']:
                print(name,r.stdout.decode(errors='replace')[-1500:])
elif action=='audit-inputs':
    for name in ['validation.status','timing.txt','pbs.stdout','pbs.stderr',
                 'reloadFile2DOpenaccMultiblock-latest.meta','run.started','run.finished']:
        r=subprocess.run(SSH+[f'cat {REMOTE}/results/{name}'],capture_output=True,check=True)
        (LOCAL/'results'/name).write_bytes(r.stdout)
        print(name,r.stdout.decode(errors='replace')[-700:])
else: raise SystemExit('unknown action')
