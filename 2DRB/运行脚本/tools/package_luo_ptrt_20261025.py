"""Prepare four immutable cases and two small GPU/host comparisons."""
from pathlib import Path
import difflib
import hashlib
import json
import math
import shutil
import tarfile

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'运行脚本/run_cases/LuoPTRT_Ra9_Ra10_20261025_r2'
MASTER=ROOT/'PTRT/2DRBOpenaccLuoPTRT.F90'
REMOTE='/data2/XLLi/PTRT/LuoPTRT_Ra9_Ra10_20261025_r2'
NODE='node05'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p,text):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(text,encoding='utf-8',newline='\n')


def change(s,old,new):
    assert s.count(old)==1,(old,s.count(old))
    return s.replace(old,new)


def prepare():
    source=MASTER.read_text(encoding='utf-8')
    baseline=ROOT/'Xs/D2Q9TRT_D2Q5LuoTRT/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90'
    protected=json.loads((OUT/'protected_sources.json').read_text(encoding='utf-8'))
    assert all(sha(ROOT/p)==h for p,h in protected.items())
    validation=json.loads((OUT/'local_validation.json').read_text(encoding='utf-8'))
    assert validation['pass'] and validation['source_sha256']==sha(MASTER),'Validate the current source first'
    (OUT/'source_master'/baseline.name).write_bytes(baseline.read_bytes())
    cases={}
    order=['P1_Ra1e9','P2_Ra1e9','P1_Ra1e10','P2_Ra1e10']
    for name in order:
        strategy=int(name[1]); n=2048 if name.endswith('1e10') else 1024
        ra=1e10 if n==2048 else 1e9
        text=source
        if strategy==2: text=change(text,'ptrtSqStrategy=1','ptrtSqStrategy=2')
        if n==2048:
            text=change(text,'nx=1024, ny=1024','nx=2048, ny=2048')
            text=change(text,'Rayleigh=1.0d9','Rayleigh=1.0d10')
        path=OUT/name/'source/solver.F90'
        assert not (OUT/name/'results/run.started').exists(),'Never replace a running case'
        write(path,text); (OUT/name/'results').mkdir(exist_ok=True)
        write(OUT/name/'source/master_to_case.diff',''.join(difflib.unified_diff(source.splitlines(True),
            text.splitlines(True),fromfile='PTRT/2DRBOpenaccLuoPTRT.F90',tofile=name+'/source/solver.F90')))
        nu=.1*n*math.sqrt(.7/(3*ra)); sn=1/(.5+3*nu)
        sq=8*(2-sn)/(8-sn) if strategy==1 else 1.6
        cases[name]=dict(nx=n,ny=n,Ra=ra,Pr=.7,Ma=.1,chi_nu=0,chi_b=0,viscosity=nu,diffusivity=nu/.7,
            s_e=sn,s_nu=sn,s_q=sq,tau_nu=1/sn,tau_q=1/sq,Lambda_s=(1/sn-.5)*(1/sq-.5),
            thermal_paraA=20*math.sqrt(3)*(nu/.7)-4,target_tff=1000,from_scratch=True,
            history_interval_tff=1,collision_diagnostic_interval_tff=.1,average_first=[500,750],
            average_final=[750,1000],Nu_reference=93.88 if n==2048 else 53.51,
            Re_reference=61797 if n==2048 else 15101,source_sha256=sha(path),
            release_requires=name[:2]+'_Ra1e9' if n==2048 else 'GPU_validation')
    for strategy in [1,2]:
        text=change(source,'nx=1024, ny=1024','nx=32, ny=32')
        text=change(text,'Rayleigh=1.0d9','Rayleigh=1.0d5')
        text=change(text,'unsteadyRunDuration=1000.0d0','unsteadyRunDuration=2.0d0')
        text=change(text,'outputSnapshotInterval=20.0d0','outputSnapshotInterval=2.0d0')
        if strategy==2: text=change(text,'ptrtSqStrategy=1','ptrtSqStrategy=2')
        write(OUT/f'GPU_P{strategy}/source/solver.F90',text)
        (OUT/f'GPU_P{strategy}/results').mkdir(exist_ok=True)
        (OUT/f'GPU_P{strategy}/host_reference').mkdir(exist_ok=True)
    cfg=dict(remote_dir=REMOTE,node=NODE,master=str(MASTER),master_sha256=sha(MASTER),
        parent=str(baseline),parent_sha256=sha(baseline),order=order,cases=cases,
        references=dict(PTRT='https://arxiv.org/abs/2602.06686v1',
            Luo='https://doi.org/10.1016/j.jcp.2014.06.047',RB='https://doi.org/10.1017/jfm.2017.19'),
        release_policy=dict(half_window_percent=1.0,wall_budget_relative=.05,Nu_reference_relative=.10,
            Re_reference_relative=.20,require_full_horizon_and_finite=True),
        note='1024/2048 fluid nodes explicitly requested by user; literature grids are 1025/2049.')
    write(OUT/'campaign.json',json.dumps(cfg,indent=2,ensure_ascii=False)+'\n')
    write(OUT/'master_to_ptrt.diff',''.join(difflib.unified_diff(baseline.read_text(encoding='utf-8-sig').splitlines(True),
        source.splitlines(True),fromfile=baseline.relative_to(ROOT).as_posix(),tofile=MASTER.relative_to(ROOT).as_posix())))
    for src,dst in [('运行脚本/tools/assess_luo_ptrt.py','assess_case.py'),
                    ('运行脚本/tools/check_ptrt_gpu_smoke.py','check_gpu_smoke.py'),
                    ('运行脚本/tools/run_luo_ptrt_case.sh','run_case.sh'),
                    ('后处理/audit_flow_history.py','audit_flow_history.py')]:
        write(OUT/dst,(ROOT/src).read_text(encoding='utf-8-sig'))
    for name in ['GPU_validation',*order]:
        (OUT/name/'results').mkdir(parents=True,exist_ok=True)
        duration='00:30:00' if name=='GPU_validation' else '96:00:00'
        job_name='PTRT_GPUtest' if name=='GPU_validation' else 'PTRT_'+name
        write(OUT/name/'run.pbs',f'''#!/bin/bash
#PBS -N {job_name}
#PBS -q batch
#PBS -l nodes={NODE}:ppn=1
#PBS -l walltime={duration}
#PBS -o {REMOTE}/{name}/results/pbs.stdout
#PBS -e {REMOTE}/{name}/results/pbs.stderr
set -euo pipefail
exec bash {REMOTE}/run_case.sh {name}
''')
    write(OUT/'submit.sh',f'''#!/bin/bash
set -euo pipefail
cd {REMOTE}
test ! -s submission.tsv
sha256sum -c upload.sha256 > upload_verification.log
previous=$(qsub GPU_validation/run.pbs)
printf 'GPU_validation\t%s\tNONE\n' "$previous" > submission.tsv
for name in P1_Ra1e9 P2_Ra1e9 P1_Ra1e10 P2_Ra1e10; do
    dependency=afterany:$previous
    job=$(qsub -W depend="$dependency" "$name/run.pbs")
    printf '%s\t%s\t%s\n' "$name" "$job" "$dependency" >> submission.tsv
    previous=$job
done
cat submission.tsv
''')
    print(json.dumps(dict(cases=cases,remote_dir=REMOTE,node=NODE),indent=2))


def archive():
    cfg=json.loads((OUT/'campaign.json').read_text(encoding='utf-8'))
    assert sha(MASTER)==cfg['master_sha256']
    for name,case in cfg['cases'].items(): assert sha(OUT/name/'source/solver.F90')==case['source_sha256']
    for strategy in [1,2]:
        assert list((OUT/f'GPU_P{strategy}/host_reference').glob('*Snapshot*.bin')),'Run host smoke first'
    paths=[p for p in OUT.rglob('*') if p.is_file() and 'results' not in p.relative_to(OUT).parts
           and p.name not in ['campaign.tar.gz','upload.sha256','submission.json','submission.tsv','submission_qstat.txt']
           and p.suffix not in ['.exe','.mod','.pyc']]
    write(OUT/'upload.sha256',''.join(sha(p)+'  '+p.relative_to(OUT).as_posix()+'\n' for p in sorted(paths)))
    with tarfile.open(OUT/'campaign.tar.gz','w:gz') as tar:
        for p in paths+[OUT/'upload.sha256']: tar.add(p,arcname=p.relative_to(OUT).as_posix(),recursive=False)
        for name in ['GPU_validation','GPU_P1','GPU_P2',*cfg['order']]:
            info=tarfile.TarInfo(name+'/results'); info.type=tarfile.DIRTYPE; info.mode=0o755; tar.addfile(info)
    print(json.dumps(dict(archive=str(OUT/'campaign.tar.gz'),sha256=sha(OUT/'campaign.tar.gz'),files=len(paths)+1)))


if __name__=='__main__':
    import sys
    archive() if '--archive' in sys.argv else prepare()
