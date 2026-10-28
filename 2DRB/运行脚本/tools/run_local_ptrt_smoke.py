"""Compile all formal selections and save two complete host-field references."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'运行脚本/run_cases/LuoPTRT_Ra9_Ra10_20261025'
COMPILER='C:/msys64/ucrt64/bin/gfortran.exe'
os.environ['PATH']=str(Path(COMPILER).parent)+os.pathsep+os.environ['PATH']
os.environ['ACC_DEVICE_TYPE']='host'
cfg=json.loads((OUT/'campaign.json').read_text(encoding='utf-8'))
report={'syntax':{},'smoke':{}}
with tempfile.TemporaryDirectory(prefix='ptrt_full_',delete=False) as tmp:
    build=Path(tmp)
    for name in [*cfg['order'],'GPU_P1','GPU_P2','legacy_fallback']:
        source=OUT/(name if name!='legacy_fallback' else 'P1_Ra1e9')/'source/solver.F90'
        text=source.read_text(encoding='utf-8')
        if name=='legacy_fallback': text=text.replace('enablePTRT=.true.','enablePTRT=.false.')
        path=build/(name+'.F90'); path.write_text(text,encoding='utf-8',newline='\n')
        flags=[COMPILER,'-cpp','-fopenacc','-ffree-line-length-none']
        res=subprocess.run(flags+['-fsyntax-only',str(path)],cwd=build,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        if res.returncode: raise RuntimeError(name+'\n'+res.stdout)
        report['syntax'][name]='PASS'
        if name.startswith('GPU_P'):
            exe=build/(name+'.exe')
            res=subprocess.run(flags+['-O2','-fcheck=bounds',str(path),'-o',str(exe)],cwd=build,
                stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            if res.returncode: raise RuntimeError(res.stdout)
            result=build/name; result.mkdir()
            res=subprocess.run([str(exe)],cwd=result,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
            (result/'solver.stdout').write_bytes(res.stdout)
            assert res.returncode==0,(name,res.returncode,res.stdout[-1000:])
            reference=OUT/name/'host_reference'
            snapshots=list(result.glob('*Snapshot*.bin')); assert len(snapshots)==1
            for p in [*snapshots,result/'solver.stdout',result/'SimulationSettings2DOpenaccLBMCDE_D2Q5.txt']:
                shutil.copyfile(p,reference/p.name)
            report['smoke'][name]={'pass':True,'grid':[32,32],'Ra':1e5,'completed_tff':2,
                'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                'snapshot_sha256':hashlib.sha256(snapshots[0].read_bytes()).hexdigest()}
    report['temporary_build']=str(build)
(OUT/'full_source_host_validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
