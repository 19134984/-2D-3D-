"""Package the prepared MRT campaign and verify all protected solver hashes."""
from pathlib import Path
import hashlib
import json
import tarfile

root = Path(__file__).resolve().parents[2]
out = root/'运行脚本/run_cases/Uniform_LuoMRT_Ra9_Ra10_T1000_20261025'
cfg = json.loads((out/'campaign.json').read_text(encoding='utf-8'))
remote = cfg['remote_dir']
protected = json.loads((out/'protected_sources.json').read_text(encoding='utf-8'))
for rel, expected in protected.items():
    assert hashlib.sha256((root/rel).read_bytes()).hexdigest() == expected, rel
for name in cfg['order']:
    expected = cfg['cases'][name]['source_sha256']
    assert hashlib.sha256((out/name/'source/solver.F90').read_bytes()).hexdigest() == expected
    pbs = f'''#!/bin/bash
#PBS -N LuoMRT_{name}
#PBS -q batch
#PBS -l nodes=node07:ppn=1
#PBS -l walltime=96:00:00
#PBS -o {remote}/{name}/results/pbs.stdout
#PBS -e {remote}/{name}/results/pbs.stderr
set -euo pipefail
exec bash {remote}/run_case.sh {name}
'''
    (out/name/'run.pbs').write_text(pbs,encoding='utf-8',newline='\n')
paths = [p for p in out.rglob('*') if p.is_file() and
         'results' not in p.relative_to(out).parts and
         p.name not in ('upload.sha256','campaign.tar.gz','submission.txt','submission_qstat.txt')
         and p.suffix not in ('.pyc',)]
manifest = ''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(out).as_posix()+'\n'
                   for p in sorted(paths))
(out/'upload.sha256').write_text(manifest,encoding='ascii',newline='\n')
with tarfile.open(out/'campaign.tar.gz','w:gz') as archive:
    for p in paths+[out/'upload.sha256']:
        archive.add(p,arcname=p.relative_to(out).as_posix(),recursive=False)
    for name in ['Ra1e9','Ra1e10','smoke']:
        info = tarfile.TarInfo(name+'/results')
        info.type = tarfile.DIRTYPE
        info.mode = 0o755
        archive.addfile(info)
print(json.dumps(dict(archive=str(out/'campaign.tar.gz'),files=len(paths)+1,
    archive_sha256=hashlib.sha256((out/'campaign.tar.gz').read_bytes()).hexdigest(),
    protected_sources_unchanged=True),ensure_ascii=False))
