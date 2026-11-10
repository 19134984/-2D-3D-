"""Read-only P100 inspection; execute on master via Python stdin, no solver runs."""
from pathlib import Path
import datetime
import hashlib
import json
import math
import re
import subprocess
import time

root=Path('/data2/XLLi/FPC/FPC_P100_node07_20261029')
cfg=json.loads((root/'campaign.json').read_text())
now=datetime.datetime.now().astimezone()
report=dict(checked_at=now.isoformat(),cases=[],completed=0)

def read_rows(path):
    if not path.exists():
        return []
    return [tuple(map(float,line.split())) for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]

for case in cfg['cases']:
    folder=root/case['name']
    result=folder/'results'
    job=(folder/'job.id').read_text().strip()
    query=subprocess.run(['qstat','-f',job],stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                         universal_newlines=True)
    state=re.search(r'job_state\s*=\s*(\w)',query.stdout)
    entry=dict(name=case['name'],Ra=case['Ra'],nx=case['nx'],job_id=job,
               queue_state=state.group(1) if state else None)
    entry['source_hash_matches']=hashlib.sha256((folder/'source/solver.F90').read_bytes()).hexdigest()==case['source_sha256']
    runtime=result/'runtime_source.sha256'
    entry['runtime_hash_matches']=runtime.exists() and runtime.read_text().split()[0]==case['source_sha256']
    assert entry['source_hash_matches'] and entry['runtime_hash_matches']
    entry['diagnostic_errors']={}
    for filename in ['solver.stdout','solver.stderr','pbs.stderr']:
        path=result/filename
        if path.exists():
            entry['diagnostic_errors'][filename]=[line[:300] for line in path.read_text(errors='replace').splitlines()
                if re.search(r'(?i)(?:\b(?:nan|infinity|nonfinite|cuda_error)\b|error\s+stop|^\s*error:)',line)][:8]
    for attempt in range(3):
        nu=read_rows(result/'Nu_VolAvg_2DOpenaccFPC.dat')
        reynolds=read_rows(result/'Re_VolAvg_2DOpenaccFPC.dat')
        if len(nu)==len(reynolds):
            break
        time.sleep(.05)
    count=min(len(nu),len(reynolds))
    nu=nu[:count]; reynolds=reynolds[:count]
    assert all(len(row)==2 and all(math.isfinite(v) for v in row) for row in nu+reynolds)
    assert all(a[0]==b[0] for a,b in zip(nu,reynolds))
    assert all(a[0]<b[0] for a,b in zip(nu,nu[1:]))
    if count:
        entry['history']=dict(samples=count,latest_time=nu[-1][0],latest_Nu=nu[-1][1],
                              latest_Re=reynolds[-1][1],finite_monotone_aligned=True)
    meta_path=result/'reloadFile2DOpenaccFPC-latest.meta'
    if meta_path.exists():
        meta=dict(line.split(None,1) for line in meta_path.read_text().splitlines())
        assert int(meta['nx'])==case['nx'] and int(meta['ny'])==case['ny']
        entry['latest_checkpoint_tff']=float(meta['time_tf'])
        entry['latest_checkpoint_file_exists']=(result/('reloadFile2DOpenaccFPC-'+meta['reloadFileName']+'.bin')).exists()
    status=result/'run.status'
    if status.exists():
        entry['solver_exit_code']=int(status.read_text())
        entry['job_exit_code']=int((result/'job.exit.status').read_text()) if (result/'job.exit.status').exists() else None
        entry['validation_exit_code']=int((result/'validation.status').read_text()) if (result/'validation.status').exists() else None
        assessment=result/'assessment.json'
        if assessment.exists():
            entry['assessment']=json.loads(assessment.read_text())
        entry['status']='finished'
        report['completed']+=1
    else:
        entry['status']='running' if entry['queue_state']=='R' else 'incomplete'
        started=result/'solver.started'
        if started.exists() and count and case['case_type']=='rb':
            start=datetime.datetime.strptime(started.read_text().strip(),'%Y-%m-%dT%H:%M:%S%z')
            elapsed=(now-start).total_seconds()/3600
            current=nu[-1][0]
            entry.update(elapsed_hours=elapsed,progress_percent=100*current/case['target_tff'],
                         estimated_remaining_hours=max(0,case['target_tff']-current)*elapsed/current)
    if case['case_type']=='rb' and count:
        assert all(abs(t-(i+1)*case['sample_interval_tff'])<1e-9 for i,(t,_) in enumerate(nu))
        if nu[-1][0]>=1000:
            def window(lo,hi,inclusive):
                ids=[i for i,(t,_) in enumerate(nu) if t>=lo and (t<=hi if inclusive else t<hi)]
                return dict(count=len(ids),Nu=math.fsum(nu[i][1] for i in ids)/len(ids),
                            Re=math.sqrt(math.fsum(reynolds[i][1]**2 for i in ids)/len(ids)))
            first=window(500,750,False); last=window(750,1000,True)
            difference={k:100*abs(first[k]-last[k])/abs(last[k]) for k in ['Nu','Re']}
            entry['independent_statistics']=dict(first_half=first,second_half=last,
                half_window_difference_percent=difference,half_window_passed=all(v<1 for v in difference.values()),
                reference_error_percent={k:100*(last[k]-case['reference_'+k])/case['reference_'+k] for k in ['Nu','Re']})
    report['cases'].append(entry)
report['queue']=subprocess.run(['qstat','-n','-u','xlli'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                               universal_newlines=True).stdout
report['gpu']=subprocess.run(['ssh','node07',
    'nvidia-smi --query-gpu=memory.free,utilization.gpu --format=csv; nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv'],
    stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True).stdout
print(json.dumps(report))
