"""Run on master after a successful gate; reconnect existing RB jobs without changing their files."""
from pathlib import Path
import datetime,hashlib,json,re,subprocess,time
root=Path('/data2/XLLi/DUGKS/DUGKS_Side256_P100_node05_20261031')
old=Path('/data2/XLLi/DUGKS/DUGKS_P100_node05_20261030')
audit=root/'replacement_audit'
cfg=json.loads((root/'campaign.json').read_text())
def run(args,check=True):
    p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    with (audit/'operations.log').open('a') as f:
        f.write(datetime.datetime.now().isoformat()+' '+repr(args)+'\n'+p.stdout+'\nrc='+str(p.returncode)+'\n')
    if check and p.returncode: raise RuntimeError(repr(args)+' '+p.stdout)
    return p
def fields(text):
    data={}; key=None
    for line in text.splitlines():
        if line.startswith('    ') and ' = ' in line:
            key,val=line.strip().split(' = ',1); data[key]=val
        elif key and line.startswith('\t'): data[key]+=line.strip()
    return data
def detail(job): return run(['qstat','-f',job]).stdout
def nodecheck(path):
    checks=[]
    for args in [['date','-Is'],['qstat','-a'],['pbsnodes','node05','node07']]:
        checks.append(dict(command=args,output=run(args).stdout))
    for node in ['node05','node07']:
        cmd='date -Is; nvidia-smi --query-gpu=memory.total,memory.used,memory.free,utilization.gpu --format=csv; nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv'
        checks.append(dict(node=node,output=run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10',node,cmd]).stdout))
    path.write_text(json.dumps(dict(time=datetime.datetime.now().isoformat(),checks=checks,reason='User-selected node05, shared execution authorized; new steady cases remain serial and existing RB follows them. Node07 is occupied by FPC.'),indent=2)+'\n')
    gpu=next(x['output'] for x in checks if x.get('node')=='node05')
    values=next(line for line in gpu.splitlines() if re.match(r'^\d+ MiB,',line))
    assert int(values.split(',')[2].strip().split()[0])>=1024
assert (audit/'cancellation.finished').exists()
gate=json.loads((root/'preflight/assessment.json').read_text())
assert gate['passed'] and len(gate['cases'])==3
assert (root/'preflight/job.exit.status').read_text().strip()=='0'
assert not (root/'submission.started').exists(), 'Inspect existing submission; never duplicate jobs'
before=json.loads((audit/'original_file_hashes.json').read_text())
for name,files in before.items():
    for path,expected in files.items():
        assert hashlib.sha256((old/name/path).read_bytes()).hexdigest()==expected
for j in range(7467,7471):
    f=fields(detail(str(j)+'.master'))
    assert f['job_state']=='H' and 'u' in f['Hold_Types']
(root/'submission.started').write_text(datetime.datetime.now().isoformat()+'\n')
submitted=[]
previous=None
for case in cfg['cases']:
    folder=root/case['name']
    assert hashlib.sha256((folder/'source/solver.F90').read_bytes()).hexdigest()==case['source_sha256']
    smoke=json.loads((folder/'results/smoke_validation.json').read_text())
    assert smoke['passed'] and smoke['source_sha256']==case['source_sha256']
    nodecheck(folder/'pre_submit_node_check.json')
    run(['bash','-n',str(folder/'run.pbs')])
    args=['qsub','-h']
    if previous: args+=['-W','depend=afterany:'+previous]
    args.append(str(folder/'run.pbs'))
    job=run(args).stdout.strip()
    assert re.fullmatch(r'\d+\.master',job)
    item=dict(case=case['name'],job_id=job,dependency=previous,node='node05',
              walltime=case['walltime'],directory=str(folder),source_sha256=case['source_sha256'])
    submitted.append(item)
    (root/'submission.json').write_text(json.dumps(submitted,indent=2)+'\n')
    text=detail(job); (folder/'qstat_after_submit.txt').write_text(text)
    f=fields(text)
    assert f['Resource_List.nodes']=='node05:ppn=1' and f['job_state']=='H' and 'u' in f['Hold_Types']
    previous=job
# Preserve the existing first->second RB edge explicitly while replacing only its incoming edge.
run(['qalter','-W','depend=afterany:'+previous+',beforeany:7468.master','7467.master'])
oldcfg=json.loads((old/'campaign.json').read_text())
oldjobs=json.loads((old/'submission.json').read_text())
chain=list(submitted)
for case,item in zip(oldcfg['cases'][3:],oldjobs[3:]):
    f=fields(detail(item['job_id']))
    assert ('afterany:'+previous) in f.get('depend',''), (item['job_id'],f.get('depend'))
    assert f['Resource_List.nodes']=='node05:ppn=1' and 'u' in f['Hold_Types']
    assert hashlib.sha256((old/case['name']/'source/solver.F90').read_bytes()).hexdigest()==case['source_sha256']
    chain.append(dict(case=case['name'],job_id=item['job_id'],dependency=previous,node='node05',
                      nx=case['nx'],ny=case['ny'],directory=case['remote_dir'],source_sha256=case['source_sha256']))
    previous=item['job_id']
(root/'active_chain.json').write_text(json.dumps(chain,indent=2)+'\n')
# Downstream tasks retain dependency holds. Release the first new task last.
for item in reversed(chain[1:]): run(['qrls','-h','u',item['job_id']])
run(['qrls','-h','u',chain[0]['job_id']])
for name,files in before.items():
    for path,expected in files.items():
        assert hashlib.sha256((old/name/path).read_bytes()).hexdigest()==expected
state={}
for item in chain:
    text=detail(item['job_id']); f=fields(text)
    assert 'u' not in f['Hold_Types']
    assert f['Resource_List.nodes']=='node05:ppn=1'
    if f.get('exec_host'): assert f['exec_host'].startswith('node05/')
    state[item['job_id']]=dict(case=item['case'],state=f['job_state'],dependency=f.get('depend'),exec_host=f.get('exec_host'))
    (audit/('after_'+item['job_id']+'.txt')).write_text(text)
report=dict(time=datetime.datetime.now().isoformat(),cancelled=['7464.master','7465.master','7466.master'],submitted=submitted,active_chain=chain,states=state,rb_source_config_and_pbs_hashes_unchanged=True,preflight_passed=True)
(root/'replacement_complete.json').write_text(json.dumps(report,indent=2)+'\n')
(root/'submission.finished').write_text(datetime.datetime.now().isoformat()+'\n')
print(json.dumps(report,indent=2))
