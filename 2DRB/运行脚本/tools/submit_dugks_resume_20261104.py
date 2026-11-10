"""Run on P100 master after a successful two-case GPU gate; submit once, serially."""
from pathlib import Path
import datetime,hashlib,json,re,subprocess
root=Path('/data2/XLLi/DUGKS/DUGKS_RB_resume_P100_node05_20261104')
cfg=json.loads((root/'campaign.json').read_text())
audit=root/'audit'
def run(args,check=True):
    p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    with (audit/'operations.log').open('a') as f:
        f.write(datetime.datetime.now().astimezone().isoformat()+' '+repr(args)+'\n'+p.stdout+'\nrc='+str(p.returncode)+'\n')
    if check and p.returncode: raise RuntimeError(repr(args)+' '+p.stdout)
    return p
def fields(text):
    out={}; key=None
    for line in text.splitlines():
        if line.startswith('    ') and ' = ' in line:
            key,value=line.strip().split(' = ',1); out[key]=value
        elif key and line.startswith('\t'): out[key]+=line.strip()
    return out
def nodecheck(path):
    checks=[]
    for args in [['date','-Is'],['qstat','-a'],['pbsnodes','node05','node07'],['df','-Pk',str(root)]]:
        checks.append(dict(command=args,output=run(args).stdout))
    for node in ['node05','node07']:
        cmd='date -Is; nvidia-smi --query-gpu=memory.total,memory.used,memory.free,utilization.gpu --format=csv; nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv'
        checks.append(dict(node=node,output=run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10',node,cmd]).stdout))
    path.write_text(json.dumps(dict(time=datetime.datetime.now().astimezone().isoformat(),checks=checks,reason='User explicitly requested restarting DUGKS on node05. Node05 idle at initial check, node07 runs FPC; only the two unfinished RB cases are submitted in serial.'),indent=2)+'\n')
    gpu=next(x['output'] for x in checks if x.get('node')=='node05')
    values=next(x for x in gpu.splitlines() if re.match(r'^\d+ MiB,',x))
    assert int(values.split(',')[2].strip().split()[0])>=4096
gate_dir=root/json.loads((root/'active_preflight.json').read_text())['directory'] if (root/'active_preflight.json').exists() else root/'preflight'
gate=json.loads((gate_dir/'assessment.json').read_text())
assert gate['passed'] and len(gate['cases'])==2
assert (gate_dir/'job.exit.status').read_text().strip()=='0'
assert not (root/'submission.started').exists(),'Never duplicate an existing submission'
queue=run(['qstat','-a']).stdout
assert 'DG_R8' not in queue and 'DG_R9' not in queue,'Another active DUGKS RB task may exist'
for job in ['7469.master','7470.master']:
    p=run(['qstat','-f',job],False)
    assert p.returncode!=0 or fields(p.stdout).get('job_state')=='C'
for c in cfg['cases']:
    folder=root/c['name']
    assert hashlib.sha256((folder/'source/solver.F90').read_bytes()).hexdigest()==c['source_sha256']
    smoke=json.loads((folder/'results/smoke_validation.json').read_text())
    assert smoke['passed'] and smoke['source_sha256']==c['source_sha256']
    if c['restart']:
        restart=json.loads((folder/'results/restart_preparation.json').read_text())
        assert restart['checkpoint_finite'] and restart['parameters_match'] and restart['histories_reconstructed_from_checkpoint']
        assert hashlib.sha256((folder/'results'/restart['checkpoint_name']).read_bytes()).hexdigest()==c['expected_seed_sha256']
(root/'submission.started').write_text(datetime.datetime.now().astimezone().isoformat()+'\n')
submitted=[]; previous=None
for c in cfg['cases']:
    folder=root/c['name']; nodecheck(folder/'pre_submit_node_check.json')
    run(['bash','-n',str(folder/'run.pbs')])
    args=['qsub','-h']
    if previous: args+=['-W','depend=afterany:'+previous]
    args.append(str(folder/'run.pbs'))
    job=run(args).stdout.strip(); assert re.fullmatch(r'\d+\.master',job)
    item=dict(case=c['name'],job_id=job,dependency=previous,node='node05',directory=str(folder),
              nx=c['nx'],ny=c['ny'],Ra=c['Ra'],restart=c['restart'],restart_tff=c.get('restart_tff',0),
              target_tff=c['target_tff'],walltime=c['walltime'],source_sha256=c['source_sha256'])
    submitted.append(item)
    (root/'submission.json').write_text(json.dumps(submitted,indent=2)+'\n')
    text=run(['qstat','-f',job]).stdout; f=fields(text)
    (folder/'qstat_after_submit.txt').write_text(text)
    assert f['Resource_List.nodes']=='node05:ppn=1' and f['job_state']=='H' and 'u' in f['Hold_Types']
    if previous: assert 'afterany:'+previous in f['depend']
    previous=job
for item in reversed(submitted): run(['qrls','-h','u',item['job_id']])
states={}
for item in submitted:
    text=run(['qstat','-f',item['job_id']]).stdout; f=fields(text)
    assert 'u' not in f['Hold_Types'] and f['Resource_List.nodes']=='node05:ppn=1'
    if f.get('exec_host'): assert f['exec_host'].startswith('node05/')
    states[item['job_id']]=dict(state=f['job_state'],exec_host=f.get('exec_host'),dependency=f.get('depend'))
    (audit/('after_submit_'+item['job_id']+'.txt')).write_text(text)
report=dict(submitted_at=datetime.datetime.now().astimezone().isoformat(),submitted=submitted,states=states,
            preflight_passed=True,original_data_unchanged=True,completed_cases_not_rerun=True)
(root/'submission_complete.json').write_text(json.dumps(report,indent=2)+'\n')
(root/'active_chain.json').write_text(json.dumps(submitted,indent=2)+'\n')
(root/'submission.finished').write_text(datetime.datetime.now().astimezone().isoformat()+'\n')
print(json.dumps(report,indent=2))
