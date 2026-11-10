"""Replace the two authorized 4096 jobs with fresh original-Luo 2048 cases."""
from pathlib import Path
import datetime
import hashlib
import io
import json
import os
import re
import shlex
import subprocess
import sys
import tarfile

ROOT=Path(__file__).resolve().parents[2]
NAME='Multilevel_original_L3_N2048_20261029'
OUT=ROOT/'运行脚本/run_cases'/NAME
REMOTE='/data2/XLLi/Multiblock/'+NAME
OLD=ROOT/'运行脚本/run_cases/Multilevel_original_L3_20261027_v3'
SSH=['C:/msys64/usr/bin/sshpass.exe','-e','C:/msys64/usr/bin/ssh.exe','-o','ConnectTimeout=15','xlli@10.169.174.118']
OUT.mkdir(exist_ok=True)
def remote(command,data=None,check=True):
    assert os.environ.get('SSHPASS'),'Use temporary session credential'
    result=subprocess.run(SSH+[command],input=data,capture_output=True)
    if check and result.returncode:
        raise RuntimeError(result.stdout.decode(errors='replace')+result.stderr.decode(errors='replace'))
    return result
def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes(text if isinstance(text,bytes) else text.encode('utf-8'))
def save(name,value):
    write(OUT/name,json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def probe(label):
    command='''date -Is
qstat -n
pbsnodes node05 node07
df -h /data2/XLLi
for n in node05 node07; do
 echo GPU_NODE=$n
 ssh -o BatchMode=yes -o ConnectTimeout=10 "$n" 'nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.free --format=csv; nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv; sleep 1; nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.free --format=csv'
done
'''
    r=remote(command);text=r.stdout.decode();nodes={}
    for n in ['node05','node07']:
        section=text.split('GPU_NODE='+n+'\n')[1].split('GPU_NODE=')[0]
        match=re.search(r'^0, Tesla P100.*?, (\d+) %, (\d+) MiB, (\d+) MiB',section,re.M)
        second=re.search(r'^(\d+) %, (\d+) MiB, (\d+) MiB',section,re.M)
        assert match and second,(n,section)
        nodes[n]=dict(utilization_pct=int(match[1]),used_MiB=int(match[2]),free_MiB=int(match[3]),second_utilization_pct=int(second[1]))
    record=dict(timestamp_utc=now(),nodes=nodes,output=text)
    save('node_checks/'+label+'.json',record)
    print(json.dumps(dict(label=label,nodes=nodes)),flush=True)
    return record
def parameter(text,name,value):
    result,count=re.subn(r'(\b'+name+r'\s*=\s*)[\d.]+(?:[dDeE][+-]?\d+)?',lambda m:m[1]+str(value),text,count=1)
    assert count==1,name
    return result

action=sys.argv[1]
if action=='cancel-old':
    records=[]
    # Remove the waiting successor first, preventing its automatic start.
    for exponent,job in [(9,'7428.master'),(8,'7427.master')]:
        case=f'unsteady_Ra1e{exponent}_N4096_L3'
        prior=json.loads((OLD/case/'manifest.json').read_text(encoding='utf-8'))
        r=remote('qstat -f '+job,check=False)
        before=r.stdout.decode()
        write(OUT/'old_jobs'/job/'qstat_before.txt',r.stdout+r.stderr)
        normalized=re.sub(r'\n\s+','',before)
        assert re.search(r'Job_Owner\s*=\s*xlli@',before),before
        assert 'ML3_U_R'+str(exponent) in before and prior['remote'] in normalized,before
        state=re.search(r'job_state\s*=\s*(\S+)',before)[1]
        command='date -Is; '+ 'c='+shlex.quote(prior['remote']+'/results')+'; tail -3 "$c/NuRe_2DOpenaccMultiblock.dat" 2>/dev/null; cat "$c/reloadFile2DOpenaccMultilevel-latest.meta" 2>/dev/null; true'
        write(OUT/'old_jobs'/job/'last_history_and_checkpoint.txt',remote(command).stdout)
        if state not in ['C','E']:
            deleted=remote('qdel '+job)
            write(OUT/'old_jobs'/job/'qdel.txt',deleted.stdout+deleted.stderr)
        records.append(dict(job=job,case=case,old_remote=prior['remote'],state_before=state,cancel_requested_at=now(),source_and_results_preserved=True))
        print(json.dumps(records[-1]),flush=True)
    save('cancellation.json',records)
    probe('after_cancellation')
elif action=='verify-cancelled':
    for job in ['7428.master','7427.master']:
        r=remote('qstat -f '+job,check=False)
        write(OUT/'old_jobs'/job/'qstat_after.txt',r.stdout+r.stderr)
        t=r.stdout.decode();state=re.search(r'job_state\s*=\s*(\S+)',t)
        assert not state or state[1]=='C',t
        print(json.dumps(dict(job=job,state=state[1] if state else 'removed_from_queue')),flush=True)
    record=probe('confirmed_old_jobs_stopped')
    for info in record['nodes'].values():
        assert info['used_MiB']<4096,'Check unexpected remaining GPU processes'
elif action=='prepare':
    assert not (OUT/'manifest.json').exists(),'Never overwrite frozen cases'
    observation=probe('before_preparation')
    # Prefer the idle, otherwise unallocated node and use two GPUs when available.
    candidates=['node05']  # User reserves node07 for other work in this session.
    usable=[n for n in candidates if max(observation['nodes'][n]['utilization_pct'],observation['nodes'][n]['second_utilization_pct'])<20 and observation['nodes'][n]['free_MiB']>8192]
    assert usable,'No currently idle suitable GPU'
    source=ROOT/'多块网格/2DRBOpenaccMultilevel.F90'
    post=ROOT/'多块网格/2DRBOpenaccMultilevelPostprocess.F90'
    raw=source.read_bytes();postraw=post.read_bytes()
    template=raw.decode('utf-8-sig').replace('\r\n','\n')
    assert re.search(r'^#define unsteadyFlow$',template,re.M) and re.search(r'^#define RayleighBenardCell$',template,re.M)
    sha=lambda v:hashlib.sha256(v).hexdigest()
    cases=[]
    for index,exponent in enumerate([8,9]):
        node=usable[index%len(usable)]
        name=f'unsteady_Ra1e{exponent}_N2048_L3';folder=OUT/name
        folder.mkdir(exist_ok=False)
        for child in ['source','results','build']:(folder/child).mkdir()
        changes=dict(nx=2048,ny=2048,nGridLevels=3,Rayleigh=f'1.0d{exponent}',Prandtl='0.7d0',Mach='0.1d0',finestLayerCellsLeftRight=128,finestLayerCellsBottomTop=128,innerLayerCellsLeftRight=64,innerLayerCellsBottomTop=64,loadInitField=0)
        text=template
        for key,value in changes.items():text=parameter(text,key,value)
        data=text.encode('utf-8');casehash=sha(data)
        prior=json.loads((OLD/f'unsteady_Ra1e{exponent}_N4096_L3/manifest.json').read_text(encoding='utf-8'))
        script=(OLD/prior['case']/'run.sh').read_text(encoding='utf-8')
        script=script.replace(prior['remote'],REMOTE+'/'+name).replace(prior['case_source_sha256'],casehash).replace('= node05','= '+node)
        write(folder/'source/solver.F90',data);write(folder/'source/local_master.F90',raw)
        write(folder/'source/postprocess.F90',postraw.replace(b'\r\n',b'\n'))
        write(folder/'run.sh',script)
        walltime='120:00:00'
        pbs=f'''#!/bin/bash
#PBS -N ML3_2K_R{exponent}
#PBS -q batch
#PBS -l nodes={node}:ppn=1
#PBS -l walltime={walltime}
#PBS -o {REMOTE}/{name}/results/pbs.stdout
#PBS -e {REMOTE}/{name}/results/pbs.stderr
set -euo pipefail
test "$(hostname -s)" = {node}
bash {REMOTE}/{name}/run.sh
'''
        write(folder/'run.pbs',pbs)
        meta=dict(case=name,remote=REMOTE+'/'+name,node=node,ra=10**exponent,pr=.7,mach=.1,mode='unsteady',geometry='RayleighBenardCell',grids=[2048,1024,512],nGridLevels=3,fine_band_cells=128,inner_band_local_cells=64,source_parameter_changes=changes,start='from_zero',target_tff=1000,statistics_tff=[500,750,1000],sampling_tff=.5,checkpoint_tff=100,tecplot_tff=100,walltime=walltime,regularization=False,numerical_kernels_changed=False,master_sha256=sha(raw),case_source_sha256=casehash,postprocess_master_sha256=sha(postraw),postprocess_case_sha256=sha(postraw.replace(b'\r\n',b'\n')),cancelled_predecessor='7427.master' if exponent==8 else '7428.master')
        write(folder/'manifest.json',json.dumps(meta,indent=2)+'\n');cases.append(meta)
    manifest=dict(name=NAME,remote=REMOTE,created_at=now(),source=str(source),master_sha256=sha(raw),cases=cases,node_selection_reason='User explicitly requests node05 for both cases and reserves node07; verify resources then serialize on node05',fallback_serial=True)
    save('manifest.json',manifest)
    with tarfile.open(OUT/'cases.tar','w') as archive:
        for p in sorted(OUT.iterdir()):
            if p.name!='cases.tar':archive.add(p,arcname=NAME+'/'+p.name)
    print(json.dumps(manifest),flush=True)
elif action=='pin-node05':
    assert not (OUT/'submission.json').exists(),'Only update unsubmitted cases'
    assert not (OUT/'upload_verification.txt').exists(),'Only update before upload'
    m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
    for c in m['cases']:
        assert not (OUT/c['case']/'job.id').exists()
        folder=OUT/c['case'];oldnode=c['node']
        for name in ['run.sh','run.pbs']:
            text=(folder/name).read_text(encoding='utf-8').replace('= '+oldnode,'= node05').replace('nodes='+oldnode+':','nodes=node05:')
            write(folder/name,text)
        c['node']='node05';write(folder/'manifest.json',json.dumps(c,indent=2)+'\n')
    m['node_selection_reason']='User explicitly requests node05 for both cases and reserves node07; serialize with afterany dependency'
    m['fallback_serial']=True;m['node_override_at']=now()
    save('manifest.json',m)
    with tarfile.open(OUT/'cases.tar','w') as archive:
        for p in sorted(OUT.iterdir()):
            if p.name!='cases.tar':archive.add(p,arcname=NAME+'/'+p.name)
    print('Both cases pinned to node05 before upload; solver sources unchanged')
elif action=='upload':
    m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
    remote('test ! -e '+shlex.quote(REMOTE)+' && tar -xf - -C /data2/XLLi/Multiblock',(OUT/'cases.tar').read_bytes())
    lines=['set -eu']
    for c in m['cases']:
        d=c['remote'];lines += [f'dos2unix {d}/run.sh {d}/run.pbs >/dev/null 2>&1',f'bash -n {d}/run.sh',f'bash -n {d}/run.pbs',f'echo "{c["case_source_sha256"]}  {d}/source/solver.F90" | sha256sum -c -']
    r=remote('bash -s',('\n'.join(lines)+'\n').encode());write(OUT/'upload_verification.txt',r.stdout);print(r.stdout.decode())
elif action=='submit':
    m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'));submitted=[]
    if (OUT/'submission.json').exists():submitted=json.loads((OUT/'submission.json').read_text(encoding='utf-8'))
    for c in m['cases']:
        if (OUT/c['case']/'job.id').exists():continue
        for oldjob in ['7427.master','7428.master']:
            old=remote('qstat -f '+oldjob,check=False).stdout.decode();state=re.search(r'job_state\s*=\s*(\S+)',old)
            assert not state or state[1]=='C','Old job not stopped: '+oldjob
        check=probe('before_submit_'+c['case']);info=check['nodes'][c['node']]
        same_node=[s for s in submitted if s['requested_node']==c['node']]
        assert info['free_MiB']>8192
        if not same_node:assert max(info['utilization_pct'],info['second_utilization_pct'])<20,'Selected GPU no longer idle'
        dep=' -W depend=afterany:'+same_node[-1]['job'] if same_node else ''
        command='cd '+shlex.quote(c['remote'])+'; test ! -s job.id || exit 91; job=$(qsub'+dep+' run.pbs) || exit $?; printf "%s\\n" "$job" > job.id; printf "%s\\n" "$job"'
        r=remote(command,check=False);job=r.stdout.decode().strip()
        if not re.fullmatch(r'\d+\.master',job):job=remote('cat '+shlex.quote(c['remote']+'/job.id')).stdout.decode().strip()
        assert re.fullmatch(r'\d+\.master',job),r.stderr.decode()
        write(OUT/c['case']/'job.id',job+'\n')
        q=remote('qstat -f '+job).stdout.decode();write(OUT/c['case']/'qstat_submission.txt',q)
        assert re.search(r'Resource_List\.nodes\s*=\s*'+c['node']+r':ppn=1',q),q
        host=re.search(r'exec_host\s*=\s*(\S+)',q)
        if host:assert host[1].split('/')[0]==c['node'],q
        item=dict(case=c['case'],job=job,requested_node=c['node'],actual_node=host[1] if host else None,state=re.search(r'job_state\s*=\s*(\S+)',q)[1],dependency=same_node[-1]['job'] if same_node else None,remote=c['remote'])
        submitted.append(item);save('submission.json',submitted)
        remote('cat > '+shlex.quote(REMOTE+'/submission.json'),(OUT/'submission.json').read_bytes())
        print(json.dumps(item),flush=True)
elif action=='status':
    m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
    for s in json.loads((OUT/'submission.json').read_text(encoding='utf-8')):
        r=remote('qstat -f '+s['job'],check=False);q=r.stdout.decode()
        write(OUT/s['case']/'qstat_latest.txt',r.stdout+r.stderr)
        print(q or r.stderr.decode(),flush=True)
        names=['compile.status','compile.log','runtime_source.sha256','job_identity.txt','solver.stdout','run.status','wrapper.status','validation.status','SimulationSettings2DOpenaccMultilevel.txt','NuRe_2DOpenaccMultiblock.dat','FinalNuRe_2DOpenaccMultiblock.dat','postprocess.status','NuReStatistics_2DOpenaccMultiblock.dat','postprocess.stdout','pbs.stderr','reloadFile2DOpenaccMultilevel-latest.meta','solver.started','run.finished']
        command='cd '+shlex.quote(s['remote']+'/results')+'; for f in '+' '.join(names)+'; do test ! -f "$f" || printf "%s\\0" "$f"; done | tar --null -T - -cf -'
        result=remote(command)
        with tarfile.open(fileobj=io.BytesIO(result.stdout),mode='r:') as archive:
            for member in archive:
                assert member.isfile() and member.name in names
                write(OUT/s['case']/'results'/member.name,archive.extractfile(member).read())
else:raise SystemExit('Unknown action')
