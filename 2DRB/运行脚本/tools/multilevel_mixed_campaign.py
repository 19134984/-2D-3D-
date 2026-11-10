"""Freeze, GPU-check and submit the authorized Ra9 mixed [4,2] experiment."""
from pathlib import Path
import datetime,hashlib,io,json,os,re,shlex,subprocess,sys,tarfile

ROOT=Path(__file__).resolve().parents[2]
NAME='Multilevel_mixed42_Ra1e9_N4096_20261030_v2'
OUT=ROOT/'运行脚本/run_cases'/NAME
REMOTE='/data2/XLLi/Multiblock/'+NAME
SOLVER=ROOT/'多块网格/2DRBOpenaccMultilevelMixedRatio.F90'
POST=ROOT/'多块网格/2DRBOpenaccMultilevelMixedRatioPostprocess.F90'
SSH=['C:/msys64/usr/bin/sshpass.exe','-e','C:/msys64/usr/bin/ssh.exe','-o','ConnectTimeout=15','xlli@10.169.174.118']
def remote(cmd,data=None,check=True):
    assert os.environ.get('SSHPASS'),'Temporary session credential required'
    r=subprocess.run(SSH+[cmd],input=data,capture_output=True)
    if check and r.returncode:raise RuntimeError(r.stdout.decode(errors='replace')+r.stderr.decode(errors='replace'))
    return r
def write(p,t):p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(t if isinstance(t,bytes) else t.encode('utf-8'))
def save(p,v):write(p,json.dumps(v,indent=2,ensure_ascii=False)+'\n')
def sha(data):return hashlib.sha256(data).hexdigest()
def assign(s,name,value):
    s,n=re.subn(r'(\b'+name+r'\s*=\s*)[\d.]+(?:[dDeE][+-]?\d+)?',lambda m:m[1]+str(value),s,count=1)
    assert n==1,name
    return s
def probe(label,allow_shared=False):
    r=remote('''date -Is
qstat -n
pbsnodes node05 node07
df -h /data2/XLLi
for n in node05 node07; do
 echo GPU_NODE=$n
 ssh -o BatchMode=yes -o ConnectTimeout=10 "$n" 'nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.free --format=csv; nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv; sleep 1; nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.free --format=csv'
done
''')
    text=r.stdout.decode();nodes={}
    for node in ['node05','node07']:
        section=text.split('GPU_NODE='+node+'\n')[1].split('GPU_NODE=')[0]
        a=re.search(r'^0, Tesla P100.*?, (\d+) %, (\d+) MiB, (\d+) MiB',section,re.M)
        b=re.search(r'^(\d+) %, (\d+) MiB, (\d+) MiB',section,re.M)
        assert a and b,(node,section)
        nodes[node]=dict(utilization_pct=int(a[1]),second_utilization_pct=int(b[1]),used_MiB=int(a[2]),free_MiB=int(a[3]))
    record=dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),nodes=nodes,output=text,selected='node05',allow_shared=allow_shared,reason='User reserves node07; node05 has enough free GPU memory. Shared use may slow both solvers; other jobs are not changed.' if allow_shared else 'User reserves node07; use node05 only after live PBS/GPU checks')
    save(OUT/'node_checks'/f'{label}.json',record)
    assert nodes['node05']['free_MiB']>8192,'Insufficient free GPU memory'
    if not allow_shared:
        assert max(nodes['node05']['utilization_pct'],nodes['node05']['second_utilization_pct'])<20,'node05 currently busy; inspect before submission'
    print(json.dumps(dict(label=label,nodes=nodes)),flush=True)
    return record
def freeze_case(name,s,master,post,diagnostic=False,steps=None):
    folder=OUT/name
    assert not (folder/'job.id').exists() and not (folder/'results/run.started').exists()
    folder.mkdir(parents=True,exist_ok=True)
    for child in ['source','build','results']:(folder/child).mkdir(exist_ok=True)
    data=s.encode('utf-8');hash_=sha(data);dest=REMOTE+'/'+name
    write(folder/'source/solver.F90',data);write(folder/'source/local_master.F90',master);write(folder/'source/postprocess.F90',post)
    old=ROOT/'运行脚本/run_cases/Multilevel_original_L3_N2048_20261029/unsteady_Ra1e9_N2048_L3'
    prior=json.loads((old/'manifest.json').read_text(encoding='utf-8'))
    run=(old/'run.sh').read_text(encoding='utf-8').replace(prior['remote'],dest).replace(prior['case_source_sha256'],hash_)
    if diagnostic:
        run=run[:run.index('cd "$CASE/build"\nset +e\nnvfortran -cpp -O2 -Mextend')]
    write(folder/'run.sh',run)
    meta=dict(case=name,remote=dest,node='node05',master_sha256=sha(master),case_source_sha256=hash_,postprocess_sha256=sha(post),diagnostic=diagnostic,expected_fine_steps=steps)
    save(folder/'manifest.json',meta)
    return meta
def pbs(name,jobname,commands,walltime):
    folder=OUT/name;(folder/'results').mkdir(parents=True,exist_ok=True)
    dest=REMOTE+'/'+name
    write(folder/'run.pbs',f'''#!/bin/bash
#PBS -N {jobname}
#PBS -q batch
#PBS -l nodes=node05:ppn=1
#PBS -l walltime={walltime}
#PBS -o {dest}/results/pbs.stdout
#PBS -e {dest}/results/pbs.stderr
set -euo pipefail
test "$(hostname -s)" = node05
{commands}
''')
def collect(name):
    names=['compile.status','compile.log','compiler_version.txt','runtime_source.sha256','job_identity.txt','solver.stdout','run.status','wrapper.status','validation.status','SimulationSettings2DOpenaccMultilevelMixedRatio.txt','NuRe_2DOpenaccMultiblock.dat','FinalNuRe_2DOpenaccMultiblock.dat','postprocess.status','pbs.stderr']
    r=remote('cd '+shlex.quote(REMOTE+'/'+name+'/results')+'; for f in '+' '.join(names)+'; do test ! -f "$f" || printf "%s\\0" "$f"; done | tar --null -T - -cf -')
    with tarfile.open(fileobj=io.BytesIO(r.stdout),mode='r:') as archive:
        for item in archive:
            assert item.isfile() and item.name in names
            write(OUT/name/'results'/item.name,archive.extractfile(item).read())
def submit(name,label,allow_shared=False,dependency=None):
    assert not (OUT/name/'job.id').exists(),'Already submitted; inspect recorded job'
    probe(label,allow_shared)
    c=REMOTE+'/'+name
    options=''
    if dependency:
        assert re.fullmatch(r'\d+\.master',dependency),dependency
        options='-W depend=afterany:'+dependency+' '
    result=remote('cd '+shlex.quote(c)+'; test ! -s job.id || exit 91; job=$(qsub '+options+'run.pbs) || exit $?; printf "%s\\n" "$job" > job.id; printf "%s\\n" "$job"',check=False)
    job=result.stdout.decode().strip()
    if not re.fullmatch(r'\d+\.master',job):job=remote('cat '+shlex.quote(c+'/job.id')).stdout.decode().strip()
    assert re.fullmatch(r'\d+\.master',job),result.stderr.decode()
    write(OUT/name/'job.id',job+'\n')
    q=remote('qstat -f '+job).stdout.decode();write(OUT/name/'qstat_submission.txt',q)
    assert re.search(r'Resource_List\.nodes\s*=\s*node05:ppn=1',q),q
    host=re.search(r'exec_host\s*=\s*(\S+)',q)
    if host:assert host[1].split('/')[0]=='node05',q
    f=OUT/'submission.json';records=json.loads(f.read_text(encoding='utf-8')) if f.exists() else []
    records.append(dict(case=name,job=job,requested_node='node05',actual_node=host[1] if host else None,remote=c,state=re.search(r'job_state\s*=\s*(\S+)',q)[1],shared_gpu_allowed=allow_shared,dependency=dependency))
    save(f,records);remote('cat > '+shlex.quote(REMOTE+'/submission.json'),f.read_bytes())
    print(json.dumps(records[-1]),flush=True)

def main():
    action=sys.argv[1]
    if action=='prepare':
        assert not (OUT/'manifest.json').exists() and not list(OUT.rglob('job.id')),'Never replace frozen or submitted experiment'
        OUT.mkdir(exist_ok=True)
        master=SOLVER.read_bytes();post=POST.read_bytes().replace(b'\r\n',b'\n');s=master.decode('utf-8-sig').replace('\r\n','\n')
        assert 'interfaceRefineRatio(max(1,nGridLevels-1))=[4,2]' in s
        production=freeze_case('Ra1e9_N4096_ratio4_2',s,master,post)
        diagnostics=[]
        for name,n,ra,fine,inner,steps in [('small',256,6,32,16,64),('actual',4096,9,256,64,64)]:
            d=s
            for key,value in dict(nx=n,ny=n,Rayleigh=f'1.0d{ra}',finestLayerCellsLeftRight=fine,finestLayerCellsBottomTop=fine,innerLayerCellsLeftRight=inner,innerLayerCellsBottomTop=inner,outputSnapshotFile=0,outputPltFile=0,outputReloadFile=0).items():d=assign(d,key,value)
            d=d.replace('do while(reloadItcOffset+itc<itc_max)',f'do while(reloadItcOffset+itc<min(itc_max,{steps}))')
            diagnostics.append(freeze_case('preflight/'+name,d,master,post,True,steps))
        pbs('preflight','ML42_GPUcheck','export NVCOMPILER_ACC_SYNCHRONOUS=1\n'+'\n'.join('bash '+d['remote']+'/run.sh' for d in diagnostics),'00:30:00')
        pbs(production['case'],'ML42_R9_4096','bash '+production['remote']+'/run.sh','240:00:00')
        m=dict(name=NAME,remote=REMOTE,source=str(SOLVER),master_sha256=sha(master),cases=[production],diagnostics=diagnostics,grids=[4096,1024,512],interface_ratios=[4,2],cumulative_step=[1,4,8],finest_band_cells=256,middle_band_local_cells=64,ra=1e9,pr=.7,mach=.1,start='from_zero',target_tff=1000,statistics_tff=[500,750,1000],walltime='240:00:00',node='node05',regularization=False,numerical_formulas_changed=False,field_output_compatibility='same text file names and columns; new snapshot version with explicit interface ratios')
        save(OUT/'manifest.json',m)
        with tarfile.open(OUT/'cases.tar','w') as archive:
            for p in sorted(OUT.iterdir()):
                if p.name!='cases.tar':archive.add(p,arcname=NAME+'/'+p.name)
        print(json.dumps(m),flush=True)
    elif action=='upload':
        m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
        remote('test ! -e '+shlex.quote(REMOTE)+' && tar -xf - -C /data2/XLLi/Multiblock',(OUT/'cases.tar').read_bytes())
        commands=['set -eu']
        for c in m['cases']+m['diagnostics']:
            d=c['remote'];commands += [f'dos2unix {d}/run.sh >/dev/null 2>&1',f'bash -n {d}/run.sh',f'echo "{c["case_source_sha256"]}  {d}/source/solver.F90" | sha256sum -c -']
        for name in ['preflight',m['cases'][0]['case']]:commands += [f'dos2unix {REMOTE}/{name}/run.pbs >/dev/null 2>&1',f'bash -n {REMOTE}/{name}/run.pbs']
        result=remote('bash -s',('\n'.join(commands)+'\n').encode());write(OUT/'upload_verification.txt',result.stdout);print(result.stdout.decode())
    elif action=='preflight':submit('preflight','before_gpu_preflight',allow_shared=True)
    elif action=='status':
        m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
        for name in ['preflight']+[c['case'] for c in m['cases']]:
            if (OUT/name/'job.id').exists():
                job=(OUT/name/'job.id').read_text().strip();r=remote('qstat -f '+job,check=False);write(OUT/name/'qstat_latest.txt',r.stdout+r.stderr);print(r.stdout.decode() or r.stderr.decode())
        for c in m['diagnostics']+m['cases']:
            collect(c['case']);f=OUT/c['case']/'results'
            print('CASE='+c['case'])
            for name in ['compile.status','run.status','wrapper.status','solver.stdout']:
                if (f/name).exists():print(name,'\n'+'\n'.join((f/name).read_text(errors='replace').splitlines()[-5:]))
    elif action=='submit-production':
        m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
        for c in m['diagnostics']:
            folder=c['remote']+'/results'
            remote(f'test "$(cat {folder}/wrapper.status)" = 0 && test "$(cat {folder}/run.status)" = 0 && grep -q "Completed total fine steps: {c["expected_fine_steps"]};" {folder}/solver.stdout')
        proof=json.loads((OUT/'gpu_preflight_validation.json').read_text(encoding='utf-8'))
        assert proof['status']=='passed','GPU preflight including CPU comparison not yet accepted'
        assert sha(SOLVER.read_bytes())==m['master_sha256'],'Local master changed after freezing; prepare a fresh campaign'
        submit(m['cases'][0]['case'],'before_Ra1e9_production',allow_shared='--share-node05' in sys.argv,dependency=next((v.split('=',1)[1] for v in sys.argv if v.startswith('--after=')),None))
    elif action=='probe':probe('initial_'+datetime.datetime.now().strftime('%H%M%S'))
    else:raise SystemExit('Unknown action')

if __name__=='__main__':
    main()
