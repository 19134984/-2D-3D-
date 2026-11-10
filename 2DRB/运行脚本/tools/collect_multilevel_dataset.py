"""Read-only collection of every submitted multilevel-grid run in this campaign."""
from pathlib import Path
import datetime, io, json, re, shutil, subprocess, sys, tarfile, hashlib
import multilevel_mixed_campaign as transport

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/multilevel_data_20261102'
RUNS=ROOT/'运行脚本/run_cases'

def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def registry():
    groups=[]
    for folder in sorted(RUNS.glob('Multilevel*')):
        if not (folder/'manifest.json').exists():continue
        m=json.loads((folder/'manifest.json').read_text(encoding='utf-8-sig'))
        sub=json.loads((folder/'submission.json').read_text(encoding='utf-8-sig')) if (folder/'submission.json').exists() else []
        for entry in sub:
            case=entry['case']; job=entry.get('job_id',entry.get('job'))
            info=next((v for v in m.get('cases',[]) if v['case']==case),None)
            group=dict(campaign=folder.name,case=case,job=job,remote=entry['remote'],
                       kind='production' if info else 'preflight',manifest=m,case_manifest=info,
                       local_campaign=str(folder),requested_node=entry.get('requested_node'),
                       actual_node=entry.get('actual_node'))
            groups.append(group)
    return groups

def inventory():
    groups=registry();OUT.mkdir(parents=True,exist_ok=True);save(OUT/'registry.json',groups)
    payload=[dict(remote=g['remote'],job=g['job']) for g in groups]
    script='''import os,json,subprocess,datetime
groups=GROUPS
files=[]
for group in groups:
 root=group['remote']
 for folder,dirs,names in os.walk(root):
  if 'build' in dirs:dirs.remove('build')
  for name in names:
   p=os.path.join(folder,name)
   try:
    st=os.stat(p)
    files.append(dict(job=group['job'],group=root,path=p,relative=os.path.relpath(p,root),size=st.st_size,mtime=st.st_mtime))
   except FileNotFoundError:pass
q=subprocess.run(['qstat','-n'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
print(json.dumps(dict(checked_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=files,queue=q.stdout+q.stderr)))
'''.replace('GROUPS',repr(payload))
    result=transport.remote('python3 -',script.encode())
    info=json.loads(result.stdout);save(OUT/'server_inventory.json',info)
    for g in groups:
        f=[x for x in info['files'] if x['group']==g['remote']]
        print(g['job'],g['kind'],g['case'],'files',len(f),'MiB',round(sum(x['size'] for x in f)/2**20,1))
    print('TOTAL GiB',round(sum(x['size'] for x in info['files'])/2**30,2))

def fetch():
    info=json.loads((OUT/'server_inventory.json').read_text(encoding='utf-8'))
    groups=json.loads((OUT/'registry.json').read_text(encoding='utf-8'))
    selected=[];large=[]
    for item in info['files']:
        p=Path(item['relative']);name=p.name
        # Text results and provenance retain the numerical columns verbatim.
        # Large field/restart/snapshot data remain indexed separately.
        binary=(name.lower().endswith(('.bin','.plt','.szplt','.vtk','.vtu')) or name.startswith('reloadFile'))
        wanted=(not binary and item['size']<50*2**20)
        if wanted:selected.append(item)
        else:large.append(item)
    save(OUT/'large_field_file_index.json',large)
    names=[x['path'] for x in selected]
    cmd='tar -czf - --null -T -'
    result=transport.remote(cmd,b'\0'.join(p.encode() for p in names)+b'\0')
    raw=OUT/'raw';raw.mkdir(exist_ok=True)
    by_path={x['path']:x for x in selected}
    counts={}
    with tarfile.open(fileobj=io.BytesIO(result.stdout),mode='r:gz') as ar:
        for member in ar:
            assert member.isfile()
            key='/'+member.name.lstrip('/')
            assert key in by_path,key
            item=by_path[key];g=next(g for g in groups if g['remote']==item['group'])
            relative=Path(item['relative']);assert '..' not in relative.parts
            dest=raw/g['campaign']/g['case']/relative
            dest.parent.mkdir(parents=True,exist_ok=True)
            data=ar.extractfile(member).read();dest.write_bytes(data)
            item['local']=str(dest.relative_to(OUT));item['sha256']=hashlib.sha256(data).hexdigest()
            counts[g['job']]=counts.get(g['job'],0)+1
    save(OUT/'download_manifest.json',selected)
    for g in groups:
        folder=Path(g['local_campaign']);d=OUT/'provenance'/g['campaign'];d.mkdir(parents=True,exist_ok=True)
        for name in ['manifest.json','submission.json','parameter_verification.json','status_20261031.json','status_20261101.json','gpu_preflight_validation.json','local_smoke_validation.json']:
            if (folder/name).exists():shutil.copy2(folder/name,d/name)
    print('Downloaded',len(selected),'files;',round(sum(x['size'] for x in selected)/2**20,2),'MiB uncompressed')
    print('Large field/restart inventory:',len(large),'files;',round(sum(x['size'] for x in large)/2**30,2),'GiB')
    print(json.dumps(counts))

if __name__=='__main__':
    if sys.argv[1]=='inventory':inventory()
    elif sys.argv[1]=='fetch':fetch()
    else:raise SystemExit('Unknown action')
