"""Freeze and submit the requested original-Luo Ra9 4096/2048/1024 case."""
from pathlib import Path
import datetime, io, json, re, sys, tarfile
import multilevel_mixed_campaign as c

NAME='Multilevel_ratio22_Ra1e9_N4096_20261031'
c.NAME=NAME
c.OUT=c.ROOT/'运行脚本/run_cases'/NAME
c.REMOTE='/data2/XLLi/Multiblock/'+NAME
CASE='Ra1e9_N4096_ratio2_2'

def source_for_case(master):
    s=master.decode('utf-8-sig').replace('\r\n','\n')
    assert 'interfaceRefineRatio(max(1,nGridLevels-1))=[4,2]' in s
    s=s.replace('interfaceRefineRatio(max(1,nGridLevels-1))=[4,2]',
                'interfaceRefineRatio(max(1,nGridLevels-1))=[2,2]')
    for name in ['innerLayerCellsLeftRight','innerLayerCellsBottomTop']:
        s=c.assign(s,name,128)
    s=s.replace('默认4倍再2倍','本算例2倍再2倍')
    s=s.replace('壁面4096级 -> 中间1024级 -> 中心512级','壁面4096级 -> 中间2048级 -> 中心1024级')
    s=s.replace('累计格距/时间步默认为1、4、8','累计格距/时间步本算例为1、2、4')
    s=s.replace('默认第2级64格乘本级格距4，对应256个最细格距','本算例第2级128格乘本级格距2，对应256个最细格距')
    s=re.sub(r'^    ! dt与四个松弛率均取当前来源级；.*$',
             '    ! dt与四个松弛率均取当前来源级；本算例[2,2]时dt=1、2、4。',s,flags=re.M)
    s=re.sub(r'^    ! 全2比例的第2级dt=2示例中，.*$',
             '    ! 本算例第2级dt=2，括号内为m_a-meq_a+F_m,a。',s,flags=re.M)
    s=re.sub(r'^    ! 4\. 用接收级的松弛率和时间步恢复完整矩；.*$',
             '    ! 4. 用接收级的松弛率和时间步恢复完整矩；本算例第2级dt=2，第3级dt=4。',s,flags=re.M)
    return s

def prepare():
    assert not c.OUT.exists(),'Never overwrite an existing campaign'
    c.OUT.mkdir(parents=True)
    master=c.SOLVER.read_bytes(); post=c.POST.read_bytes().replace(b'\r\n',b'\n')
    s=source_for_case(master)
    case=c.freeze_case(CASE,s,master,post)
    small=s
    for key,value in dict(nx=256,ny=256,Rayleigh='1.0d6',finestLayerCellsLeftRight=32,
                          finestLayerCellsBottomTop=32,innerLayerCellsLeftRight=32,
                          innerLayerCellsBottomTop=32,outputSnapshotFile=0,
                          outputPltFile=0,outputReloadFile=0).items():
        small=c.assign(small,key,value)
    small=small.replace('do while(reloadItcOffset+itc<itc_max)',
                        'do while(reloadItcOffset+itc<min(itc_max,64))')
    diagnostic=c.freeze_case('local_smoke',small,master,post,True,64)
    c.pbs(CASE,'ML22_R9_4096','bash '+case['remote']+'/run.sh','240:00:00')
    m=dict(name=NAME,remote=c.REMOTE,source=str(c.SOLVER),master_sha256=c.sha(master),
           cases=[case],diagnostics=[diagnostic],grids=[4096,2048,1024],interface_ratios=[2,2],
           cumulative_step=[1,2,4],finest_band_cells=256,middle_band_local_cells=128,
           band_physical_thickness_in_finest_units=[256,256],ra=1e9,pr=.7,mach=.1,
           start='from_zero',target_tff=1000,statistics_tff=[500,750,1000],
           walltime='240:00:00',node='node05',regularization=False,
           numerical_formulas_changed=False,comparison='Same physical ring widths as failed [4,2] case; only interface ratios and middle local cell count change')
    c.save(c.OUT/'manifest.json',m)
    with tarfile.open(c.OUT/'cases.tar','w') as ar:
        for p in sorted(c.OUT.iterdir()):
            if p.name!='cases.tar': ar.add(p,arcname=NAME+'/'+p.name)
    print(json.dumps(m),flush=True)

def upload():
    m=json.loads((c.OUT/'manifest.json').read_text(encoding='utf-8'))
    c.remote('test ! -e '+c.shlex.quote(c.REMOTE)+' && tar -xf - -C /data2/XLLi/Multiblock',(c.OUT/'cases.tar').read_bytes())
    commands=['set -eu']
    for case in m['cases']+m['diagnostics']:
        folder=case['remote']
        commands += [f'dos2unix {folder}/run.sh >/dev/null 2>&1',f'bash -n {folder}/run.sh',
                     f'echo "{case["case_source_sha256"]}  {folder}/source/solver.F90" | sha256sum -c -']
    commands += [f'dos2unix {c.REMOTE}/{CASE}/run.pbs >/dev/null 2>&1',f'bash -n {c.REMOTE}/{CASE}/run.pbs']
    result=c.remote('bash -s',('\n'.join(commands)+'\n').encode())
    c.write(c.OUT/'upload_verification.txt',result.stdout)
    print(result.stdout.decode())

def submit():
    m=json.loads((c.OUT/'manifest.json').read_text(encoding='utf-8'))
    assert c.sha(c.SOLVER.read_bytes())==m['master_sha256'],'Local source changed after freeze'
    frozen=c.OUT/CASE/'source/solver.F90'
    assert frozen.read_text(encoding='utf-8')==source_for_case(c.SOLVER.read_bytes())
    proof=json.loads((c.OUT/'local_smoke_validation.json').read_text(encoding='utf-8'))
    assert proof['status']=='passed' and proof['steps']==64
    assert proof['production_source_sha256']==c.sha(frozen.read_bytes())
    c.submit(CASE,'before_Ra1e9_ratio22',allow_shared=True)

def status():
    job=(c.OUT/CASE/'job.id').read_text().strip()
    q=c.remote('qstat -f '+job,check=False)
    c.write(c.OUT/CASE/'qstat_latest.txt',q.stdout+q.stderr)
    print(q.stdout.decode() or q.stderr.decode())
    c.collect(CASE)
    folder=c.OUT/CASE/'results'
    for name in ['compile.status','run.status','wrapper.status','solver.stdout']:
        if (folder/name).exists():print(name,'\n'+'\n'.join((folder/name).read_text(errors='replace').splitlines()[-8:]))

if __name__=='__main__':
    action=sys.argv[1]
    if action=='prepare':prepare()
    elif action=='upload':upload()
    elif action=='submit':submit()
    elif action=='status':status()
    elif action=='probe':c.probe('initial_'+datetime.datetime.now().strftime('%H%M%S'),allow_shared=True)
    else:raise SystemExit('Unknown action')
