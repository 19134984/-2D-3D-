"""Consolidate original observations; keep failed endpoints out of valid averages."""
from pathlib import Path
import csv, datetime, hashlib, json, math, re, shutil
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/multilevel_data_20261102'
groups=json.loads((OUT/'registry.json').read_text(encoding='utf-8'))
download=json.loads((OUT/'download_manifest.json').read_text(encoding='utf-8'))
for item in download:
    p=OUT/item['local'];data=p.read_bytes()
    assert len(data)==item['size'] and hashlib.sha256(data).hexdigest()==item['sha256']
cancel=ROOT/'运行脚本/run_cases/Multilevel_original_L3_20261027_v3/cancelled_and_replaced_20261029.json'
shutil.copy2(cancel,OUT/'provenance'/cancel.name)
cancelled=set(json.loads(cancel.read_text(encoding='utf-8'))['cancelled_jobs'])

def text(p):return p.read_text(encoding='utf-8',errors='replace') if p.exists() else ''
def matrix(p,n):
    rows=[]
    for line in text(p).splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            try:values=[float(x.replace('D','E')) for x in line.split()]
            except ValueError:continue
            assert len(values)==n,(p,len(values));rows.append(values)
    return np.array(rows,dtype=float).reshape(-1,n)
def value(source,key):
    m=re.search(r'\b'+key+r'\s*=\s*([\d.]+(?:[dDeE][+-]?\d+)?)',source)
    return float(m[1].replace('d','e').replace('D','E')) if m else None
def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csvfile(name,headers,rows):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(rows)
def time_average(h,a,b):
    t=h[:,2];assert t[0]<=a<b<=t[-1]
    y=h[:,3:8].copy();y[:,1]**=2
    tx=np.r_[a,t[(t>a)&(t<b)],b]
    yy=np.column_stack([np.interp(tx,t,y[:,i]) for i in range(5)])
    avg=np.trapezoid(yy,tx,axis=0)/(b-a);avg[1]=math.sqrt(avg[1]);return avg

summary=[];endpoints=[];history=[];convergence=[];statistics=[];parameters=[];diagnostics=[];checks=[]
metrics=['t_ff','NuVolAvg','ReVolRMS','Nu_hot','Nu_cold','Nu_middle','mass','meanT','Tmin','Tmax','rhoMin','rhoMax']
for group in sorted(groups,key=lambda g:int(g['job'].split('.')[0])):
    folder=OUT/'raw'/group['campaign']/group['case'];results=folder/'results'
    if group['kind']=='preflight':
        for d in group['manifest'].get('diagnostics',[]):
            rel=d['remote'].split(group['remote']+'/',1)[-1]
            p=folder/rel;rf=p/'results';source=text(p/'source/solver.F90');f=matrix(rf/'FinalNuRe_2DOpenaccMultiblock.dat',12)
            log=text(rf/'solver.stdout');complete=re.search(r'Completed total fine steps:\s*(\d+)',log)
            compile_status=text(rf/'compile.status').strip();run_status=text(rf/'run.status').strip()
            status='短测完成' if run_status=='0' else ('编译失败' if compile_status and compile_status!='0' else '前项失败，未运行')
            item=dict(job=group['job'],case=rel,status=status,compile_status=compile_status or None,run_status=run_status or None,
                      nx=value(source,'nx'),ra=value(source,'Rayleigh'),steps=int(complete[1]) if complete else None,
                      final=f[-1].tolist() if len(f) else None,source=str(p/'source/solver.F90'))
            diagnostics.append(item)
        continue
    source=text(folder/'source/solver.F90');log=text(results/'solver.stdout')
    h=matrix(results/'NuRe_2DOpenaccMultiblock.dat',14);conv=matrix(results/'Convergence_2DOpenaccMultiblock.dat',4)
    final=matrix(results/'FinalNuRe_2DOpenaccMultiblock.dat',12)
    assert not len(h) or np.isfinite(h).all();assert not len(final) or np.isfinite(final).all()
    nx=int(value(source,'nx'));levels=int(value(source,'nGridLevels'))
    ratios=re.search(r'interfaceRefineRatio\([^\n]*?\)\s*=\s*\[([^\]]+)\]',source)
    ratios=[int(v) for v in ratios[1].split(',')] if ratios else [int(value(source,'levelRefineRatio'))]*(levels-1)
    cumulative=np.r_[1,np.cumprod(ratios)].astype(int).tolist();grids=[nx//d for d in cumulative]
    steady=bool(re.search(r'^#define steadyFlow\s*$',source,re.M));mode='侧壁稳态' if steady else 'RB非稳态'
    ra=value(source,'Rayleigh');pr=value(source,'Prandtl');mach=value(source,'Mach')
    settings=next(results.glob('SimulationSettings*.txt'),None) if results.exists() else None
    setting=text(settings) if settings else ''
    time_match=re.search(r'Fine viscosity, diffusivity, timeUnit:.*\s(\S+)\n',setting)
    timeunit=float(time_match[1]) if time_match else None
    failure=re.search(r'Invalid owned cell: step,block,level,i,j=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',log)
    completed=re.search(r'Completed total fine steps:\s*(\d+)',log)
    run_status=text(results/'run.status').strip();final_text=text(results/'FinalNuRe_2DOpenaccMultiblock.dat')
    invalid_step=int(failure[1]) if failure else None
    endt=(invalid_step/timeunit if failure else (float(final[-1,0]) if len(final) else (float(h[-1,2]) if len(h) else None)))
    if group['job'] in cancelled:status='已取消（有部分数据）' if len(h) else '已取消（未启动）'
    elif failure:status='数值异常退出'
    elif steady and '# converged = T' in final_text:status='稳态判据通过'
    elif completed and run_status=='0':status='完成目标时长'
    else:status='状态待核对'
    error_match=re.search(r'# errorU errorT =\s*(\S+)\s+(\S+)',final_text)
    errors=[float(error_match[1]),float(error_match[2])] if error_match else [None,None]
    report=dict(job=group['job'],mode=mode,Ra=ra,Pr=pr,Mach=mach,grids=grids,ratios=ratios,status=status,
                end_tff=endt,last_sample_tff=float(h[-1,2]) if len(h) else None,samples=len(h),
                primary_basis='稳态停止终点' if steady and len(final) else '无有效长期统计',Nu=None,Re=None,
                errorU=errors[0],errorT=errors[1],half_Nu=None,half_Re=None,half_check='不适用',
                start=text(results/'solver.started').strip() or None,end=text(results/'run.finished').strip() or text(results/'wrapper.finished').strip() or None,
                run_status=run_status or None,failed_fine_step=invalid_step,invalid_cell=[int(x) for x in failure.groups()[1:]] if failure else None,
                campaign=group['campaign'],case=group['case'],remote_results=group['remote']+'/results/',
                note='异常退出前瞬时值不作长期统计' if failure else ('用户取消，未完成统计区间' if group['job'] in cancelled else ''),
                source_relative=str((folder/'source/solver.F90').relative_to(OUT)))
    if steady and len(final):report['Nu'],report['Re']=final[-1,1:3].tolist()
    if not steady and completed and run_status=='0' and len(h) and h[-1,2]>=1000:
        avgs=[time_average(h,a,b) for a,b in [(500,1000),(500,750),(750,1000)]]
        diff=np.abs(avgs[2]-avgs[1])/np.maximum(np.abs(avgs[0]),1e-30)
        for (a,b,label),avg in zip([(500,1000,'whole'),(500,750,'first'),(750,1000,'last')],avgs):
            statistics.append(dict(job=group['job'],window=label,start=a,end=b,metrics=avg.tolist(),half_diff=diff.tolist() if label=='whole' else None))
        stattext=text(results/'NuReStatistics_2DOpenaccMultiblock.dat')
        maximum=0.
        for label,avg in zip(['whole','first','last'],avgs):
            found=re.search(r'^'+label+r'\s*:\s*(.*)$',stattext,re.M);assert found,(group['job'],label)
            old=np.array([float(v) for v in found[1].split()]);relative=np.max(np.abs(old-avg)/np.maximum(1,np.abs(old)));maximum=max(maximum,float(relative))
        assert maximum<1e-10,(group['job'],maximum)
        report.update(primary_basis='500–1000 TFF 时间统计',Nu=float(avgs[0][0]),Re=float(avgs[0][1]),half_Nu=float(diff[0]),half_Re=float(diff[1]),half_check='Nu/Re半窗均通过1%' if max(diff[:2])<.01 else '半窗未通过1%')
        if report['half_check']=='半窗未通过1%':report['note']='完成时长；Re前后半窗差超过1%，统计尚未通过'
        checks.append(dict(job=group['job'],independent_history_integration_vs_postprocess_max_relative_difference=maximum))
    last=(final[-1].tolist() if len(final) else ([float(h[-1,2])]+h[-1,3:].tolist() if len(h) else None))
    if last:
        endpoints.append(dict(job=group['job'],basis='稳态停止终点' if steady and len(final) else ('完成时瞬时场' if len(final) else '最后保存采样'),valid_long_term_result=bool(steady and len(final)),values=last))
    for row in h:history.append([group['job']]+row.tolist())
    for row in conv:convergence.append([group['job']]+row.tolist())
    levelparams=[]
    for match in re.finditer(r'Level (\d+) dx,Snu,Sq,Qk,Qnu,area,tau:\s*(.*)',setting):
        vals=[float(x) for x in match[2].split()];levelparams.append(dict(level=int(match[1]),dx=vals[0],Snu=vals[1],Sq=vals[2],Qk=vals[3],Qnu=vals[4],area=vals[5],tau=vals[6]))
    runtime=text(results/'runtime_source.sha256').split();localsha=hashlib.sha256((folder/'source/solver.F90').read_bytes()).hexdigest()
    if runtime:assert runtime[0]==localsha,(group['job'],runtime[0],localsha)
    wall=re.search(r'Wall seconds:\s*([\d.]+)',log)
    physics_match=re.search(r'^Physics:\s*(.*)$',setting,re.M)
    physics=[float(x) for x in physics_match[1].split()] if physics_match else None
    diffusion_match=re.search(r'^Fine viscosity, diffusivity, timeUnit:\s*(.*)$',setting,re.M)
    fine_coeff=[float(x) for x in diffusion_match[1].split()] if diffusion_match else None
    duration=(datetime.datetime.fromisoformat(report['end'])-datetime.datetime.fromisoformat(report['start'])).total_seconds() if report['start'] and report['end'] else None
    parameters.append(dict(job=group['job'],mode=mode,Ra=ra,Pr=pr,Mach=mach,nx=nx,levels=levels,grids=grids,ratios=ratios,
                           fine_band=int(value(source,'finestLayerCellsLeftRight')),middle_band=int(value(source,'innerLayerCellsLeftRight')),
                           physical_band_widths=[int(value(source,'finestLayerCellsLeftRight')),int(value(source,'innerLayerCellsLeftRight'))*cumulative[1]],
                           timeUnit=timeunit,paraA=physics[5] if physics else None,
                           fine_viscosity=fine_coeff[0] if fine_coeff else None,fine_diffusivity=fine_coeff[1] if fine_coeff else None,
                           level_parameters=levelparams,source_sha256=localsha,runtime_sha256=runtime[0] if runtime else None,
                           solver_wall_seconds=float(wall[1]) if wall else duration,started=report['start'],ended=report['end'],
                           source_relative=report['source_relative'],remote_results=report['remote_results'],
                           stored_nodes_fraction=float(re.search(r'Stored nodes / full fine nodes =\s*(\S+)',setting)[1]) if 'Stored nodes / full fine nodes' in setting else None))
    summary.append(report)

dataset=dict(checked_at=json.loads((OUT/'server_inventory.json').read_text())['checked_at'],summary=summary,endpoints=endpoints,
             statistics=statistics,parameters=parameters,diagnostics=diagnostics,history=history,convergence=convergence,
             validation=checks,quantitative_file_count=len(download),large_file_count=len(json.loads((OUT/'large_field_file_index.json').read_text())),
             protocol=dict(Nu='time average by trapezoidal integration with interpolated exact window endpoints',
                           Re='sqrt(time average of ReVolRMS squared)',half_difference='abs(last-first)/abs(whole)',
                           steady='FinalNuRe at actual converged stopping step, not last periodic history sample',
                           scope='No new solver runs. No benchmark agreement or grid independence claim.'))
save('dataset.json',dataset)
headers=['作业','类型','Ra','Pr','Ma','网格','接口比例','状态','结束或失败TFF','Nu统计口径','Nu','Re','errorU','errorT','Nu半窗相对差','Re半窗相对差','半窗检查','备注','服务器结果目录']
rows=[[x['job'],x['mode'],x['Ra'],x['Pr'],x['Mach'],'/'.join(map(str,x['grids'])),'/'.join(map(str,x['ratios'])),x['status'],x['end_tff'],x['primary_basis'],x['Nu'],x['Re'],x['errorU'],x['errorT'],x['half_Nu'],x['half_Re'],x['half_check'],x['note'],x['remote_results']] for x in summary]
csvfile('算例汇总.csv',headers,rows)
csvfile('NuRe全部历史.csv',['job','sample','itc_total']+metrics,history)
csvfile('稳态收敛全部历史.csv',['job','itc_total','t_ff','errorU','errorT'],convergence)
csvfile('终点与最后采样.csv',['job','basis','steady_converged']+metrics,[[x['job'],x['basis'],x['valid_long_term_result']]+x['values'] for x in endpoints])
csvfile('非稳态时间统计.csv',['job','window','t_start','t_end','NuVolAvg','Re_rms_space_time','Nu_hot','Nu_cold','Nu_middle'],[[x['job'],x['window'],x['start'],x['end']]+x['metrics'] for x in statistics])
large=json.loads((OUT/'large_field_file_index.json').read_text(encoding='utf-8'))
csvfile('全量场与重启文件索引.csv',['job','server_path','relative_path','bytes','GiB'],[[x['job'],x['path'],x['relative'],x['size'],x['size']/2**30] for x in large])
for x in summary:print(x['job'],x['status'],'/'.join(map(str,x['grids'])),x['end_tff'],x['Nu'],x['Re'],x['half_check'])
print('History rows',len(history),'Convergence rows',len(convergence),'diagnostic attempts',len(diagnostics))
