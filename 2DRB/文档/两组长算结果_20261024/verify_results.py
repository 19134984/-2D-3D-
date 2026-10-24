"""对下载的两组结果运行已审阅Fortran后处理，并独立复核时间积分。"""
from pathlib import Path
import tarfile,json,math,hashlib,subprocess,re,shutil
P=Path(__file__).resolve().parent
REPO=P.parent.parent
E=P/'server_evidence';E.mkdir(exist_ok=True)
with tarfile.open(P/'two_completed_cases_20261024.tar.gz') as tf:
    for m in tf.getmembers():
        f=(E/m.name).resolve()
        assert E.resolve() in f.parents and m.isfile()
        f.parent.mkdir(parents=True,exist_ok=True)
        f.write_bytes(tf.extractfile(m).read())

source=REPO/'多块网格/2DRBOpenaccMultiblockPostprocess.F90'
B=P/'build';B.mkdir(exist_ok=True)
(B/'postprocess.F90').write_text(source.read_text(encoding='utf-8'),encoding='utf-8',newline='\n')
fc='C:/msys64/ucrt64/bin/gfortran.exe'
def command(args,d,log):
    z=subprocess.run(args,cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (d/log).write_bytes(z.stdout)
    assert z.returncode==0,z.stdout.decode('utf-8',errors='replace')
command([fc,'-cpp','-ffree-line-length-none','-fsyntax-only','postprocess.F90'],B,'syntax.log')
command([fc,'-cpp','-ffree-line-length-none','-O2','-fcheck=all','postprocess.F90','-o','postprocess.exe'],B,'compile.log')

def mean(rows,lo,hi,col,square=False):
    total=[];coverage=0.
    for a,b in zip(rows,rows[1:]):
        left=max(lo,a[2]);right=min(hi,b[2])
        if right<=left:continue
        va=a[col]**2 if square else a[col];vb=b[col]**2 if square else b[col]
        x=va+(vb-va)*(left-a[2])/(b[2]-a[2]);y=va+(vb-va)*(right-a[2])/(b[2]-a[2])
        total.append(.5*(x+y)*(right-left));coverage+=right-left
    assert abs(coverage-(hi-lo))<1e-8
    return math.fsum(total)/(hi-lo)

cases=[('Ra9_RR3_F2048','diagnosis_Ra1e9_N2048_RR3_T1000_20261022/rr3',2048,1e9,53.51,15101),
       ('Ra8_CDE_coarse05_F512','Multiblock_LBMCDE_Ra8_CoarseChi05_T1000_20261023/Ra1e8_fine0_coarse05',512,1e8,25.25,3662)]
result={'checked_date':'2026-10-24','postprocess_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'cases':{}}
for name,rel,N,Ra,refnu,refre in cases:
    base=E/rel;d=base/'results';src=base/'source/solver.F90'
    h=hashlib.sha256(src.read_bytes()).hexdigest()
    assert h in (d/'runtime_source.sha256').read_text()
    assert int((d/'run.status').read_text())==0 and int((d/'compile.status').read_text())==0
    text=(d/'solver.stdout').read_text(errors='replace')
    assert 'Invalid state:' not in text and 'ERROR STOP' not in text
    rows=[list(map(float,l.split())) for l in (d/'NuRe_2DOpenaccMultiblock.dat').read_text().splitlines() if l.strip() and not l.lstrip().startswith('#')]
    assert len(rows)==2000 and all(len(r)==14 and all(math.isfinite(x) for x in r) for r in rows)
    step=int(rows[0][1]);delta=rows[0][2]
    assert all(int(r[0])==i and int(r[1])==i*step and abs(r[2]-i*delta)<1e-8 for i,r in enumerate(rows,1))
    assert rows[-1][2]>=1000
    out=P/name;out.mkdir(exist_ok=True)
    shutil.copy2(d/'NuRe_2DOpenaccMultiblock.dat',out/'NuRe_2DOpenaccMultiblock.dat')
    command([str(B/'postprocess.exe')],out,'postprocess.stdout')
    ft=(out/'NuReStatistics_2DOpenaccMultiblock.dat').read_text()
    windows={}
    for label,lo,hi in [('whole:',500,1000),('first:',500,750),('last :',750,1000)]:
        want=[mean(rows,lo,hi,k,k==4) for k in range(3,8)];want[1]=math.sqrt(want[1])
        got=[float(v) for v in re.findall(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?',ft.split(label,1)[1].splitlines()[0])]
        assert len(got)==5 and all(math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-11) for a,b in zip(want,got))
        windows['{}-{}'.format(lo,hi)]=want
    whole=windows['500-1000'];first=windows['500-750'];last=windows['750-1000']
    if (d/'completion.json').exists():completion=json.loads((d/'completion.json').read_text());final=completion['final_fields'];ck=completion['final_checkpoint']
    else:completion=json.loads((d/'assessment.json').read_text());final=completion['final_fields'];ck=completion['checkpoint']
    assert ck['all_values_finite'] and final['all_nodes_finite'] and final['rho_min']>0
    r=dict(Ra=Ra,N_fine_equivalent=N,ratio=2,last_tff=rows[-1][2],source_sha256=h,
           reference_Nu=refnu,reference_Re=refre,Nu=whole[0],Re_space_time_RMS=whole[1],Nu_hot=whole[2],Nu_cold=whole[3],Nu_middle=whole[4],
           Nu_relative_error_percent=100*abs(whole[0]/refnu-1),Re_relative_error_percent=100*abs(whole[1]/refre-1),
           Nu_half_difference_percent=100*abs(last[0]-first[0])/abs(whole[0]),Re_half_difference_percent=100*abs(last[1]-first[1])/abs(whole[1]),
           max_mass_relative_drift_percent=100*max(abs(x[8]/N**2-1) for x in rows),final_mass_relative_change_percent=100*(rows[-1][8]/N**2-1),
           wall_volume_max_difference_percent=100*max(abs(whole[2]/whole[0]-1),abs(whole[3]/whole[0]-1)),windows=windows,
           source_and_runtime_hash_match=True,history_finite_contiguous=True,Fortran_postprocess_independently_verified=True,
           server_final_fields=final,server_checkpoint_verification=ck)
    result['cases'][name]=r
(P/'verified_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
