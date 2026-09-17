"""Read-only audit of the completed first multiblock case and its final checkpoint."""
from pathlib import Path
import subprocess,sys,json,re,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
CASE=ROOT/'运行脚本/Multiblock_sideheated_20260916_v3/Ra1e6/ratio2/sideheated_N256_20260916_v3'
R=CASE/'results'
name=(R/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().strip()
assert re.fullmatch(r'reloadFile2DOpenaccMultiblock-\d+\.bin',name)
if '--fetch' in sys.argv:
    ssh=['C:/msys64/usr/bin/sshpass.exe','-e','C:/msys64/usr/bin/ssh.exe','-o','ConnectTimeout=15','xlli@10.169.174.118']
    with (R/name).open('wb') as out:
        subprocess.run(ssh+[f'cat /data2/XLLi/Multiblock/Ra1e6/ratio2/sideheated_N256_20260916_v3/results/{name}'],stdout=out,check=True)
blob=(R/name).read_bytes()
pos=16
assert blob[:16].decode().strip()=='MB2DRESTART0010'
def take(dtype,n):
    global pos
    a=np.frombuffer(blob,dtype=dtype,count=n,offset=pos)
    pos+=a.nbytes
    return a
grid=take('<i4',9)
nx,ny,ratio,left,right,bottom,top,overlap,nb=map(int,grid)
model=take('<i4',12)
physics=take('<f8',16)
counters=take('<i4',7)
errors=take('<f8',2)
assert (nx,ny,ratio,left,right,bottom,top)==(256,256,2,32,33,32,33)
assert np.allclose(physics[:3],[1e6,.71,.1],rtol=0,atol=0)
assert int(counters[0])==834000
blocks=[]
for b in range(nb):
    ni,nj,ilo,ihi,jlo,jhi,nh=map(int,take('<i4',7))
    x0,y0,h,*box=take('<f8',7)
    fields={}
    for key,nc in [('f',9),('g',5),('u',1),('v',1),('T',1),('rho',1),('Fx',1),('Fy',1),('p',20*(nh+1)),('up',1),('vp',1),('Tp',1)]:
        arr=take('<f8',ni*nj*nc)
        assert np.isfinite(arr).all(),(b,key)
        fields[key]=arr.reshape((ni,nj,nc),order='F') if nc>1 else arr.reshape((ni,nj),order='F')
    x=x0+(np.arange(ni)+.5)*h
    y=y0+(np.arange(nj)+.5)*h
    def area(bounds):
        xa,xb,ya,yb=bounds
        wx=np.maximum(0,np.minimum(x+h/2,xb)-np.maximum(x-h/2,xa))
        wy=np.maximum(0,np.minimum(y+h/2,yb)-np.maximum(y-h/2,ya))
        return wx[:,None]*wy[None,:]
    weights=area(box)
    if b==1:
        weights=weights-area(blocks[0]['box'])
    blocks.append(dict(box=box,fields=fields,weights=weights))
assert pos==len(blob),(pos,len(blob))
area_total=sum(b['weights'].sum() for b in blocks)
assert abs(area_total-nx*ny)<1e-8
nu=.1*nx*np.sqrt(3*.71/1e6)/3
alpha=nu/.71
time_unit=nx/(.1/np.sqrt(3))
conv=sum(np.sum(b['weights']*b['fields']['u']*b['fields']['T']) for b in blocks)/area_total
vel2=sum(np.sum(b['weights']*(b['fields']['u']**2+b['fields']['v']**2)) for b in blocks)/area_total
mass=sum(np.sum(b['weights']*b['fields']['rho']) for b in blocks)
final={'t_ff':int(counters[0])/time_unit,'NuVolAvg':1+conv*nx/alpha,'ReVolRMS':np.sqrt(vel2)*nx/nu}
temp=blocks[1]['fields']['T']
final['Nu_hot']=float(np.sum((4.0-9*temp[0,:]+temp[1,:])/3))
final['Nu_cold']=float(np.sum((4.0+9*temp[-1,:]-temp[-2,:])/3))
history=np.loadtxt(R/'NuRe_2DOpenaccMultiblock.dat')
checks=np.loadtxt(R/'Convergence_2DOpenaccMultiblock.dat')
assert np.isfinite(history).all() and np.isfinite(checks).all()
assert np.all(np.diff(checks[:,0])==2000)
assert np.all(np.diff(history[:,0])>0)
expected=np.ceil(np.arange(1,len(history)+1)*10*time_unit/ratio)*ratio/time_unit
assert np.allclose(history[:,0],expected,rtol=0,atol=1e-10)
assert int(counters[1])-1==len(history)
assert np.array_equal(checks[-1,1:],errors)
assert np.all(errors<=1e-7)
for file in ['solver.stdout','pbs.stderr']:
    s=(R/file).read_text(errors='replace')
    assert not re.search(r'(?i)(?<![a-z])(nan|infinity|inf)(?![a-z])|CUDA_ERROR|Error:|ERROR STOP|FATAL',s)
assert (R/'compile.status').read_text().strip()=='0'
assert (R/'run.status').read_text().strip()=='0'
casehash=hashlib.sha256((CASE/'source/solver.F90').read_bytes()).hexdigest()
assert (R/'runtime_source.sha256').read_text().split()[0]==casehash
refs={'paper_257':{'NuVolAvg':8.8282,'ReVolRMS':99.6110,'Nu_hot':8.8509},
      'paper_extrapolated':{'NuVolAvg':8.8186,'ReVolRMS':99.3542,'Nu_hot':8.8206}}
deviations={name:{k:100*(final[k]-val)/val for k,val in row.items()} for name,row in refs.items()}
report={'job':'6979.master','state':'completed_steady_threshold_passed','fine_steps':int(counters[0]),
        'errors':errors.tolist(),'final_checkpoint_metrics':final,'last_history_tff':history[-1,0],
        'history_rows':len(history),'convergence_rows':len(checks),'checkpoint_all_fields_finite':True,
        'history_continuous_finite':True,'mass_relative_drift_percent':100*(mass/area_total-1),
        'hot_cold_difference_percent':100*abs(final['Nu_hot']-final['Nu_cold'])/final['Nu_hot'],
        'paper_signed_errors_percent':deviations,'source_sha256':casehash,
        'checkpoint_sha256':hashlib.sha256(blob).hexdigest()}
(CASE/'验收结果.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
