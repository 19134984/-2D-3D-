"""Compare a coupled side-heated transient with the unchanged uniform-grid parent.
Imports the verification helpers; all generated executables and states remain in TEMP.
"""
import json
import numpy as np
import verify_multiblock as v

src=v.SOURCE.read_text(encoding='utf-8-sig')
parent=v.PARENT.read_text(encoding='utf-8-sig')
n=96
steps=2000
pd,pe=v.compile_source('parent_transient',v.variant(parent,side=True),
    v.uniform_driver(True).replace('step=1,40',f'step=1,{steps}'),ra=10000)
md,mexe=v.compile_source('mb_transient',v.variant(src,side=True),v.uniform_driver().replace('step=1,40',f'step=1,{steps//2}')
    .replace(v.state_writer(), '    call calNuRe()\n    call output_SnapshotFile()\n'),ra=10000)
v.run([pe],pd,timeout=180)
v.run([mexe],md,timeout=180)
raw=np.fromfile(pd/'state.bin',dtype='<f8')
u,vel,tt,rho=[raw[k*n*n:(k+1)*n*n].reshape((n,n),order='F') for k in range(14,18)]
ref=[u,vel,tt,rho]
num=np.zeros(4); denom=np.zeros(4); maxerr=np.zeros(4)
with next(md.glob('*Snapshot-*.bin')).open('rb') as f:
    magic=f.read(16)
    assert magic == b'MB2DSNAPSHOT0002'
    nb,nx,ny,itc=np.fromfile(f,dtype='<i4',count=4)
    tff,L=np.fromfile(f,dtype='<f8',count=2)
    for b in range(nb):
        ni,nj=np.fromfile(f,dtype='<i4',count=2)
        first_x,first_y,h,*owned_box=np.fromfile(f,dtype='<f8',count=7)
        dx=np.fromfile(f,dtype='<f8',count=ni); dy=np.fromfile(f,dtype='<f8',count=nj)
        area_weights=dx[:,None]*dy[None,:]
        arr=np.fromfile(f,dtype='<f8',count=4*ni*nj).reshape((ni,nj,4),order='F')
        xx=first_x+np.arange(ni)*h-.5
        yy=first_y+np.arange(nj)*h-.5
        ix=np.floor(xx).astype(int); iy=np.floor(yy).astype(int)
        wx=xx-ix; wy=yy-iy
        ix2=np.minimum(ix+1,n-1); iy2=np.minimum(iy+1,n-1)
        for a,field in enumerate(ref):
            interp=((1-wx)[:,None]*(1-wy)[None,:]*field[np.ix_(ix,iy)]+
                    wx[:,None]*(1-wy)[None,:]*field[np.ix_(ix2,iy)]+
                    (1-wx)[:,None]*wy[None,:]*field[np.ix_(ix,iy2)]+
                    wx[:,None]*wy[None,:]*field[np.ix_(ix2,iy2)])
            num[a]+=np.sum((arr[:,:,a]-interp)**2*area_weights)
            denom[a]+=np.sum(interp**2*area_weights)
            maxerr[a]=max(maxerr[a],np.max(np.abs(arr[:,:,a]-interp)))
nuFine=.1*n*np.sqrt(3*.7/10000)/3
diff=nuFine/.7
uniformNu=1+np.mean(u*tt)*n/diff
uniformRe=np.sqrt(np.mean(u*u+vel*vel))*n/nuFine
uniformHot=np.mean((8*.5-9*tt[0,:]+tt[1,:])/3)*n
uniformCold=np.mean((-8*(-.5)+9*tt[-1,:]-tt[-2,:])/3)*n
coarse=np.loadtxt(md/'NuRe_2DOpenaccMultiblock.dat')
metrics=np.array([uniformNu,uniformRe,uniformHot,uniformCold])
actual=coarse[1:5]
result={'case':'side-heated, EnableUseG, Ra=1e4, Pr=0.7, Ma=0.1',
        'layout':'aligned nested nodes, snapshot v2, clipped integration weights',
        'fine_equivalent_grid':[n,n],'fine_steps':steps,'t_ff':float(tff),
        'field_order':['u','v','T','rho'],'field_relative_l2':np.sqrt(num/denom).tolist(),
        'field_max_abs':maxerr.tolist(),'metric_order':['NuVolAvg','ReVolRMS','Nu_hot','Nu_cold'],
        'uniform':metrics.tolist(),'multiblock':actual.tolist(),
        'relative_metric_difference':(np.abs(actual-metrics)/np.abs(metrics)).tolist(),
        'source_sha256':v.hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
        'parent_sha256':v.hashlib.sha256(v.PARENT.read_bytes()).hexdigest(),
        'device':'OpenACC host','build_directory':str(v.BUILD),
        'scope':'Transient comparison only; this is not a converged benchmark or GPU validation.'}
assert np.all(np.isfinite(actual)) and np.all(np.isfinite(maxerr))
(v.HERE/'transient_comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
