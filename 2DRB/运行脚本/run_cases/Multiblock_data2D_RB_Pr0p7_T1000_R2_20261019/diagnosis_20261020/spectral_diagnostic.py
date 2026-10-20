"""Frozen-velocity D2Q9 MRT Fourier diagnostic, not a full coupled-flow stability proof."""
import numpy as np,json
from pathlib import Path
M=np.array([[1,1,1,1,1,1,1,1,1],[-4,-1,-1,-1,-1,2,2,2,2],[4,-2,-2,-2,-2,1,1,1,1],[0,1,0,-1,0,1,-1,-1,1],[0,-2,0,2,0,1,-1,-1,1],[0,0,1,0,-1,1,1,-1,-1],[0,0,-2,0,2,1,1,-1,-1],[0,1,-1,1,-1,0,0,0,0],[0,0,0,0,0,1,-1,1,-1]],float)
Mi=np.linalg.inv(M);ex=np.array([0,1,0,-1,0,1,-1,-1,1]);ey=np.array([0,0,1,0,-1,1,1,-1,-1])
k=np.linspace(-np.pi,np.pi,65);kx,ky=np.meshgrid(k,k);phase=np.exp(-1j*(kx[:,:,None]*ex+ky[:,:,None]*ey)).reshape(-1,9)
out=[]
for exponent,n in [(6,128),(7,256),(8,512),(9,1024),(10,2048)]:
 sn=1/(.5+.1*n*np.sqrt(3*.7/10**exponent));sq=8*(2-sn)/(8-sn);snc=1/(.5+(1/sn-.5)/2)
 for mode,sqc in [('current',1/(.5+(1/sq-.5)/2)),('magic_relation',8*(2-snc)/(8-snc))]:
  for u,v in [(.02,0),(.03,0),(.02,.02)]:
   S=np.diag([0,snc,snc,0,sqc,0,sqc,snc,snc]);E=np.zeros((9,9))
   E[:,0]=[1,-2-3*(u*u+v*v),1+3*(u*u+v*v),0,0,0,0,-u*u+v*v,-u*v]
   E[:,3]=[0,6*u,-6*u,1,-1,0,0,2*u,v];E[:,5]=[0,6*v,-6*v,0,0,1,-1,-2*v,u]
   C=Mi@(np.eye(9)-S+S@E)@M;eig=np.linalg.eigvals(phase[:,:,None]*C[None,:,:]);which=np.unravel_index(np.argmax(abs(eig)),eig.shape);radius=float(abs(eig[which]));whichk=np.unravel_index(which[0],kx.shape)
   row=dict(Ra=10**exponent,n=n,mode=mode,SnuCoarse=snc,SqCoarse=sqc,ux=u,uy=v,spectral_radius=radius,kx=float(kx[whichk]),ky=float(ky[whichk]));out.append(row)
   if exponent>=8 and (u,v)==(.03,0):print(json.dumps(row))
Path(__file__).with_name('spectral_comparison.json').write_text(json.dumps(out,indent=2))
