"""Independent population-space checks against the actual Fortran P-TRT kernel.

The small physical cases call the production collision(), not a Python copy.
All test grids/viscosities are confined to temporary validation executables.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'PTRT/2DRBOpenaccLuoPTRT.F90'
PARENT = ROOT/'Xs/D2Q9TRT_D2Q5LuoTRT/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90'
OUT = ROOT/'运行脚本/run_cases/LuoPTRT_Ra9_Ra10_20261025'
CX=np.array([0,1,0,-1,0,1,-1,-1,1])
CY=np.array([0,0,1,0,-1,1,1,-1,-1])
W=np.array([4/9,*([1/9]*4),*([1/36]*4)])
OPP=np.array([0,3,4,1,2,7,8,5,6])
# The task's Luo ordering; the code uses rows [0,3,8,1,6,2,7,4,5].
M=np.array([[1]*9,CX,CY,[-4,-1,-1,-1,-1,2,2,2,2],
            [0,1,-1,1,-1,0,0,0,0],[0,0,0,0,0,1,-1,1,-1],
            [0,-2,0,2,0,1,-1,-1,1],[0,0,-2,0,2,1,1,-1,-1],
            [4,-2,-2,-2,-2,1,1,1,1]],dtype=float)
PERM=[0,3,8,1,6,2,7,4,5]
MC=M[PERM]
MI=np.linalg.inv(MC)


def routine(s,name):
    pat=r'^\s*subroutine '+re.escape(name)+r'\([^\n]*\).*?^\s*end subroutine '+re.escape(name)+r'\b[^\n]*'
    matches=list(re.finditer(pat,s,re.M|re.S|re.I))
    assert len(matches)==1,name
    return matches[0].group(0)


def feq(dr,rho,u):
    eu=u[0]*CX+u[1]*CY
    u2=u@u
    return W*(dr+rho*(3*eu+4.5*eu**2-1.5*u2+4.5*eu**3-4.5*eu*u2))


def independent(f,rho,u,force,snu,sq):
    eq=feq(f.sum(),rho,u)
    d=f-eq
    phi=(CX**2-1/3)*(CY**2-1/3)
    ghost=W*phi*(phi@d)/(4/81)
    dplus=.5*(d+d[OPP])-ghost
    g1=3*W*(CX*force[0]+CY*force[1])
    g2=W*(9*(CX*u[0]+CY*u[1])*(CX*force[0]+CY*force[1])-3*(u@force))
    dminus=.5*(d-d[OPP])+.5*g1
    return eq+(1-snu)*dplus+(1-sq)*dminus+.5*g1+(1-.5*snu)*g2


DRIVER='''program verify_actual_kernel
use commondata
implicit none
integer :: mode,count,l,i,j,k,step,steps,kind,strategy,srcI,srcJ,nmode,opp(0:8)
real(8) :: mm(0:8),mp(0:8),fv(0:8),gv(0:4),rr,ux,uy,fx0,fy0,sn,sq0
real(8) :: amp,wave,x,y,dr,eu,us2,exactU,exactV,err,norm,force,mean0,mean1,mmass
read(*,*) mode
opp=[0,3,4,1,2,7,8,5,6]
if(mode==1) then
 read(*,*) count
 do l=1,count
  read(*,*) mm,rr,ux,uy,fx0,fy0,sn,sq0
  call ptrt_moment_collision(mm,rr,ux,uy,fx0,fy0,sn,sq0,mp)
  write(*,'(9(ES26.17E3,1X))') mp
 enddo
 stop
endif
allocate(f(0:nx+1,0:ny+1,0:8),f_post(0:nx+1,0:ny+1,0:8))
allocate(g(0:nx+1,0:ny+1,0:4),g_post(0:nx+1,0:ny+1,0:4))
allocate(u(0:nx+1,0:ny+1),v(0:nx+1,0:ny+1),rho(0:nx+1,0:ny+1),T(0:nx+1,0:ny+1))
allocate(Fx(0:nx+1,0:ny+1),Fy(0:nx+1,0:ny+1))
allocate(Sxx(0:nx+1,0:ny+1),Sxy(0:nx+1,0:ny+1),Syy(0:nx+1,0:ny+1),Sdiv(0:nx+1,0:ny+1))
f=0; f_post=0; g=0; g_post=0; u=0; v=0; rho=rho0; T=Tref; Fx=0; Fy=0
Sxx=0; Sxy=0; Syy=0; Sdiv=0
omega=[4.0d0/9,(1.0d0/9,k=1,4),(1.0d0/36,k=5,8)]
if(mode==2) then
 read(*,*) count
 do l=1,count
  read(*,*) fv,fy0
  do k=0,8
   f(:,:,k)=fv(k)
  enddo
  T=Tref+fy0/(rho0*gBeta)
  ! Deliberately wrong cached macros: collision must reconstruct from f/current force.
  u=99; v=-99; rho=7
  call collision()
  !$acc wait(1)
  write(*,'(9(ES26.17E3,1X))') f_post(1,1,:)
 enddo
 stop
endif
if(mode==3) then
 read(*,*) count
 do l=1,count
  read(*,*) gv,ux,uy
  do k=0,4
   g(:,:,k)=gv(k)
  enddo
  T=sum(gv); u=ux; v=uy
  call collisionT()
  !$acc wait(1)
  write(*,'(5(ES26.17E3,1X))') g_post(1,1,:)
 enddo
 stop
endif
read(*,*) kind
amp=1.0d-3
wave=2*pi/dble(nx)
force=0
if(kind==1) force=1.0d-8
if(kind==2) force=8*viscosity*0.01d0/dble(nx*nx)
T=Tref+force/(rho0*gBeta)
do j=1,ny
 do i=1,nx
  x=wave*(i-0.5d0); y=wave*(j-0.5d0)
  ux=0; uy=0; dr=0
  if(kind==3) then
   ux=-amp*cos(x)*sin(y); uy=amp*sin(x)*cos(y)
   dr=-rho0*amp*amp*(cos(2*x)+cos(2*y))/(4*cs2)
  endif
  u(i,j)=ux; v(i,j)=uy; rho(i,j)=rho0+dr
  do k=0,8
   eu=ex(k)*ux+ey(k)*uy; us2=ux*ux+uy*uy
   f(i,j,k)=omega(k)*(dr+rho0*(3*eu+4.5d0*eu*eu-1.5d0*us2+ &
       4.5d0*eu**3-4.5d0*eu*us2))-1.5d0*omega(k)*ey(k)*force
  enddo
 enddo
enddo
steps=100
if(kind==2) steps=ceiling(12*dble(nx*nx)/(viscosity*pi*pi))
if(kind==3) steps=nint(0.1d0/(viscosity*wave*wave))
mean0=sum(f(1:nx,1:ny,:))
do step=1,steps
 call collision()
 !$acc wait(1)
 do k=0,8
  do j=1,ny
   do i=1,nx
    srcI=i-ex(k); srcJ=j-ey(k)
    if(kind==2.AND.(srcI.LT.1.OR.srcI.GT.nx)) then
     f(i,j,k)=f_post(i,j,opp(k))
    else
     srcI=modulo(srcI-1,nx)+1; srcJ=modulo(srcJ-1,ny)+1
     f(i,j,k)=f_post(srcI,srcJ,k)
    endif
   enddo
  enddo
 enddo
 call macro()
 !$acc wait(1)
enddo
err=0; norm=0
do j=1,ny
 do i=1,nx
  exactU=0; exactV=steps*force/rho0
  if(kind==2) then
   x=i-0.5d0
   exactV=force*x*(dble(nx)-x)/(2*viscosity*rho0)
   ! Exact startup solution: the residual transient at 12 decay times is not spatial error.
   do nmode=1,31,2
    exactV=exactV-4*force*dble(nx*nx)/(rho0*viscosity*pi**3*dble(nmode**3))* &
        sin(nmode*pi*x/dble(nx))*exp(-viscosity*(nmode*pi/dble(nx))**2*steps)
   enddo
  elseif(kind==3) then
   x=wave*(i-0.5d0); y=wave*(j-0.5d0)
   exactU=-amp*cos(x)*sin(y)*exp(-2*viscosity*wave*wave*steps)
   exactV=amp*sin(x)*cos(y)*exp(-2*viscosity*wave*wave*steps)
  endif
  err=err+(u(i,j)-exactU)**2+(v(i,j)-exactV)**2
  norm=norm+exactU**2+exactV**2
 enddo
enddo
mean1=sum(f(1:nx,1:ny,:))
write(*,'(I0,1X,4(ES26.17E3,1X))') steps,sqrt(err/norm),abs(mean1-mean0),v(1,1),v(nx/2,ny/2)
end program verify_actual_kernel
'''


def run(exe,mode,rows):
    payload=str(mode)+'\n'+str(len(rows))+'\n'+'\n'.join(' '.join(format(float(v),'.17g') for v in row) for row in rows)+'\n'
    p=subprocess.run([str(exe)],input=payload,text=True,capture_output=True,check=True,
                     env={**os.environ,'ACC_DEVICE_TYPE':'host'})
    return np.array([[float(x) for x in line.split()] for line in p.stdout.splitlines() if line.strip()])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler',default='C:/msys64/ucrt64/bin/gfortran.exe' if os.name=='nt' else 'gfortran')
    args=parser.parse_args()
    src=SOURCE.read_text(encoding='utf-8')
    parent=PARENT.read_text(encoding='utf-8-sig')
    report={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'tests':{}}
    for name in ['collisionT','streamingT','bouncebackT','macroT','streaming','bounceback']:
        assert routine(src,name)==routine(parent,name),name+' changed'
    report['tests']['unchanged_thermal_and_stream_boundary_routines']='PASS: byte-identical text'
    protected=json.loads((OUT/'protected_sources.json').read_text(encoding='utf-8'))
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in protected.items())
    report['protected_sources_unchanged']=True
    report['tests']['matrix_inverse_max_error']=float(np.max(np.abs(MC@MI-np.eye(9))))
    rng=np.random.default_rng(20261025)
    eq_err=0.0
    for _ in range(500):
        dr=rng.uniform(-.02,.02); rr=rng.uniform(.7,1.3); vel=rng.uniform(-.15,.15,2)
        ux,uy=vel
        expected=[dr,rr*ux,rr*uy,-2*dr+3*rr*(vel@vel),rr*(ux*ux-uy*uy),rr*ux*uy,
                  rr*ux*(3*uy*uy-1),rr*uy*(3*ux*ux-1),dr-3*rr*(vel@vel)]
        eq_err=max(eq_err,float(np.max(np.abs(M@feq(dr,rr,vel)-expected))))
    assert eq_err<1e-12
    report['tests']['cubic_equilibrium_max_error']=eq_err
    code_module=re.search(r'^\s*module commondata\b.*?^\s*end module commondata\b',src,re.M|re.S).group(0)
    funcs='\n'.join(routine(src,n) for n in ['ptrt_moment_collision','collision','collisionT','macro'])
    os.environ['PATH']=str(Path(args.compiler).parent)+os.pathsep+os.environ['PATH']
    with tempfile.TemporaryDirectory(prefix='luo_ptrt_',delete=False) as temp:
        build=Path(temp)
        print('Validation build:',build,flush=True)
        def compile_case(n,strategy,thermal_parent=False):
            mod=code_module.replace('nx=1024, ny=1024',f'nx={n}, ny={n}')
            mod=mod.replace('ptrtSqStrategy=1',f'ptrtSqStrategy={strategy}')
            # Explicitly confined to physical validation; high-Ra case parameters are untouched.
            mod=mod.replace('viscosity=Mach*lengthUnit*dsqrt(Prandtl/(3.0d0*Rayleigh))','viscosity=0.02d0')
            fun=funcs
            if thermal_parent:
                fun=fun.replace(routine(src,'collisionT'),routine(parent,'collisionT'))
            label=f'n{n}_p{strategy}'+('_parentT' if thermal_parent else '')
            path=build/(label+'.F90'); exe=build/(label+('.exe' if os.name=='nt' else ''))
            prefix='#define unsteadyFlow\n#define FLOW_ODD_ORIGINAL_MAGIC\n#define RayleighBenardCell\n'
            path.write_text(prefix+mod+'\n'+DRIVER+'\n'+fun,encoding='utf-8',newline='\n')
            cmd=[args.compiler,'-cpp','-fopenacc','-ffree-line-length-none','-O2','-fcheck=bounds',str(path),'-o',str(exe)]
            res=subprocess.run(cmd,cwd=build,capture_output=True,text=True)
            if res.returncode: raise RuntimeError(res.stderr)
            return exe
        exe=compile_case(16,1)
        rows=[]; expected=[]; pre=[]; forces=[]
        for index in range(1200):
            dr=rng.uniform(-.01,.01); rr=rng.uniform(.8,1.2); vel=rng.uniform(-.1,.1,2)
            force=rng.uniform(-.001,.001,2)
            if index%3==0: force[:]=0
            elif index%3==1: force[0]=0
            sn=rng.uniform(.2,1.999)
            sq0=8*(2-sn)/(8-sn) if index%2 else 1.6
            mom=MC@feq(dr,rr,vel)+rng.uniform(-.002,.002,9)
            mom[0]=dr; mom[3]=rr*vel[0]-.5*force[0]; mom[5]=rr*vel[1]-.5*force[1]
            fp=MI@mom
            rows.append([*mom,rr,*vel,*force,sn,sq0])
            expected.append(independent(fp,rr,vel,force,sn,sq0)); pre.append(mom); forces.append(force)
        got=run(exe,1,rows); ref=np.array(expected)@MC.T
        diff=float(np.max(np.abs(got-ref)))
        mass=float(np.max(np.abs(got[:,0]-np.array(pre)[:,0])))
        impulse=float(np.max(np.abs(got[:,[3,5]]-np.array(pre)[:,[3,5]]-forces)))
        ghost=float(np.max(np.abs(got[:,0]+got[:,1]+got[:,2])))
        assert max(diff,mass,impulse,ghost)<1e-12
        report['tests']['actual_fortran_kernel_1200_states']=dict(projection_error=diff,mass_error=mass,
            impulse_error=impulse,ghost_error=ghost,includes='zero/vertical/general force, P1 small Sq, P2')
        rows=[]; expected=[]
        sn=1/(.5+3*.02); sq0=8*(2-sn)/(8-sn)
        for index in range(80):
            vel=rng.uniform(-.04,.04,2); force=np.array([0.0,rng.uniform(-1e-5,1e-5)])
            ff=feq(rng.uniform(-.003,.003),1,vel)+rng.uniform(-.0001,.0001,9)
            actual_vel=np.array([CX@ff,CY@ff])+.5*force
            rows.append([*ff,force[1]])
            expected.append(independent(ff,1,actual_vel,force,sn,sq0))
        diff=float(np.max(np.abs(run(exe,2,rows)-expected)))
        assert diff<1e-12
        report['tests']['production_collision_and_force_no_double_count_max_error']=diff
        thermalrows=[[*rng.uniform(-.1,.1,5),*rng.uniform(-.1,.1,2)] for _ in range(80)]
        thermalA=run(exe,3,thermalrows)
        thermalB=run(compile_case(16,1,True),3,thermalrows)
        diff=float(np.max(np.abs(thermalA-thermalB)))
        assert diff==0
        report['tests']['actual_D2Q5_same_velocity_max_difference']=diff
        phys={}
        for strategy in [1,2]:
            for n in [16,32]:
                phys[f'P{strategy}_N{n}']={}
                candidate=exe if (strategy,n)==(1,16) else compile_case(n,strategy)
                for kind,name in [(1,'uniform_force_periodic'),(2,'poiseuille_halfway'),(3,'taylor_green')]:
                    res=subprocess.run([str(candidate)],input=f'4\n{kind}\n',text=True,capture_output=True,check=True,
                        env={**os.environ,'ACC_DEVICE_TYPE':'host'},timeout=120)
                    vals=[float(x) for x in res.stdout.split()]
                    steps,err,mass,wall,center=vals
                    print(f'P{strategy} N{n} {name}: steps={steps:g} relative_L2={err:.8g} mass={mass:.4g}',flush=True)
                    assert np.isfinite(vals).all()
                    if kind==1: assert err<1e-9 and mass<1e-10
                    else: assert err<(.04 if n==16 else .01),(strategy,n,name,vals)
                    phys[f'P{strategy}_N{n}'][name]=dict(steps=int(steps),relative_L2_error=err,
                        total_mass_change=mass,first_cell_v=wall,center_v=center)
            for name in ['poiseuille_halfway','taylor_green']:
                a=phys[f'P{strategy}_N16'][name]['relative_L2_error']
                b=phys[f'P{strategy}_N32'][name]['relative_L2_error']
                assert b<a,(strategy,name,a,b)
                if name=='taylor_green' or strategy==2:
                    assert math.log2(a/b)>1.7,(strategy,name,a,b)
                phys[f'P{strategy}_convergence_{name}']=dict(error_reduction=a/b,observed_order=math.log2(a/b))
        report['tests']['physical_cases_actual_Fortran_collision']=phys
    report['pass']=True
    report['note']='Local OpenACC host validation; P100 GPU smoke and high-Ra completion are separate checks.'
    (OUT/'local_validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
