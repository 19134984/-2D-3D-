"""Use actual Fei kernels in periodic drivers: viscosity and the paper's moving temperature wave.

The wave is prescribed-velocity transport, not a high-Ra coupled RB validation.
"""
from __future__ import annotations
import argparse,datetime as dt,hashlib,json,math,re,shutil,tempfile,time
from pathlib import Path
from verify_local import compiler_environment,replace_main,replace_parameter,run_command
HERE=Path(__file__).resolve().parent

DRIVER=r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,k,n,steps
    integer,parameter :: mode=__MODE__
    real(kind=8),parameter :: advection=__ADV__,amplitude=__AMP__
    real(kind=8) :: feq(0:8),geq(0:4),waveNumber,y,quantity,as,ac
    call acc_init(acc_device_host)
    call initial()
    steps=__STEPS__
    waveNumber=2d0*pi/dble(ny)
    rho=1d0;u=0d0;v=advection;T=1d0;Fx=0d0;Fy=0d0
    do j=1,ny
      do i=1,nx
        y=dble(j)-0.5d0
        if(mode==0) u(i,j)=amplitude*sin(waveNumber*y)
        if(mode==1) T(i,j)=1d0+amplitude*sin(waveNumber*y)
        call fei_flow_equilibrium(rho(i,j),u(i,j),v(i,j),feq)
        f(i,j,:)=feq
        call fei_thermal_equilibrium(T(i,j),u(i,j),v(i,j),geq)
        g(i,j,:)=geq
      enddo
    enddo
    call enter_data_2d_openacc()
    call sample(0)
    do n=1,steps
      if(mode==0) then
        call collision()
        call streaming()
        call periodic_links(.true.)
        call macro()
      else
        call collisionT()
        call streamingT()
        call periodic_links(.false.)
        call macroT()
      endif
      if(mod(n,100)==0.or.n==steps) call sample(n)
    enddo
    call exit_data_2d_openacc()
contains
    subroutine periodic_links(flow)
      logical,intent(in) :: flow
      integer :: ii,jj,aa,ip,jp,last
      last=4
      if(flow)last=8
      !$acc parallel loop collapse(2) default(none) present(f,f_post,g,g_post,ex,ey) &
      !$acc& firstprivate(flow,last) private(aa,ip,jp) async(1)
      do jj=1,ny
        do ii=1,nx
          do aa=1,last
            ip=ii-ex(aa);jp=jj-ey(aa)
            if(ip.GE.1.AND.ip.LE.nx.AND.jp.GE.1.AND.jp.LE.ny) cycle
            ip=modulo(ip-1,nx)+1;jp=modulo(jp-1,ny)+1
            if(flow) then
              f(ii,jj,aa)=f_post(ip,jp,aa)
            else
              g(ii,jj,aa)=g_post(ip,jp,aa)
            endif
          enddo
        enddo
      enddo
    end subroutine periodic_links
    subroutine sample(step)
      integer,intent(in) :: step
      !$acc wait(1)
      !$acc update self(rho,u,v,T)
      as=0d0;ac=0d0
      do j=1,ny
        do i=1,nx
          y=dble(j)-0.5d0;quantity=u(i,j)
          if(mode==1)quantity=T(i,j)-1d0
          as=as+quantity*sin(waveNumber*y);ac=ac+quantity*cos(waveNumber*y)
        enddo
      enddo
      as=2d0*as/dble(nx*ny);ac=2d0*ac/dble(nx*ny)
      write(*,'(a,1x,i8,5(1x,es25.16e3))') 'SAMPLE',step,as,ac,sum(rho),sum(T),minval(rho)
    end subroutine sample
end program main
"""

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=HERE/'2DRBOpenaccFeiCLBM.F90')
    p.add_argument('--compiler',type=Path,required=True)
    p.add_argument('--report',type=Path,default=HERE/'transport_validation_local.json')
    p.add_argument('--timeout',type=float,default=180)
    args=p.parse_args();raw=args.source.read_bytes();source=raw.decode('utf-8')
    root=Path(tempfile.mkdtemp(prefix='fei-clbm-transport-')).resolve();env=compiler_environment(args.compiler)
    report={'passed':False,'source_sha256':hashlib.sha256(raw).hexdigest(),'verified_at_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'backend':'GNU OpenACC host; NVIDIA GPU not tested','build_directory':str(root),'cases':[],'qualification':'Periodic production kernels; no high-Ra, walls or long-time turbulence validation'}
    start=time.monotonic()
    try:
        for name,mode,adv,steps,amp in [('shear_viscosity',0,0.0,4000,1e-5),('temperature_rest',1,0.0,10132,.01),('paper_moving_temperature_Ma0p3',1,.3/math.sqrt(3),10132,.01)]:
            clone=source
            for key,value in [('nx','6'),('ny','100'),('viscosity','.05d0'),('diffusivity','.05d0'),('tauf','.65d0'),('taug','.6d0'),('Prandtl','1d0'),('gBeta','0d0'),('timeUnit','1d0'),('velocityUnit','1d0')]:clone=replace_parameter(clone,key,value)
            driver=DRIVER.replace('__MODE__',str(mode)).replace('__ADV__',f'{adv:.16f}d0').replace('__AMP__',f'{amp:.16f}d0').replace('__STEPS__',str(steps))
            directory=root/name;directory.mkdir();(directory/'solver.F90').write_text(replace_main(clone,driver))
            flags=['-cpp','-fopenacc','-O2','-ffree-line-length-none','-fcheck=all','-ffpe-trap=invalid,zero,overflow']
            run_command([str(args.compiler),*flags,'-fsyntax-only','solver.F90'],directory,env,args.timeout)
            run_command([str(args.compiler),*flags,'solver.F90','-o','solver.exe'],directory,env,args.timeout)
            result=run_command([str(directory/'solver.exe')],directory,env,args.timeout)
            (directory/'run.log').write_text(result.stdout+result.stderr)
            rows=[[float(x.replace('D','E')) for x in line.split()[1:]] for line in result.stdout.splitlines() if line.startswith('SAMPLE ')]
            assert rows and all(all(map(math.isfinite,x)) for x in rows)
            k=2*math.pi/100
            samples=[(r[0],math.log(math.hypot(r[1],r[2]))) for r in rows if r[0]>=200]
            mx=sum(x for x,y in samples)/len(samples);my=sum(y for x,y in samples)/len(samples)
            slope=sum((x-mx)*(y-my) for x,y in samples)/sum((x-mx)**2 for x,y in samples)
            measured=-slope/(k*k);error=abs(measured/.05-1)
            phase=max(abs(math.atan2(math.sin(math.atan2(r[2],r[1])+k*adv*r[0]),math.cos(math.atan2(r[2],r[1])+k*adv*r[0]))) for r in rows)
            end=rows[-1];reference=amp*math.exp(-.05*k*k*end[0]);expected=(reference*math.cos(k*adv*end[0]),-reference*math.sin(k*adv*end[0]))
            endpoint=math.hypot(end[1]-expected[0],end[2]-expected[1])/reference
            mass=max(abs(r[3]-600)/600 for r in rows);heat=max(abs(r[4]-600)/600 for r in rows)
            checks={'coefficient_within_one_percent':error<.01,'wave_endpoint_within_one_percent':endpoint<.01,'mass_conserved':mass<1e-11,'temperature_conserved':heat<1e-11,'positive_density':min(r[5] for r in rows)>0}
            case={'name':name,'passed':all(checks.values()),'checks':checks,'grid':[6,100],'steps':steps,'target_coefficient':.05,'measured_coefficient':measured,'relative_coefficient_error':error,'maximum_phase_error_radians':phase,'relative_analytical_wave_endpoint_error':endpoint,'maximum_relative_mass_drift':mass,'maximum_relative_temperature_drift':heat,'prescribed_advection':adv,'Ma':adv*math.sqrt(3),'final_t_star':.05*k*k*steps}
            report['cases'].append(case)
            assert case['passed'],name+' failed'
        report['source_unchanged']=args.source.read_bytes()==raw
        report['passed']=report['source_unchanged']
    except Exception as e:report['failure']=str(e)
    report['elapsed_seconds']=round(time.monotonic()-start,3)
    args.report.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
