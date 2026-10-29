"""Actual Fei D2Q9/D2Q5 central-moment kernels: equations, walls, output and restart.
All temporary numerical settings are applied to isolated source copies.
GNU OpenACC host verification does not certify NVIDIA GPU execution.
"""
from __future__ import annotations
import argparse,datetime as dt,hashlib,json,math,os,re,shutil,struct,subprocess,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
NX,NY=12,10
def replace_parameter(source: str, name: str, value: str) -> str:
    pattern = rf"(?im)(\b{re.escape(name)}\s*=\s*)[^,!\n]+"
    result, count = re.subn(pattern, lambda m: m[1] + value, source, count=1)
    if count != 1:
        raise AssertionError(f"Parameter {name} not found")
    return result


def select_macros(source: str, *, steady: bool, rb: bool) -> str:
    for macro, enabled in {"steadyFlow": steady, "unsteadyFlow": not steady,
                           "SideHeatedCell": not rb, "RayleighBenardCell": rb}.items():
        pattern = rf"(?m)^[ \t]*!?[ \t]*#define[ \t]+{macro}[ \t]*(?:![^\n]*)?$"
        source, count = re.subn(pattern, ("#define " if enabled else "!#define ") + macro, source)
        if count != 1:
            raise AssertionError(f"Expected one case macro {macro}, found {count}")
    return source


def replace_main(source: str, main: str) -> str:
    pattern = r"(?ims)^\s*program\s+main\s*$.*?^\s*end\s+program\s+main\s*$"
    result, count = re.subn(pattern, lambda _: "\n" + main + "\n", source, count=1)
    if count != 1:
        raise AssertionError("Expected exactly one production program main")
    return result


def small_source(source: str, *, steady: bool, rb: bool, steps: int,
                 restart: bool = False) -> str:
    source = select_macros(source, steady=steady, rb=rb)
    for name, value in (("nx", str(NX)), ("ny", str(NY)), ("Rayleigh", "1.0d4"),
                        ("loadInitField", "1" if restart else "0"),
                        ("steadyCheckIntervalItc", "11")):
        source = replace_parameter(source, name, value)
    for name, value in (("outputSnapshotInterval", "7.0d0/timeUnit"),
                        ("reloadFileInterval", "16.0d0/timeUnit"),
                        ("outputPltFileInterval", "16.0d0/timeUnit")):
        pattern = rf"(?im)^(\s*real\(kind=8\),\s*parameter\s*::\s*{name}\s*=\s*)[^!\n]+"
        source, count = re.subn(pattern, lambda m: m[1] + value + "   ", source)
        if not count:
            raise AssertionError(f"Output parameter {name} not found")
    source = replace_parameter(source, "unsteadyRunDuration", f"{steps}.0d0/timeUnit")
    if steady:
        source = replace_parameter(source, "itc_max", str(steps))
    return source


def compiler_environment(compiler: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = str(compiler.parent) + os.pathsep + env.get("PATH", "")
    env["ACC_DEVICE_TYPE"] = "host"
    env["OMP_NUM_THREADS"] = "1"
    env["LC_ALL"] = "C"
    # macOS 系统 Python 可能清除 DYLD_*；可用此普通变量显式传入便携运行库路径。
    if env.get("FENG_FORTRAN_RUNTIME_LIBRARY_PATH"):
        env["DYLD_LIBRARY_PATH"] = env["FENG_FORTRAN_RUNTIME_LIBRARY_PATH"]
    return env


def run_command(command: list[str], directory: Path, env: dict[str, str],
                timeout: float, *, expect_success: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, cwd=directory, env=env, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=timeout)
    if expect_success and result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(command)}\n"
                           + (result.stdout + result.stderr)[-12000:])
    return result


def build_run(source: str, directory: Path, compiler: Path, env: dict[str, str],
              timeout: float, *, syntax_only: bool = False) -> str:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "solver.F90").write_text(source, encoding="utf-8")
    flags = ["-cpp", "-ffree-line-length-none", "-fopenacc", "-O0", "-g",
             "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
    syntax = run_command([str(compiler), *flags, "-fsyntax-only", "solver.F90"], directory, env, timeout)
    (directory / "syntax.log").write_text(syntax.stdout + syntax.stderr, encoding="utf-8")
    if syntax_only:
        return ""
    compiled = run_command([str(compiler), *flags, "solver.F90", "-o", "solver.exe"], directory, env, timeout)
    (directory / "compile.log").write_text(compiled.stdout + compiled.stderr, encoding="utf-8")
    result = run_command([str(directory / "solver.exe")], directory, env, timeout)
    output = result.stdout + result.stderr
    (directory / "run.log").write_text(output, encoding="utf-8")
    return output


def read_records(path: Path) -> list[bytes]:
    data = path.read_bytes()
    records: list[bytes] = []
    offset = 0
    while offset < len(data):
        size, = struct.unpack_from("<i", data, offset)
        if size < 0 or offset + size + 8 > len(data):
            raise AssertionError(f"Invalid sequential-unformatted record: {path}")
        end, = struct.unpack_from("<i", data, offset + size + 4)
        if size != end:
            raise AssertionError(f"Record marker mismatch: {path}")
        records.append(data[offset + 4:offset + size + 4])
        offset += size + 8
    return records


def metadata(directory: Path) -> tuple[Path, dict[str, str]]:
    files = list(directory.glob("reload*-latest.meta"))
    if len(files) != 1:
        raise AssertionError("Expected one latest checkpoint metadata")
    entries = {}
    for line in files[0].read_text().splitlines():
        pair = line.split(maxsplit=1)
        if len(pair) == 2:
            entries[pair[0]] = pair[1].strip()
    path = files[0].with_name(files[0].name.replace("latest.meta", entries["reloadFileName"] + ".bin"))
    return path, entries


MATH_MAIN = r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,k,n
    integer,parameter :: cx(0:8)=[0,1,0,-1,0,1,-1,-1,1],cy(0:8)=[0,0,1,0,-1,1,1,-1,-1]
    real(kind=8),parameter :: wi(0:8)=[4d0/9,1d0/9,1d0/9,1d0/9,1d0/9,1d0/36,1d0/36,1d0/36,1d0/36]
    real(kind=8) :: pf(0:8),qf(0:8),pt(0:4),qt(0:4),cm(0:8),cmPost(0:8),eq(0:8),src(0:8)
    real(kind=8) :: matrix(0:8,0:8),expectedCM(0:8),ct(0:4),ctPost(0:4),ctExpected(0:4)
    real(kind=8) :: ux,uy,forceX,forceY,density,scalar,qx,qy,err,expected,initialHeat
    call acc_init(acc_device_host)
    call initial()
    call require(thermalCs2==0.5d0,'paper_cT_squared_is_one_half')
    call require(size(g,3)==5.and.size(f,3)==9,'separate_D2Q5_D2Q9_shapes')
    do n=1,80
      pf=wi*(1d0+0.21d0*[(sin(0.77d0*dble(n)+1.3d0*dble(k)),k=0,8)])
      density=sum(pf)
      forceX=0.007d0*cos(dble(n));forceY=0.009d0*sin(dble(n))
      ux=(sum(pf*dble(cx))+0.5d0*forceX)/density
      uy=(sum(pf*dble(cy))+0.5d0*forceY)/density
      call direct_central(pf,ux,uy,cm)
      matrix=0d0
      matrix(0,0)=1d0;matrix(1,1)=1d0;matrix(2,2)=1d0
      matrix(3,3)=0.5d0*(feiBulkRate+Snu);matrix(4,4)=matrix(3,3)
      matrix(3,4)=0.5d0*(feiBulkRate-Snu);matrix(4,3)=matrix(3,4)
      matrix(5,5)=Snu;matrix(6,6)=feiThirdRate;matrix(7,7)=feiThirdRate;matrix(8,8)=feiFourthRate
      eq=0d0;eq(0)=density;eq(3)=density/3d0;eq(4)=density/3d0;eq(8)=density/9d0
      src=0d0;src(1)=forceX;src(2)=forceY;src(6)=forceY/3d0;src(7)=forceX/3d0
      expectedCM=cm-matmul(matrix,cm-eq)+src-0.5d0*matmul(matrix,src)
      call fei_flow_collision(pf,ux,uy,forceX,forceY,qf)
      call direct_central(qf,ux,uy,cmPost)
      call measure('flow_all_central_moments_Eq8',maxval(abs(cmPost-expectedCM)),8d-14)
      call measure('flow_mass',abs(sum(qf)-sum(pf)),5d-14)
      call measure('flow_x_force_increment',abs(sum(qf*dble(cx))-sum(pf*dble(cx))-forceX),5d-14)
      call measure('flow_y_force_increment',abs(sum(qf*dble(cy))-sum(pf*dble(cy))-forceY),5d-14)
      call fei_flow_equilibrium(density,ux,uy,qf)
      call direct_central(qf,ux,uy,cmPost)
      call measure('continuous_flow_central_equilibrium',maxval(abs(cmPost-eq)),8d-14)
      scalar=(-1d0)**n*(0.2d0+0.01d0*dble(n));ux=0.12d0*cos(dble(n));uy=0.09d0*sin(dble(n))
      call fei_thermal_equilibrium(scalar,ux,uy,pt)
      call direct_thermal(pt,ux,uy,ct)
      ctExpected=[scalar,0d0,0d0,scalar*0.5d0,scalar*0.5d0]
      call measure('thermal_continuous_central_equilibrium',maxval(abs(ct-ctExpected)),5d-14)
      pt=pt+0.013d0*[(sin(0.53d0*dble(n)+dble(k)),k=0,4)]
      call direct_thermal(pt,ux,uy,ct)
      ctExpected=ct
      ctExpected(1:2)=(1d0-feiThermalRate)*ct(1:2)
      ctExpected(3:4)=(1d0-feiThermalHighRate)*ct(3:4)+feiThermalHighRate*ct(0)*0.5d0
      call fei_thermal_collision(pt,ux,uy,qt)
      call direct_thermal(qt,ux,uy,ctPost)
      call measure('thermal_all_central_moments_Eq19',maxval(abs(ctPost-ctExpected)),5d-14)
      call measure('thermal_scalar_conservation',abs(sum(qt)-sum(pt)),5d-14)
    enddo
    ! Actual wrapper at nonunit density: Fei Eq3 uses reference rho0=1 force density.
    density=1.7d0;forceY=gBeta*0.2d0
    rho=density;T=Tref+0.2d0;u=0d0;v=0.5d0*forceY/density
    call fei_flow_equilibrium(density,0d0,0d0,pf)
    do k=0,8
      f(:,:,k)=pf(k)
    enddo
    call enter_data_2d_openacc()
    call collision()
    !$acc wait(1)
    !$acc update self(f_post,Fy)
    call measure('reference_density_force_wrapper',maxval(abs(Fy-forceY)),5d-14)
    call measure('reference_density_momentum_wrapper', &
         abs(sum(f_post(1,1,:)*dble(cy))-forceY),5d-14)
    f=f_post(1:nx,1:ny,:)
    !$acc update device(f)
    call macro()
    !$acc wait(1)
    !$acc update self(v)
    call measure('new_force_half_step_wrapper',maxval(abs(v-1.5d0*forceY/density)),5d-14)
    call exit_data_2d_openacc()
    ! Actual device-resident wrappers and D2Q5 walls: signed linear conduction.
    rho=1d0;u=0d0;v=0d0
    do j=1,ny
      do i=1,nx
#ifdef SideHeatedCell
        T(i,j)=Thot+(dble(i)-0.5d0)/dble(nx)*(Tcold-Thot)
#else
        T(i,j)=Thot+(dble(j)-0.5d0)/dble(ny)*(Tcold-Thot)
#endif
        call fei_thermal_equilibrium(T(i,j),0d0,0d0,pt)
        g(i,j,:)=pt
      enddo
    enddo
    initialHeat=sum(T)
    call enter_data_2d_openacc()
    do n=1,600
      call collisionT()
      call streamingT()
      call bouncebackT()
      call macroT()
    enddo
    !$acc wait(1)
    !$acc update self(T,g)
    err=0d0
    do j=1,ny
      do i=1,nx
#ifdef SideHeatedCell
        expected=Thot+(dble(i)-0.5d0)/dble(nx)*(Tcold-Thot)
#else
        expected=Thot+(dble(j)-0.5d0)/dble(ny)*(Tcold-Thot)
#endif
        err=max(err,abs(T(i,j)-expected))
      enddo
    enddo
    call measure('linear_conduction_temperature',err,2d-10)
    call measure('linear_conduction_heat_balance',abs(sum(T)-initialHeat),2d-9)
#ifdef SideHeatedCell
    call SideHeatedcalc_Nu_global()
    call SideHeatedcalc_Nu_wall_avg()
#else
    call RBcalc_Nu_global()
    call RBcalc_Nu_wall_avg()
#endif
    call measure('conduction_Nu_volume',abs(Nu_global-1d0),2d-9)
    call measure('conduction_Nu_hot',abs(Nu_hot-1d0),2d-9)
    call measure('conduction_Nu_cold',abs(Nu_cold-1d0),2d-9)
    call exit_data_2d_openacc()
    write(*,'(a)') 'ALL_MATHEMATICAL_CHECKS_PASSED'
contains
    subroutine direct_central(pop,x,y,m)
      real(kind=8),intent(in) :: pop(0:8),x,y
      real(kind=8),intent(out) :: m(0:8)
      real(kind=8) :: a,b
      integer :: ii
      m=0d0
      do ii=0,8
        a=dble(cx(ii))-x;b=dble(cy(ii))-y
        m=m+pop(ii)*[1d0,a,b,a*a,b*b,a*b,a*a*b,a*b*b,a*a*b*b]
      enddo
    end subroutine direct_central
    subroutine direct_thermal(pop,x,y,m)
      real(kind=8),intent(in) :: pop(0:4),x,y
      real(kind=8),intent(out) :: m(0:4)
      integer :: ii
      m=0d0
      do ii=0,4
        m=m+pop(ii)*[1d0,dble(cx(ii))-x,dble(cy(ii))-y,(dble(cx(ii))-x)**2,(dble(cy(ii))-y)**2]
      enddo
    end subroutine direct_thermal
    subroutine measure(label,residual,limit)
      character(len=*),intent(in) :: label
      real(kind=8),intent(in) :: residual,limit
      write(*,'(a,1x,a,2(1x,es24.16e3))') 'CHECK',label,residual,limit
      if(.not.(residual<=limit)) error stop 'Mathematical residual exceeds tolerance'
    end subroutine measure
    subroutine require(condition,label)
      logical,intent(in) :: condition
      character(len=*),intent(in) :: label
      if(.not.condition) then
        write(*,'(a)') 'FAILED '//label
        error stop 42
      endif
    end subroutine require
end program main
"""

def inspect_run(directory):
    state,meta=metadata(directory)
    records=read_records(state)
    assert len(records[0])==48 and records[0][:32].rstrip()==b'FEI_CLBM_2D_V1'
    assert struct.unpack('<4i',records[0][32:])==(NX,NY,9,5)
    assert [len(x) for x in records[1:3]]==[NX*NY*9*8,NX*NY*5*8]
    values=[struct.unpack(f'<{len(x)//8}d',x) for x in records[1:3]]
    assert all(math.isfinite(v) for record in values for v in record)
    mass_error=abs(sum(values[0])-NX*NY)/(NX*NY)
    assert mass_error<1e-11
    histories={}
    for stem in ['Nu_VolAvg','Re_VolAvg']:
        paths=list(directory.glob(stem+'*.dat'));assert len(paths)==1
        rows=[tuple(map(float,line.split())) for line in paths[0].read_text().splitlines() if line.strip()]
        assert len(rows)==int(meta['nureSampleCount']) and all(len(x)==2 and all(map(math.isfinite,x)) for x in rows)
        assert all(a[0]<b[0] for a,b in zip(rows,rows[1:]))
        histories[paths[0].name]=paths[0].read_bytes()
    snapshots={p.name:p.read_bytes() for p in directory.glob('*Snapshot-*.bin')}
    assert snapshots
    for name,data in snapshots.items():
        assert [len(x) for x in read_records(directory/name)]==[NX*NY*8]*4
    return state,meta,histories,snapshots,{'itc_total':int(meta['itc_total']),'relative_mass_error':mass_error,'checkpoint_record_bytes':[len(x) for x in records]}

def integration_checks(source,root,compiler,env,timeout):
    report={}
    for name,steady,rb in [('side_steady',True,False),('rb_unsteady',False,True)]:
        full=root/name;split=root/(name+'_restart')
        build_run(small_source(source,steady=steady,rb=rb,steps=48),full,compiler,env,timeout)
        build_run(small_source(source,steady=steady,rb=rb,steps=32),split,compiler,env,timeout)
        _,_,old_history,_,_=inspect_run(split)
        seed=root/(name+'_seed32');shutil.copytree(split,seed)
        build_run(small_source(source,steady=steady,rb=rb,steps=48,restart=True),split,compiler,env,timeout)
        fullstate,fmeta,fh,fs,fmetrics=inspect_run(full)
        splitstate,smeta,sh,ss,smetrics=inspect_run(split)
        assert fullstate.read_bytes()==splitstate.read_bytes(),name+' checkpoint mismatch'
        assert all(sh[n].startswith(data) for n,data in old_history.items())
        assert all(sh[n].splitlines()[-1]==data.splitlines()[-1] for n,data in fh.items())
        assert all(ss.get(n)==data for n,data in fs.items())
        if not steady:assert ss==fs and sh==fh
        protected={p.name:p.read_bytes() for p in split.iterdir() if p.suffix=='.bin' or p.name.endswith('.meta') or p.name in sh}
        run_command([str(split/'solver.exe')],split,env,timeout)
        assert all((split/n).read_bytes()==v for n,v in protected.items())
        rejected={}
        for corruption in ['missing_Nu','extra_Nu','wrong_time_Nu','wrong_model_meta']:
            d=root/(name+'_'+corruption);shutil.copytree(seed,d)
            path=d/'Nu_VolAvg_2DOpenaccFeiCLBM.dat'
            if corruption=='missing_Nu':path.unlink()
            elif corruption=='extra_Nu':path.write_text(path.read_text()+path.read_text().splitlines()[-1]+'\n')
            elif corruption=='wrong_time_Nu':
                rows=path.read_text().splitlines();row=rows[0].split();row[0]='9.9E9';rows[0]=' '.join(row);path.write_text('\n'.join(rows)+'\n')
            else:
                meta=d/next(d.glob('*latest.meta')).name;meta.write_text(meta.read_text().replace('FEI_CLBM_2D_V1','WRONG_SCHEME'))
            r=run_command([str(split/'solver.exe')],d,env,timeout,expect_success=False)
            assert r.returncode!=0 and ('FeiCLBM' in r.stdout+r.stderr),corruption
            rejected[corruption]=True
        report[name]={'continuous':fmetrics,'restart':smetrics,'checkpoint_bitwise_equal':True,'common_snapshots_bitwise_equal':True,'history_preserved':True,'completed_target_preserved':True,'rejected_corrupt_restarts':rejected}
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=HERE/'2DRBOpenaccFeiCLBM.F90')
    parser.add_argument('--compiler',type=Path,required=True)
    parser.add_argument('--report',type=Path,default=HERE/'validation_local.json')
    parser.add_argument('--timeout',type=float,default=120)
    args=parser.parse_args();raw=args.source.read_bytes();source=raw.decode('utf-8')
    root=Path(tempfile.mkdtemp(prefix='fei-clbm-local-')).resolve();env=compiler_environment(args.compiler)
    report={'passed':False,'verified_at_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'source_sha256':hashlib.sha256(raw).hexdigest(),'backend':'GNU Fortran OpenACC host; NVIDIA GPU not tested','build_directory':str(root),'high_Ra_validated':False,'grid_independence_validated':False}
    started=time.monotonic()
    try:
        report['four_macro_compilations']={}
        for steady in [True,False]:
            for rb in [True,False]:
                label=('steady' if steady else 'unsteady')+('_rb' if rb else '_side')
                build_run(small_source(source,steady=steady,rb=rb,steps=48),root/label,args.compiler,env,args.timeout,syntax_only=True)
                report['four_macro_compilations'][label]=True
        report['mathematical_checks']={}
        for name,rb,custom in [('side_default',False,False),('rb_default',True,False),('side_nondefault_rates',False,True)]:
            clone=small_source(source,steady=True,rb=rb,steps=48)
            if custom:
                for n,v in [('feiBulkRate','.91d0'),('feiThirdRate','.8d0'),('feiFourthRate','1.23d0'),('feiThermalHighRate','.78d0')]:clone=replace_parameter(clone,n,v)
            output=build_run(replace_main(clone,MATH_MAIN),root/name,args.compiler,env,args.timeout)
            assert 'ALL_MATHEMATICAL_CHECKS_PASSED' in output
            metrics={}
            for line in output.splitlines():
                if line.startswith('CHECK '):
                    _,label,residual,limit=line.split();val=float(residual.replace('D','E'))
                    metrics[label]={'max_abs_residual':max(val,metrics.get(label,{}).get('max_abs_residual',0.0)),'tolerance':float(limit.replace('D','E'))}
            report['mathematical_checks'][name]=metrics
        report['integration_checks']=integration_checks(source,root,args.compiler,env,args.timeout)
        report['source_unchanged']=args.source.read_bytes()==raw
        report['passed']=report['source_unchanged']
    except Exception as error:report['failure']=str(error)
    report['elapsed_seconds']=round(time.monotonic()-started,3)
    args.report.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k in ['passed','source_sha256','elapsed_seconds','failure','build_directory']},ensure_ascii=False,indent=2))
    return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
