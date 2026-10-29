"""核验实际 Wei 双BGK Fortran：碰撞公式、显式浮力、导热、输出与严格续算。

只在隔离临时目录中替换生产源码的测试参数或主程序，生产文件保持不变。
需要支持 OpenACC 的 GNU Fortran；OpenACC host 不代表 NVIDIA GPU 验收。
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent
NX, NY = 12, 10


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



MATH_MAIN = r"""
program main
  use commondata
  use openacc
  implicit none
  integer :: i,j,k,n
  integer,parameter :: cx(0:8)=[0,1,0,-1,0,1,-1,-1,1],cy(0:8)=[0,0,1,0,-1,1,1,-1,-1]
  real(kind=8),parameter :: wi(0:8)=[4d0/9,1d0/9,1d0/9,1d0/9,1d0/9,1d0/36,1d0/36,1d0/36,1d0/36]
  real(kind=8) :: p(0:8),q(0:8),eq(0:8),oracle(0:8),ux,uy,scalar,cu,fytest,rhotest,initialHeat,expected,err
  call acc_init(acc_device_host)
  call initial()
  call measure('D2Q9_weights',maxval(abs(omegaT-wi)),1d-15)
  do n=1,40
    ux=.09d0*cos(dble(n));uy=.07d0*sin(dble(n));scalar=.1d0+.02d0*dble(n)
    call wei_equilibrium(scalar,ux,uy,q)
    do k=0,8
      cu=dble(cx(k))*ux+dble(cy(k))*uy
      eq(k)=wi(k)*scalar*(1d0+3d0*cu+4.5d0*cu**2-1.5d0*(ux**2+uy**2))
      p(k)=wi(k)*scalar*(1d0+.3d0*sin(dble(k+n)))
    enddo
    call measure('polynomial_equilibrium',maxval(abs(q-eq)),2d-15)
    call measure('equilibrium_mass',abs(sum(q)-scalar),2d-15)
    call measure('equilibrium_jx',abs(sum(q*dble(cx))-scalar*ux),2d-15)
    fytest=.00013d0
    oracle=p+(eq-p)/tauf+3d0*wi*fytest*dble(cy)
    call wei_flow_collision(p,scalar,ux,uy,fytest,q)
    call measure('Eq1_BGK_and_Eq3_force',maxval(abs(q-oracle)),2d-15)
    oracle=p+(eq-p)/taug
    call wei_thermal_collision(p,scalar,ux,uy,q)
    call measure('Eq2_thermal_BGK',maxval(abs(q-oracle)),2d-15)
  enddo
  p=wi*(1d0+.05d0*[(sin(dble(k+1)),k=0,8)]);p=p/sum(p)
  rhotest=sum(p);ux=sum(p*dble(cx))/rhotest;uy=sum(p*dble(cy))/rhotest
  rho=rhotest;u=ux;v=uy;T=.3d0;fytest=rhotest*gBeta*.3d0
  do k=0,8
    f(:,:,k)=p(k)
  enddo
  call enter_data_2d_openacc()
  call collision()
  !$acc wait(1)
  !$acc update self(f_post)
  q=f_post(1,1,:)
  call measure('production_flow_mass',abs(sum(q)-rhotest),5d-14)
  call measure('production_force_jx',abs(sum(q*dble(cx))-rhotest*ux),5d-14)
  call measure('production_force_jy',abs(sum(q*dble(cy))-rhotest*uy-fytest),5d-14)
  if(abs(fytest)<1d-10) error stop 'Nontrivial force test required'
  call macro()
  !$acc wait(1)
  !$acc update self(u,v)
  call measure('raw_momentum_u',abs(u(1,1)-ux),5d-15)
  call measure('raw_momentum_v_NO_half_force',abs(v(1,1)-uy),5d-15)
  u=0d0;v=0d0
  do j=1,ny
    do i=1,nx
#ifdef SideHeatedCell
      T(i,j)=Thot+(dble(i)-.5d0)/dble(nx)*(Tcold-Thot)
#else
      T(i,j)=Thot+(dble(j)-.5d0)/dble(ny)*(Tcold-Thot)
#endif
      g(i,j,:)=wi*T(i,j)
    enddo
  enddo
  initialHeat=sum(T)
  !$acc update device(u,v,T,g)
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
      expected=Thot+(dble(i)-.5d0)/dble(nx)*(Tcold-Thot)
#else
      expected=Thot+(dble(j)-.5d0)/dble(ny)*(Tcold-Thot)
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
  subroutine measure(label,residual,limit)
    character(len=*),intent(in) :: label
    real(kind=8),intent(in) :: residual,limit
    write(*,'(a,1x,a,2(1x,es24.16e3))') 'CHECK',label,residual,limit
    if(.not.(residual<=limit)) error stop 42
  end subroutine measure
end program main
"""
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



def inspect_run(directory: Path) -> dict:
    state,meta=metadata(directory)
    records=read_records(state)
    if records[0][:32].rstrip().decode()!='WEI_BGK_2D_V1':
        raise AssertionError('Wrong scheme header')
    if struct.unpack('<4i',records[0][32:])!=(NX,NY,9,9):
        raise AssertionError('Wrong checkpoint grid/Q')
    if any(len(r)!=NX*NY*9*8 for r in records[1:3]):
        raise AssertionError('Incomplete f/g populations')
    flow,thermal=[struct.unpack(f'<{NX*NY*9}d',r) for r in records[1:3]]
    if not all(math.isfinite(x) for x in (*flow,*thermal)):
        raise AssertionError('Nonfinite f/g')
    mass_error=abs(sum(flow)-NX*NY)/(NX*NY)
    if mass_error>1e-11:raise AssertionError(f'Flow mass drift {mass_error}')
    for prefix in ('Nu_VolAvg','Re_VolAvg'):
        paths=list(directory.glob(prefix+'*.dat'))
        if len(paths)!=1:raise AssertionError('Missing Nu/Re series')
        rows=[tuple(map(float,l.split())) for l in paths[0].read_text().splitlines() if l.strip()]
        if len(rows)!=int(meta['nureSampleCount']):raise AssertionError('History count mismatch')
        if any(len(r)!=2 or not all(math.isfinite(x) for x in r) for r in rows):raise AssertionError('Invalid history')
        if any(a[0]>=b[0] for a,b in zip(rows,rows[1:])):raise AssertionError('Unordered history')
    snapshots=list(directory.glob('*Snapshot-*.bin'))
    for path in snapshots:
        rs=read_records(path)
        if [len(r) for r in rs]!=[NX*NY*8]*4:raise AssertionError('Wrong snapshot shape')
    return {'itc_total':int(meta['itc_total']),'relative_mass_error':mass_error,'record_sizes':[len(r) for r in records], 'snapshots':len(snapshots)}


def integration_checks(source,root,compiler,env,timeout):
    report={}
    for name,steady,rb in [('side_steady',True,False),('rb_unsteady',False,True)]:
        full=root/name;split=root/(name+'_restart')
        build_run(small_source(source,steady=steady,rb=rb,steps=48),full,compiler,env,timeout)
        report[name]=inspect_run(full)
        build_run(small_source(source,steady=steady,rb=rb,steps=32),split,compiler,env,timeout)
        seed=root/(name+'_seed');shutil.copytree(split,seed)
        old={p.name:p.read_bytes() for p in split.glob('*VolAvg*.dat')}
        build_run(small_source(source,steady=steady,rb=rb,steps=48,restart=True),split,compiler,env,timeout)
        inspect_run(split)
        if metadata(full)[0].read_bytes()!=metadata(split)[0].read_bytes():raise AssertionError(name+' checkpoint mismatch')
        for n,data in old.items():
            if not (split/n).read_bytes().startswith(data):raise AssertionError('History overwritten')
            if (split/n).read_bytes().splitlines()[-1]!=(full/n).read_bytes().splitlines()[-1]:raise AssertionError('Final Nu/Re row mismatch')
        snapshots={p.name:p.read_bytes() for p in full.glob('*Snapshot-*.bin')}
        other={p.name:p.read_bytes() for p in split.glob('*Snapshot-*.bin')}
        if any(other.get(n)!=x for n,x in snapshots.items()):raise AssertionError('Snapshot mismatch')
        protected={p.name:p.read_bytes() for p in split.iterdir() if p.suffix in ('.bin','.meta','.dat')}
        run_command([str(split/'solver.exe')],split,env,timeout)
        if any((split/n).read_bytes()!=x for n,x in protected.items()):raise AssertionError('Completed restart changed state')
        invalid={}
        for label in ['missing_nu','ahead_nu','bad_time','bad_parameters','bad_binary']:
            d=root/(name+'_'+label);shutil.copytree(seed,d)
            if label=='missing_nu':(d/'Nu_VolAvg_2DOpenaccWeiBGK.dat').unlink()
            elif label=='ahead_nu':
                p=d/'Nu_VolAvg_2DOpenaccWeiBGK.dat';p.write_text(p.read_text()+p.read_text().splitlines()[-1]+'\n')
            elif label=='bad_time':
                p=d/'Nu_VolAvg_2DOpenaccWeiBGK.dat';lines=p.read_text().splitlines();row=lines[0].split();row[0]='999999';lines[0]=' '.join(row);p.write_text('\n'.join(lines)+'\n')
            elif label=='bad_parameters':
                p=list(d.glob('*latest.meta'))[0];lines=p.read_text().splitlines();lines=[('physicalParameters '+' '.join(['9.9e9']+x.split()[2:])) if x.startswith('physicalParameters') else x for x in lines];p.write_text('\n'.join(lines)+'\n')
            else:
                p=metadata(d)[0];data=bytearray(p.read_bytes());data[4]=ord('X');p.write_bytes(data)
            r=run_command([str(split/'solver.exe')],d,env,timeout,expect_success=False)
            if r.returncode==0 or 'WeiBGK' not in r.stdout+r.stderr:raise AssertionError('Bad restart accepted '+label)
            invalid[label]=True
        report[name+'_strict_restart']={'checkpoint_bitwise_equal':True,'shared_snapshots_bitwise_equal':True,'history_preserved':True,'completed_horizon_preserved':True,'invalid_restarts_rejected':invalid}
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=HERE/'2DRBOpenaccWeiBGK.F90')
    parser.add_argument('--compiler',type=Path,required=True)
    parser.add_argument('--report',type=Path,default=HERE/'validation_local.json')
    parser.add_argument('--timeout',type=float,default=120)
    args=parser.parse_args();raw=args.source.read_bytes();source=raw.decode()
    root=Path(tempfile.mkdtemp(prefix='wei-bgk-local-')).resolve();env=compiler_environment(args.compiler)
    report={'passed':False,'source_sha256':hashlib.sha256(raw).hexdigest(),'backend':'GNU Fortran OpenACC host; not NVIDIA GPU','build_directory':str(root),'high_Ra_validated':False,'grid_independence_validated':False}
    started=time.monotonic()
    try:
        report['four_macro_compilations']={}
        for steady in [True,False]:
            for rb in [False,True]:
                name=('steady' if steady else 'unsteady')+('_rb' if rb else '_side')
                build_run(small_source(source,steady=steady,rb=rb,steps=48),root/name,args.compiler,env,args.timeout,syntax_only=True)
                report['four_macro_compilations'][name]=True
        report['mathematical_checks']={}
        for name,rb in [('side',False),('rb',True)]:
            output=build_run(replace_main(small_source(source,steady=True,rb=rb,steps=48),MATH_MAIN),root/('math_'+name),args.compiler,env,args.timeout)
            if 'ALL_MATHEMATICAL_CHECKS_PASSED' not in output:raise AssertionError('Math driver unfinished')
            metrics={}
            for line in output.splitlines():
                if line.startswith('CHECK '):
                    _,label,value,limit=line.split();value=float(value);metrics[label]={'max_abs_residual':max(value,metrics.get(label,{}).get('max_abs_residual',0)),'tolerance':float(limit)}
            report['mathematical_checks'][name]=metrics
        report['integration_checks']=integration_checks(source,root,args.compiler,env,args.timeout)
        report['passed']=True
    except Exception as ex:report['failure']=str(ex)
    report['source_unchanged']=raw==args.source.read_bytes();report['passed'] &= report['source_unchanged']
    report['elapsed_seconds']=time.monotonic()-started;args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['passed','source_sha256','elapsed_seconds','failure'] if k in report},indent=2))
    return 0 if report['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
