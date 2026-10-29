"""核验实际 Luo-MRT/Zhang-HR Fortran：碰撞公式、滤波、导热、输出与严格续算。

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
ZHANG_HISTORY = "ZhangHistory2DOpenaccZhangHR.txt"


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
    integer, parameter :: cx(0:8)=[0,1,0,-1,0,1,-1,-1,1]
    integer, parameter :: cy(0:8)=[0,0,1,0,-1,1,1,-1,-1]
    real(kind=8), parameter :: wi(0:8)=[4d0/9,1d0/9,1d0/9,1d0/9,1d0/9, &
                                       1d0/36,1d0/36,1d0/36,1d0/36]
    real(kind=8) :: p(0:8),q(0:8),eq(0:8),neq(0:8),oracle(0:8),noneq(0:8)
    real(kind=8) :: ux,uy,scalar,cu,gradx,grady,dtx,dty,aprX,aprY,afdX,afdY,ahX,ahY
    real(kind=8) :: beforeMass,beforeJx,beforeJy,forceY,initialHeat,expected,err,secondMoment
    real(kind=8) :: left,right,down,upward,el,er,ed,eu,center,filtered,sigma,sigmaOracle
    real(kind=8) :: raw(nx,ny),savedg(nx,ny,0:8),filterOracle(nx,ny),sigmaField(nx,ny),rangeValue
    call acc_init(acc_device_host)
    call initial()
    call require(all(ex==cx).and.all(ey==cy),'D2Q9_direction_order')
    call require(maxval(abs(omegaT-wi))<1d-15,'D2Q9_weight_definition')
    ! Eq22: first-order thermal equilibrium, including signed scalars.
    do n=1,9
      ux=0.11d0*cos(dble(n));uy=0.08d0*sin(dble(n));scalar=-0.9d0+0.2d0*dble(n)
      call zhang_thermal_equilibrium(scalar,ux,uy,p)
      call measure('linear_equilibrium_mass',abs(sum(p)-scalar),5d-14)
      call measure('linear_equilibrium_flux_x',abs(sum(p*dble(cx))-scalar*ux),5d-14)
      call measure('linear_equilibrium_flux_y',abs(sum(p*dble(cy))-scalar*uy),5d-14)
      secondMoment=sum(p*dble(cx*cx))
      call measure('linear_equilibrium_xx',abs(secondMoment-scalar/3d0),5d-14)
    enddo
    ! Independent Eq23/27/29/30 oracle. scalar is explicitly independent of sum(pop)
    ! so the macro-only filter coupling is tested as well as the filter-off branch.
    do n=1,160
      scalar=0.8d0*cos(dble(n)*0.2d0);ux=0.08d0*cos(dble(n));uy=0.06d0*sin(dble(n))
      p=wi*(scalar+0.2d0*sin(dble(n)*0.719d0+[(dble(k)*1.713d0,k=0,8)]))
      if(mod(n,2)==0) scalar=sum(p)
      gradx=0.014d0*sin(dble(n));grady=0.019d0*cos(dble(n))
      dtx=0.007d0*cos(dble(n)*0.3d0);dty=0.009d0*sin(dble(n)*0.8d0)
      do k=0,8
        cu=dble(cx(k))*ux+dble(cy(k))*uy
        eq(k)=scalar*wi(k)*(1d0+3d0*cu)
      enddo
      neq=p-eq;aprX=sum(dble(cx)*neq);aprY=sum(dble(cy)*neq)
      afdX=-taug*gradx/3d0-0.5d0*dtx;afdY=-taug*grady/3d0-0.5d0*dty
      ahX=zhangBeta*aprX+(1d0-zhangBeta)*afdX
      ahY=zhangBeta*aprY+(1d0-zhangBeta)*afdY
      do k=0,8
        if(zhangThermalModel==0) then
          noneq(k)=neq(k)
        elseif(zhangThermalModel==1) then
          noneq(k)=3d0*wi(k)*(dble(cx(k))*aprX+dble(cy(k))*aprY)
        else
          noneq(k)=3d0*wi(k)*(dble(cx(k))*ahX+dble(cy(k))*ahY)
        endif
        oracle(k)=eq(k)+(1d0-1d0/taug)*noneq(k) &
            +(1d0-0.5d0/taug)*3d0*wi(k)*(dble(cx(k))*dtx+dble(cy(k))*dty)
      enddo
      call zhang_thermal_collision(p,scalar,ux,uy,gradx,grady,dtx,dty,q)
      call measure('Eq23_27_29_30_collision',maxval(abs(q-oracle)),3d-13)
      if(zhangThermalModel/=0.or.mod(n,2)==0) then
        call measure('collision_scalar_mass',abs(sum(q)-scalar),3d-13)
      endif
    enddo
    ! Ghost helper: boundary semantics must match each main case on every edge/corner.
    do j=1,ny
      do i=1,nx
        center=-0.17d0;left=0.31d0;right=-0.43d0;down=0.59d0;upward=-0.61d0
        el=left;er=right;ed=down;eu=upward
#ifdef SideHeatedCell
        if(i==1) el=2d0*Thot-center
        if(i==nx) er=2d0*Tcold-center
        if(j==1) ed=center
        if(j==ny) eu=center
#else
        if(i==1) el=center
        if(i==nx) er=center
        if(j==1) ed=2d0*Thot-center
        if(j==ny) eu=2d0*Tcold-center
#endif
        call zhang_boundary_neighbors(i,j,center,left,right,down,upward)
        call measure('physical_wall_ghosts',maxval(abs([left-el,right-er,down-ed,upward-eu])),1d-15)
      enddo
    enddo
    call zhang_filter_value(-0.3d0,-0.3d0,-0.3d0,-0.3d0,-0.3d0,0d0,sigma,filtered)
    call measure('filter_constant_range_zero',max(abs(sigma),abs(filtered+0.3d0)),1d-15)
    center=-0.4d0;left=1d0;right=-0.1d0;down=0.9d0;upward=0.2d0
    sigmaOracle=zhangFilterSigma0*(1d0-exp(-((right-center)**2+(upward-center)**2)/(2d0*zhangFilterEpsilon)**2))
    expected=center-sigmaOracle*(center-0.25d0*(left+right+down+upward))
    call zhang_filter_value(center,left,right,down,upward,2d0,sigma,filtered)
    call measure('Eq32_dynamic_sigma',abs(sigma-sigmaOracle),3d-14)
    call measure('Eq31_filter_value',abs(filtered-expected),3d-14)
    ! Flow forcing is checked from direct moments, without relying on thermal equilibrium.
    p=wi*(1d0+0.05d0*[(sin(dble(k+1)),k=0,8)]);p=p/sum(p)
    beforeMass=sum(p);beforeJx=sum(p*dble(cx));beforeJy=sum(p*dble(cy))
    rho=beforeMass;T=Tref+0.2d0;forceY=beforeMass*gBeta*0.2d0
    u=beforeJx/beforeMass;v=(beforeJy+0.5d0*forceY)/beforeMass
    do k=0,8
      f(:,:,k)=p(k)
    enddo
    call enter_data_2d_openacc()
    call collision()
    !$acc wait(1)
    !$acc update self(f_post)
    q=f_post(1,1,:)
    call measure('Luo_collision_mass',abs(sum(q)-beforeMass),5d-14)
    call measure('Luo_collision_x_momentum',abs(sum(q*dble(cx))-beforeJx),5d-14)
    call measure('Luo_collision_y_force_increment',abs(sum(q*dble(cy))-beforeJy-forceY),5d-14)
    call require(abs(forceY)>1d-10,'force_test_is_nontrivial')
    ! Wrapper integration: independently construct central gradients and one-step Cu
    ! differences on a non-linear manufactured state, including every wall/corner.
    do j=1,ny
      do i=1,nx
        T(i,j)=0.3d0*sin(0.7d0*dble(i))+0.17d0*cos(0.6d0*dble(j))
        u(i,j)=0.02d0*cos(0.3d0*dble(i));v(i,j)=0.03d0*sin(0.4d0*dble(j))
        zhangFluxPreviousX(i,j)=0.004d0*sin(dble(i+j))
        zhangFluxPreviousY(i,j)=0.003d0*cos(dble(i-j))
        g(i,j,:)=wi*T(i,j)
      enddo
    enddo
    do j=1,ny
      do i=1,nx
        center=T(i,j);left=T(max(1,i-1),j);right=T(min(nx,i+1),j)
        down=T(i,max(1,j-1));upward=T(i,min(ny,j+1))
#ifdef SideHeatedCell
        if(i==1) left=2d0*Thot-center
        if(i==nx) right=2d0*Tcold-center
        if(j==1) down=center
        if(j==ny) upward=center
#else
        if(i==1) left=center
        if(i==nx) right=center
        if(j==1) down=2d0*Thot-center
        if(j==ny) upward=2d0*Tcold-center
#endif
        gradx=0.5d0*(right-left);grady=0.5d0*(upward-down)
        dtx=T(i,j)*u(i,j)-zhangFluxPreviousX(i,j)
        dty=T(i,j)*v(i,j)-zhangFluxPreviousY(i,j)
        p=g(i,j,:)
        call zhang_thermal_collision(p,T(i,j),u(i,j),v(i,j),gradx,grady,dtx,dty,q)
        savedg(i,j,:)=q
      enddo
    enddo
    !$acc update device(T,u,v,g,zhangFluxPreviousX,zhangFluxPreviousY)
    call collisionT()
    !$acc wait(1)
    !$acc update self(g_post,zhangFluxPreviousX,zhangFluxPreviousY)
    call measure('wrapper_central_gradient_and_dt_Cu',maxval(abs(g_post(1:nx,1:ny,:)-savedg)),3d-14)
    call measure('dt_Cu_history_updated_after_collision_x',maxval(abs(zhangFluxPreviousX-T*u)),3d-14)
    call measure('dt_Cu_history_updated_after_collision_y',maxval(abs(zhangFluxPreviousY-T*v)),3d-14)
    ! Signed-temperature conduction, actual production collision/streaming/walls/macros/Nu.
    u=0d0;v=0d0;zhangFluxPreviousX=0d0;zhangFluxPreviousY=0d0
    do j=1,ny
      do i=1,nx
#ifdef SideHeatedCell
        T(i,j)=Thot+(dble(i)-0.5d0)/dble(nx)*(Tcold-Thot)
#else
        T(i,j)=Thot+(dble(j)-0.5d0)/dble(ny)*(Tcold-Thot)
#endif
        g(i,j,:)=wi*T(i,j)
      enddo
    enddo
    initialHeat=sum(T)
    !$acc update device(u,v,T,g,zhangFluxPreviousX,zhangFluxPreviousY)
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
    ! macroT must filter only T, leave g untouched, and accumulate the reported changes.
    do j=1,ny
      do i=1,nx
        raw(i,j)=0.3d0*sin(0.7d0*dble(i))+0.17d0*cos(0.6d0*dble(j))
        g(i,j,:)=wi*raw(i,j)
      enddo
    enddo
    savedg=g;raw=sum(g,dim=3);rangeValue=maxval(raw)-minval(raw)
    do j=1,ny
      do i=1,nx
        center=raw(i,j);left=raw(max(1,i-1),j);right=raw(min(nx,i+1),j)
        down=raw(i,max(1,j-1));upward=raw(i,min(ny,j+1))
#ifdef SideHeatedCell
        if(i==1) left=2d0*Thot-center
        if(i==nx) right=2d0*Tcold-center
        if(j==1) down=center
        if(j==ny) upward=center
#else
        if(i==1) left=center
        if(i==nx) right=center
        if(j==1) down=2d0*Thot-center
        if(j==ny) upward=2d0*Tcold-center
#endif
        sigma=0d0
        if(zhangEnableFilter==1) sigma=zhangFilterSigma0*(1d0-exp( &
          -((right-center)**2+(upward-center)**2)/(zhangFilterEpsilon*rangeValue)**2))
        sigmaField(i,j)=sigma
        filterOracle(i,j)=center-sigma*(center-0.25d0*(left+right+down+upward))
      enddo
    enddo
    zhangFilterAccum=0d0;zhangFilterAbsAccum=0d0
    !$acc update device(g,zhangFilterAccum,zhangFilterAbsAccum)
    call macroT()
    !$acc wait(1)
    !$acc update self(T,g,zhangSigmaLast,zhangFilterAccum,zhangFilterAbsAccum)
    call measure('macro_only_filter_keeps_g',maxval(abs(g-savedg)),0d0)
    call measure('macro_filter_stencil',maxval(abs(T-filterOracle)),3d-14)
    call measure('macro_filter_sigma',maxval(abs(zhangSigmaLast-sigmaField)),3d-14)
    call measure('filter_signed_accumulation',maxval(abs(zhangFilterAccum-(filterOracle-raw))),3d-14)
    call measure('filter_abs_accumulation',maxval(abs(zhangFilterAbsAccum-abs(filterOracle-raw))),3d-14)
    call exit_data_2d_openacc()
    write(*,'(a)') 'ALL_MATHEMATICAL_CHECKS_PASSED'
contains
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


def inspect_zhang_history(directory: Path, meta: dict[str, str]) -> int:
    lines = (directory / ZHANG_HISTORY).read_text().splitlines()
    if not lines or not lines[0].startswith("# itc t_ff sigma_max "):
        raise AssertionError("Zhang history header is missing")
    rows = [line.split() for line in lines[1:] if line.strip()]
    if len(rows) != int(meta["nureSampleCount"]):
        raise AssertionError("Zhang history count differs from Nu/Re/checkpoint")
    total_step = int(meta["itc_total"])
    time_unit = total_step / float(meta["time_tf"])
    previous_step, previous_time, previous_abs = -1, -1.0, -1.0
    for row in rows:
        if len(row) != 10:
            raise AssertionError("Invalid Zhang history column count")
        step = int(row[0])
        values = [float(value) for value in row[1:]]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("Nonfinite Zhang history")
        time_tf, sigma_max, signed_sum, abs_sum, raw_heat, heat, t_min, t_max, mass = values
        if not (previous_step < step <= total_step and time_tf > previous_time):
            raise AssertionError("Zhang history steps/times are not ordered within checkpoint")
        if abs(time_tf - step / time_unit) > 1e-12 * max(1.0, abs(time_tf)):
            raise AssertionError("Zhang history step and time disagree")
        if not (0.0 <= sigma_max <= 1.0 and previous_abs <= abs_sum and abs_sum >= abs(signed_sum)-1e-10):
            raise AssertionError("Invalid cumulative Zhang filter diagnostics")
        if t_min > t_max or abs(mass / (NX*NY) - 1.0) > 1e-11:
            raise AssertionError("Invalid Zhang temperature range/flow mass")
        previous_step, previous_time, previous_abs = step, time_tf, abs_sum
    return len(rows)


def rejected_zhang_restarts(seed: Path, executable: Path, root: Path,
                              env: dict[str, str], timeout: float) -> dict:
    """破坏真实32步检查点配套历史，验证生产读取器拒绝而非继续拼接。"""
    zhang_bytes = (seed / ZHANG_HISTORY).read_bytes()
    lines = zhang_bytes.decode().splitlines(keepends=True)
    row_indices = [index for index, line in enumerate(lines)
                   if line.strip() and not line.lstrip().startswith("#")]
    if not row_indices:
        raise AssertionError("Restart rejection seed must contain zhang samples")
    _, meta = metadata(seed)
    total_step = int(meta["itc_total"])
    time_unit = total_step / float(meta["time_tf"])
    cases = {}
    for label in ("missing", "short", "ahead", "duplicate_tail", "wrong_time", "wrong_step"):
        directory = root / label
        shutil.copytree(seed, directory)
        path = directory / ZHANG_HISTORY
        changed = list(lines)
        if label == "missing":
            path.unlink()
        else:
            if label == "short":
                del changed[row_indices[-1]]
            elif label == "duplicate_tail":
                changed.append(lines[row_indices[-1]])
            elif label == "ahead":
                row = lines[row_indices[-1]].split()
                row[0] = str(total_step + 1)
                row[1] = format((total_step + 1) / time_unit, ".16E")
                changed.append(" ".join(row) + "\n")
            else:
                index = row_indices[0]
                row = changed[index].split()
                if label == "wrong_time":
                    row[1] = format(float(row[1]) + 1.0 / time_unit, ".16E")
                else:
                    row[0] = str(int(row[0]) + 1)
                changed[index] = " ".join(row) + "\n"
            path.write_text("".join(changed), encoding="utf-8")
        protected = {p.name: p.read_bytes() for p in directory.iterdir()
                     if p.suffix == ".bin" or p.name.endswith("-latest.meta")
                     or p.name == ZHANG_HISTORY
                     or p.name.startswith(("Nu_VolAvg", "Re_VolAvg"))}
        result = run_command([str(executable)], directory, env, timeout, expect_success=False)
        output = result.stdout + result.stderr
        (directory / "rejected_restart.log").write_text(output, encoding="utf-8")
        if result.returncode == 0 or "ZhangHR restart diagnostics" not in output:
            raise AssertionError(f"Zhang {label} restart was not refused by the production validator: "
                                 + output[-4000:])
        if any((directory / filename).read_bytes() != content
               for filename, content in protected.items()):
            raise AssertionError(f"Rejected zhang {label} restart changed saved state/history")
        if label == "missing" and path.exists():
            raise AssertionError("Rejected missing zhang history was silently recreated")
        cases[label] = {"rejected_by_zhang_validator": True,
                        "checkpoint_and_histories_preserved": True}
    return cases


def inspect_run(directory: Path, expected_magic: str, *, filter_on: bool = False) -> dict:
    state, meta = metadata(directory)
    records = read_records(state)
    if len(records[0]) != 48 or records[0][:32].rstrip().decode() != expected_magic:
        raise AssertionError("Wrong scheme identifier in checkpoint")
    if struct.unpack("<4i", records[0][32:]) != (NX, NY, 9, 9):
        raise AssertionError("Wrong checkpoint mesh or velocity counts")
    if any(len(r) != NX * NY * 9 * 8 for r in records[1:3]):
        raise AssertionError("Checkpoint does not contain complete f/g")
    flow, thermal = [struct.unpack(f"<{NX*NY*9}d", r) for r in records[1:3]]
    if not all(math.isfinite(v) for v in (*flow, *thermal)):
        raise AssertionError("Nonfinite f/g")
    if len(records[3]) != 6 * NX * NY * 8:
        raise AssertionError("Checkpoint must retain T, two flux histories and three filter fields")
    states = struct.unpack(f"<{6*NX*NY}d", records[3])
    if not all(math.isfinite(value) for value in states):
        raise AssertionError("Nonfinite macro/filter/history checkpoint state")
    field_size = NX * NY
    saved_temperature = states[:field_size]
    raw_temperature = [sum(thermal[k*field_size+i] for k in range(9)) for i in range(field_size)]
    filter_delta = max(abs(a-b) for a, b in zip(saved_temperature, raw_temperature))
    if not filter_on and filter_delta > 1e-14:
        raise AssertionError("Unfiltered checkpoint T differs from sum(g)")
    if filter_on and (filter_delta <= 1e-14 or max(states[4*field_size:5*field_size]) <= 0):
        raise AssertionError("Filtered strict-restart case did not exercise an independent T/filter state")
    mass_error = abs(sum(flow) - NX * NY) / (NX * NY)
    if mass_error > 1e-11:
        raise AssertionError(f"Closed-cavity mass error: {mass_error}")
    snapshots = sorted(directory.glob("*Snapshot-*.bin"))
    if not snapshots:
        raise AssertionError("No post-processing snapshots")
    for snapshot in snapshots:
        fields = read_records(snapshot)
        if [len(r) for r in fields] != [NX * NY * 8] * 4:
            raise AssertionError("Snapshot must be four u/v/T/rho fields")
        if not all(math.isfinite(x) for r in fields for x in struct.unpack(f"<{NX*NY}d", r)):
            raise AssertionError("Nonfinite macro snapshot")
    counts = {}
    for prefix in ("Nu_VolAvg", "Re_VolAvg"):
        paths = list(directory.glob(prefix + "*.dat"))
        if len(paths) != 1:
            raise AssertionError(f"Expected one {prefix} series")
        rows = [tuple(map(float, line.split())) for line in paths[0].read_text().splitlines()
                if line.strip() and not line.lstrip().startswith("#")]
        if not rows or any(len(r) != 2 or not all(math.isfinite(x) for x in r) for r in rows):
            raise AssertionError("Nu/Re contains invalid samples")
        if any(a[0] >= b[0] for a, b in zip(rows, rows[1:])):
            raise AssertionError("Nu/Re times are not strictly increasing")
        counts[prefix] = len(rows)
    if any(n != int(meta["nureSampleCount"]) for n in counts.values()):
        raise AssertionError("Saved sample count differs from Nu/Re history")
    return {"itc_total": int(meta["itc_total"]), "relative_mass_error": mass_error,
            "finite_Nu_Re_samples": counts, "snapshot_count": len(snapshots),
            "zhang_history_samples": inspect_zhang_history(directory, meta),
            "maximum_filtered_T_minus_sum_g": filter_delta,
            "checkpoint_record_bytes": [len(r) for r in records]}


def integration_checks(source: str, root: Path, compiler: Path,
                       env: dict[str, str], timeout: float) -> dict:
    match = re.search(r"reloadMagic\s*=\s*['\"]([^'\"]+)['\"]", source, re.I)
    if not match:
        raise AssertionError("Scheme-specific reloadMagic missing")
    magic = match[1]
    results = {}
    for name, steady, rb, filter_on in (("side_steady", True, False, False),
                                        ("rb_unsteady", False, True, False),
                                        ("rb_unsteady_filtered", False, True, True)):
        case_source = replace_parameter(source, "zhangEnableFilter", "1" if filter_on else "0")
        full = root / name
        build_run(small_source(case_source, steady=steady, rb=rb, steps=48), full, compiler, env, timeout)
        results[name] = inspect_run(full, magic, filter_on=filter_on)
        split = root / (name + "_restart")
        build_run(small_source(case_source, steady=steady, rb=rb, steps=32), split, compiler, env, timeout)
        seed = root / (name + "_restart_seed32")
        shutil.copytree(split, seed)
        histories = {p.name: p.read_bytes() for prefix in ("Nu_VolAvg", "Re_VolAvg")
                     for p in split.glob(prefix + "*.dat")}
        histories[ZHANG_HISTORY] = (split / ZHANG_HISTORY).read_bytes()
        build_run(small_source(case_source, steady=steady, rb=rb, steps=48, restart=True),
                  split, compiler, env, timeout)
        after = inspect_run(split, magic, filter_on=filter_on)
        full_state, fm = metadata(full)
        split_state, sm = metadata(split)
        if int(fm["itc_total"]) != 48 or int(sm["itc_total"]) != 48:
            raise AssertionError("Strict restart changed the total horizon")
        if full_state.read_bytes() != split_state.read_bytes():
            raise AssertionError(f"{name}: continuous/split checkpoint bytes differ")
        for filename, content in histories.items():
            current = (split / filename).read_bytes()
            if not current.startswith(content):
                raise AssertionError(f"{name}: restart overwrote history")
            if current.splitlines()[-1] != (full / filename).read_bytes().splitlines()[-1]:
                raise AssertionError(f"{name}: final history row differs: {filename}")
        if not steady and (split / ZHANG_HISTORY).read_bytes() != (full / ZHANG_HISTORY).read_bytes():
            raise AssertionError("RB continuous/split zhang history bytes differ")
        full_snapshots = {p.name: p.read_bytes() for p in full.glob("*Snapshot-*.bin")}
        split_snapshots = {p.name: p.read_bytes() for p in split.glob("*Snapshot-*.bin")}
        if (not steady and full_snapshots.keys() != split_snapshots.keys()) or any(
                split_snapshots.get(n) != data for n, data in full_snapshots.items()):
            raise AssertionError(f"{name}: strict restart snapshot timeline differs")
        protected = {p.name: p.read_bytes() for p in split.iterdir()
                     if p.name in histories or p.name.endswith("-latest.meta")}
        run_command([str(split / "solver.exe")], split, env, timeout)
        if any((split / n).read_bytes() != data for n, data in protected.items()):
            raise AssertionError("Completed-target restart changed checkpoint/history")
        results[name + "_strict_restart"] = {"checkpoint_bitwise_equal": True,
            "snapshot_bitwise_equal": True, "Nu_Re_history_preserved": True,
            "zhang_history_preserved": True, "zhang_final_row_bitwise_equal": True,
            "zhang_full_history_bitwise_equal": not steady,
            "zhang_history_comparison_note": ("Steady stop at32 adds a final sample; compare old prefix and final48 row"
                                                if steady else "Every RB sample and header are bitwise equal"),
            "completed_target_restart_preserved": True, "after": after}
        results[name + "_invalid_zhang_restarts"] = rejected_zhang_restarts(
            seed, split / "solver.exe", root / (name + "_invalid_zhang"), env, timeout)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "2DRBOpenaccZhangHR.F90")
    parser.add_argument("--compiler", type=Path)
    parser.add_argument("--report", type=Path, default=HERE / "validation_local.json")
    parser.add_argument("--temp-root", type=Path)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--math-only", action="store_true")
    args = parser.parse_args()
    compiler = args.compiler or Path(shutil.which("gfortran") or "C:/msys64/ucrt64/bin/gfortran.exe")
    if not compiler.is_file():
        parser.error(f"Compiler not found: {compiler}")
    source_bytes = args.source.read_bytes()
    source = source_bytes.decode("utf-8-sig").replace("\r\n", "\n")
    if args.temp_root:
        args.temp_root.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="zhang-hr-local-", dir=args.temp_root)).resolve()
    if not str(root).isascii():
        parser.error("Use an ASCII --temp-root for portable Fortran builds")
    env = compiler_environment(compiler)
    report = {"passed": False, "verified_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "source": str(args.source.resolve()), "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
              "compiler": str(compiler), "backend": "GNU Fortran OpenACC host; NVIDIA GPU not tested",
              "build_directory": str(root), "high_Ra_accuracy_validated": False,
              "grid_independence_validated": False,
              "scope": "Actual production kernels; independent HR/forcing/filter/conduction and finite-step output/restart checks"}
    started = time.monotonic()
    try:
        report["compiler_version"] = run_command([str(compiler), "--version"], root, env, args.timeout).stdout.splitlines()[0]
        report["four_macro_compilations"] = {}
        for steady in (True, False):
            for rb in (False, True):
                label = ("steady_" if steady else "unsteady_") + ("rb" if rb else "side")
                build_run(small_source(source, steady=steady, rb=rb, steps=48), root / ("syntax_" + label),
                          compiler, env, args.timeout, syntax_only=True)
                report["four_macro_compilations"][label] = True
        report["mathematical_checks"] = {}
        for label, rb, model, beta, filter_on in (
                ("HR_side", False, 2, ".95d0", 0), ("HR_RB", True, 2, ".95d0", 0),
                ("BGK", False, 0, ".95d0", 0), ("PR1", False, 1, ".95d0", 0),
                ("HR_beta0", False, 2, "0d0", 0), ("HR_beta1", False, 2, "1d0", 0),
                ("HR_filter", False, 2, ".95d0", 1)):
            clone = small_source(source, steady=True, rb=rb, steps=48)
            clone = replace_parameter(clone, "zhangThermalModel", str(model))
            clone = replace_parameter(clone, "zhangBeta", beta)
            clone = replace_parameter(clone, "zhangEnableFilter", str(filter_on))
            output = build_run(replace_main(clone, MATH_MAIN), root / ("math_" + label), compiler, env, args.timeout)
            if "ALL_MATHEMATICAL_CHECKS_PASSED" not in output:
                raise AssertionError("Mathematical driver did not finish")
            metrics = {}
            for line in output.splitlines():
                if line.startswith("CHECK "):
                    _, name, residual, limit = line.split()
                    residual = float(residual.replace("D", "E"))
                    metrics[name] = {"max_abs_residual": max(residual, metrics.get(name, {}).get("max_abs_residual", 0.0)),
                                     "tolerance": float(limit.replace("D", "E"))}
            report["mathematical_checks"][label] = metrics
        if not args.math_only:
            report["integration_checks"] = integration_checks(source, root, compiler, env, args.timeout)
        report["source_unchanged"] = args.source.read_bytes() == source_bytes
        report["passed"] = report["source_unchanged"]
    except Exception as error:
        report["failure"] = str(error)
    finally:
        report["elapsed_seconds"] = round(time.monotonic() - started, 3)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({k: report[k] for k in ("passed", "source_sha256", "build_directory", "elapsed_seconds", "failure")
                          if k in report}, ensure_ascii=False, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
