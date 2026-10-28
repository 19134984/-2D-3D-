"""验证实际混合 RLB4 + RK4 源码，不修改参数原件、不访问服务器。

需要 Python 3、NumPy 和 gfortran；使用 OpenACC host，所有构建/运行都在 ASCII 临时目录完成。
检查矩约束、热方程离散算子和 RK4 模态推进，再检查小网格输出与严格续算。
这些检查不构成 NVIDIA GPU、高 Ra 精度或网格独立性验证。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import tempfile
import time

import numpy as np


HERE = Path(__file__).resolve().parent
NX, NY = 12, 10


def replace_parameter(source: str, name: str, value: str) -> str:
    pattern = rf"(?im)(\b{re.escape(name)}\s*=\s*)[^,!\n]+"
    changed, count = re.subn(pattern, lambda m: m[1] + value, source, count=1)
    if count != 1:
        raise AssertionError(f"Parameter {name} not found")
    return changed


def small_source(source: str, *, steady: bool, rb: bool, steps: int,
                 restart: bool = False) -> str:
    for macro, enabled in {"steadyFlow": steady, "unsteadyFlow": not steady,
                           "SideHeatedCell": not rb, "RayleighBenardCell": rb}.items():
        pattern = rf"(?m)^[ \t]*!?[ \t]*#define[ \t]+{macro}[ \t]*(?:![^\n]*)?$"
        source, count = re.subn(pattern, ("#define " if enabled else "!#define ") + macro, source)
        if count == 0:
            raise AssertionError(f"Case-selection macro {macro} not found")
    for name, value in (("nx", str(NX)), ("ny", str(NY)), ("Rayleigh", "1.0d4"),
                        ("loadInitField", "1" if restart else "0"),
                        ("steadyCheckIntervalItc", "13")):
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


def replace_main(source: str, main: str) -> str:
    source, count = re.subn(r"(?ims)^\s*program\s+main\s*$.*?^\s*end\s+program\s+main\s*$",
                            lambda _: "\n" + main + "\n", source, count=1)
    if count != 1:
        raise AssertionError("Expected one actual program main")
    return source


MATH_MAIN = r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,k,n
    integer,parameter :: cx(0:8)=[0,1,0,-1,0,1,-1,-1,1]
    integer,parameter :: cy(0:8)=[0,0,1,0,-1,1,1,-1,-1]
    integer,parameter :: px(0:8)=[0,1,0,2,0,1,2,1,2]
    integer,parameter :: py(0:8)=[0,0,1,0,2,1,1,2,2]
    real(kind=8),parameter :: eigenX(nx)=[__EIGEN_X__]
    real(kind=8),parameter :: eigenY(ny)=[__EIGEN_Y__]
    real(kind=8),parameter :: decayX=__DECAY_X__,decayY=__DECAY_Y__
    real(kind=8) :: pop(0:8),post(0:8),eq(0:8),mom(0:8),before(0:8),target(0:8)
    real(kind=8) :: stage(0:nx+1,0:ny+1),rhs(nx,ny),baseT(nx,ny),mode(nx,ny)
    real(kind=8) :: rr,ux,uy,forceX,forceY,err,expected,x,y,initialHeat,wallHeat
    real(kind=8) :: lambda,dt,z,amplification,amp,stepGain,fullError,halfError

    call acc_init(acc_device_host)
    call initial()
    call require(all(ex==cx).and.all(ey==cy),'direction_order')

    ! Oracle 是离散速度上的直接求和，独立于展开后的 Hermite 公式。
    do n=1,12
      rr=0.9d0+0.02d0*dble(n)
      ux=0.13d0*cos(dble(n));uy=0.09d0*sin(dble(n))
      call rlb_equilibrium(rr,ux,uy,eq)
      call moment_definition(eq,ux,uy,mom)
      target=0.0d0
      target(0)=rr;target(3)=cs2*rr;target(4)=cs2*rr;target(8)=cs2*cs2*rr
      call measure('fourth_order_equilibrium_central_moments',maxval(abs(mom-target)),4.0d-15)
      call rlb_collide_node(eq,rr,ux,uy,0.0d0,0.0d0,post)
      call measure('equilibrium_fixed_point',maxval(abs(post-eq)),4.0d-15)
    enddo

    ! 论文的宏观速度定义为原始一阶矩，没有 F/2；检查实际 initial() 的约定。
    err=0.0d0
    do j=1,ny
      do i=1,nx
        call moment_definition(f(i,j,:),0.0d0,0.0d0,mom)
        err=max(err,abs(mom(0)-rho(i,j)),abs(mom(1)-rho(i,j)*u(i,j)))
        err=max(err,abs(mom(2)-rho(i,j)*v(i,j)))
      enddo
    enddo
    call measure('initial_raw_mass_momentum',err,3.0d-14)

    ! 加入保持质量和动量的非平衡扰动，检查正则化的守恒和二阶应力松弛。
    rr=1.17d0;ux=0.073d0;uy=-0.041d0
    call rlb_equilibrium(rr,ux,uy,eq)
    pop=eq
    pop(0)=pop(0)-0.004d0
    pop(1)=pop(1)+0.003d0;pop(3)=pop(3)+0.003d0
    pop(2)=pop(2)-0.001d0;pop(4)=pop(4)-0.001d0
    pop(5)=pop(5)+0.0007d0;pop(7)=pop(7)+0.0007d0
    pop(6)=pop(6)-0.0007d0;pop(8)=pop(8)-0.0007d0
    call moment_definition(pop,ux,uy,before)
    call rlb_collide_node(pop,rr,ux,uy,0.0d0,0.0d0,post)
    call moment_definition(post,ux,uy,mom)
    call measure('regularized_zero_force_mass',abs(mom(0)-before(0)),4.0d-14)
    call measure('regularized_zero_force_momentum',maxval(abs(mom(1:2)-before(1:2))),4.0d-14)
    target=0.0d0;target(3)=cs2*rr;target(4)=cs2*rr
    do k=3,5
      expected=target(k)+(1.0d0-rlbRelaxationRate)*(before(k)-target(k))
      call measure('regularized_second_order_stress',abs(mom(k)-expected),4.0d-14)
    enddo
    call measure('regularized_third_order_central_noneq',maxval(abs(mom(6:7))),4.0d-14)
    expected=cs2*cs2*rr+(1.0d0-rlbRelaxationRate)*cs2*(before(3)+before(4)-2.0d0*cs2*rr)
    call measure('regularized_fourth_order_central_closure',abs(mom(8)-expected),4.0d-14)

    forceX=1.3d-4;forceY=-2.1d-4
    call moment_definition(pop,0.0d0,0.0d0,before)
    call rlb_collide_node(pop,rr,ux,uy,forceX,forceY,post)
    call moment_definition(post,0.0d0,0.0d0,mom)
    call measure('forced_mass',abs(mom(0)-before(0)),4.0d-14)
    call measure('forced_x_momentum_increment',abs(mom(1)-before(1)-forceX),4.0d-14)
    call measure('forced_y_momentum_increment',abs(mom(2)-before(2)-forceY),4.0d-14)

    ! 正交的纯高阶非平衡扰动不改变二阶应力，应被正则化清除。
    call rlb_equilibrium(rr,0.0d0,0.0d0,eq)
    do n=1,2
      do k=0,8
        if(n.EQ.1) then
          pop(k)=eq(k)+0.01d0*omega(k)*(dble(cx(k)*cx(k))-cs2)*dble(cy(k))
        else
          pop(k)=eq(k)+0.01d0*omega(k)*(dble(cx(k)*cx(k))-cs2)*(dble(cy(k)*cy(k))-cs2)
        endif
      enddo
      call rlb_collide_node(pop,rr,0.0d0,0.0d0,0.0d0,0.0d0,post)
      call measure('pure_high_order_ghost_removed',maxval(abs(post-eq)),4.0d-15)
    enddo

    ! 均匀场的内点差分为零；边界仍保留实际冷热壁，不人为改成周期边界。
    stage=0.5d0*(Thot+Tcold);u=0.073d0;v=-0.041d0
    call thermal_rhs(stage,rhs)
    !$acc wait(1)
    call measure('uniform_interior_rhs',maxval(abs(rhs(2:nx-1,2:ny-1))),4.0d-14)

    ! 常速下二次多项式的中心二阶导数精确，独立核查两个对流和扩散分量。
    do j=1,ny
      do i=1,nx
        x=dble(i)-0.5d0;y=dble(j)-0.5d0
        stage(i,j)=0.002d0*x*x-0.001d0*x*y+0.003d0*y*y+0.01d0*x-0.02d0*y+1.3d0
      enddo
    enddo
    call thermal_rhs(stage,rhs)
    !$acc wait(1)
    err=0.0d0
    do j=2,ny-1
      do i=2,nx-1
        x=dble(i)-0.5d0;y=dble(j)-0.5d0
        expected=-0.073d0*(0.004d0*x-0.001d0*y+0.01d0) &
                 +0.041d0*(-0.001d0*x+0.006d0*y-0.02d0)+0.01d0*diffusivity
        err=max(err,abs(rhs(i,j)-expected))
      enddo
    enddo
    call measure('quadratic_advection_diffusion_rhs',err,4.0d-14)

    u=0.0d0;v=0.0d0
    amp=0.07d0
    do j=1,ny
      do i=1,nx
        x=(dble(i)-0.5d0)/dble(nx);y=(dble(j)-0.5d0)/dble(ny)
#ifdef SideHeatedCell
        baseT(i,j)=Thot+x*(Tcold-Thot)
        mode(i,j)=eigenX(i)
#else
        baseT(i,j)=Thot+y*(Tcold-Thot)
        mode(i,j)=eigenY(j)
#endif
      enddo
    enddo

    ! 线性纯导热应为离散稳态；其全场及两壁 Nu 均为 1。
    T=baseT;initialHeat=sum(T)
    stage(1:nx,1:ny)=T
    call thermal_rhs(stage,rhs)
    !$acc wait(1)
    call measure('linear_conduction_rhs',maxval(abs(rhs)),4.0d-14)
    do n=1,80
      call advance_temperature_rk4(1.0d0)
    enddo
    !$acc wait(1)
    call measure('linear_conduction_temperature',maxval(abs(T-baseT)),2.0d-12)
    call measure('linear_conduction_heat_balance',abs(sum(T)-initialHeat),2.0d-11)
#ifdef SideHeatedCell
    call SideHeatedcalc_Nu_global()
    call SideHeatedcalc_Nu_wall_avg()
#else
    call RBcalc_Nu_global()
    call RBcalc_Nu_wall_avg()
#endif
    call measure('linear_conduction_Nu_volume',abs(Nu_global-1.0d0),1.0d-10)
    call measure('linear_conduction_Nu_hot',abs(Nu_hot-1.0d0),1.0d-10)
    call measure('linear_conduction_Nu_cold',abs(Nu_cold-1.0d0),1.0d-10)

    ! 独立矩阵求得的 Dirichlet 扩散特征模态；绝热方向取常数，不需要周期扩展。
    T=baseT+amp*mode
    stage(1:nx,1:ny)=T
    call thermal_rhs(stage,rhs)
    !$acc wait(1)
#ifdef SideHeatedCell
    lambda=diffusivity*decayX
#else
    lambda=diffusivity*decayY
#endif
    call measure('thermal_discrete_eigenmode_rhs',maxval(abs(rhs-lambda*amp*mode)),4.0d-14)

    ! 积分守恒 oracle：全域扩散变化等于冷热壁的净导热通量，绝热壁没有通量。
#ifdef SideHeatedCell
    wallHeat=diffusivity*sum((8.0d0*Thot-9.0d0*T(1,:)+T(2,:) &
                             +8.0d0*Tcold-9.0d0*T(nx,:)+T(nx-1,:))/3.0d0)
#else
    wallHeat=diffusivity*sum((8.0d0*Thot-9.0d0*T(:,1)+T(:,2) &
                             +8.0d0*Tcold-9.0d0*T(:,ny)+T(:,ny-1))/3.0d0)
#endif
    call measure('thermal_integral_boundary_flux_balance',abs(sum(rhs)-wallHeat),4.0d-13)
    call require(abs(wallHeat)>1.0d-9,'boundary_flux_test_nontrivial')

    ! RK4 在该线性模态上的精确稳定性多项式，能区分 Euler/RK2/RK3 和阶段时序错误。
    dt=0.2d0/diffusivity
    z=lambda*dt/dble(thermalSubsteps)
    stepGain=1.0d0+z+z*z/2.0d0+z*z*z/6.0d0+z**4/24.0d0
    amplification=stepGain**thermalSubsteps
    call advance_temperature_rk4(dt)
    !$acc wait(1)
    call measure('rk4_full_step_mode_amplification',maxval(abs(T-baseT-amp*amplification*mode)),5.0d-13)
    fullError=maxval(abs(T-baseT-amp*exp(lambda*dt)*mode))

    T=baseT+amp*mode
    z=lambda*dt/(2.0d0*dble(thermalSubsteps))
    stepGain=1.0d0+z+z*z/2.0d0+z*z*z/6.0d0+z**4/24.0d0
    amplification=stepGain**(2*thermalSubsteps)
    call advance_temperature_rk4(0.5d0*dt)
    call advance_temperature_rk4(0.5d0*dt)
    !$acc wait(1)
    call measure('rk4_two_half_steps_mode_amplification',maxval(abs(T-baseT-amp*amplification*mode)),5.0d-13)
    halfError=maxval(abs(T-baseT-amp*exp(lambda*dt)*mode))
    call require(fullError>1.0d-11,'rk4_order_test_above_roundoff')
    call measure('rk4_exponential_error_refinement_ratio',halfError/fullError,0.10d0)
    write(*,'(a)') 'ALL_MATHEMATICAL_CHECKS_PASSED'
contains
    subroutine moment_definition(distribution,sx,sy,result)
      real(kind=8),intent(in) :: distribution(0:8),sx,sy
      real(kind=8),intent(out) :: result(0:8)
      integer :: d,l
      result=0.0d0
      do l=0,8
        do d=0,8
          result(l)=result(l)+distribution(d)*(dble(cx(d))-sx)**px(l)*(dble(cy(d))-sy)**py(l)
        enddo
      enddo
    end subroutine
    subroutine measure(label,residual,limit)
      character(len=*),intent(in) :: label
      real(kind=8),intent(in) :: residual,limit
      write(*,'(a,1x,a,2(1x,es24.16e3))') 'CHECK',label,residual,limit
      if(.not.(residual<=limit)) error stop 41
    end subroutine
    subroutine require(condition,label)
      logical,intent(in) :: condition
      character(len=*),intent(in) :: label
      if(.not.condition) then
        write(*,*) 'FAILED ',label
        error stop 42
      endif
    end subroutine
end program main
"""


def mathematical_main() -> str:
    """从独立一维半格边界扩散矩阵构造模态，再交给实际 Fortran RHS/RK4。

    Dirichlet 边界的离散壁面通量采用三点二阶公式；其余内点是中心差分。
    用相似变换把该三对角矩阵转为对称矩阵，再选第三个衰减模态，以使
    RK4 的时间误差明显高于浮点舍入。没有复制实际热 RHS 的逐点实现。
    """
    result = MATH_MAIN
    for axis, cells in (("X", NX), ("Y", NY)):
        operator = np.diag(np.full(cells, -2.0))
        operator += np.diag(np.ones(cells-1), 1) + np.diag(np.ones(cells-1), -1)
        operator[0, 0] = operator[-1, -1] = -4.0
        operator[0, 1] = operator[-1, -2] = 4.0/3.0
        scaling = np.ones(cells)
        scaling[0] = scaling[-1] = math.sqrt(4.0/3.0)
        symmetric = operator * scaling[None, :] / scaling[:, None]
        np.testing.assert_allclose(symmetric, symmetric.T, rtol=0, atol=1e-15)
        values, vectors = np.linalg.eigh(symmetric)
        decay = values[-3]
        vector = scaling * vectors[:, -3]
        vector /= max(abs(vector))
        if vector[0] < 0:
            vector = -vector
        np.testing.assert_allclose(operator @ vector, decay * vector, rtol=0, atol=8e-15)
        array_text = ",".join(f"{float(value):.17e}".replace("e", "d") for value in vector)
        result = result.replace(f"__EIGEN_{axis}__", array_text)
        result = result.replace(f"__DECAY_{axis}__", f"{float(decay):.17e}".replace("e", "d"))
    return result


FLAGS = ["-cpp", "-ffree-line-length-none", "-fopenacc", "-O0", "-g",
         "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]


def command(args: list[str], directory: Path, env: dict[str, str],
            timeout: float) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(args, cwd=directory, env=env, capture_output=True, text=True,
                            encoding="utf-8", errors="replace", timeout=timeout)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(args)}\n"
                           + (result.stdout + result.stderr)[-12000:])
    return result


def build(source: str, directory: Path, compiler: str,
          env: dict[str, str], timeout: float, *, run: bool = True) -> str:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "solver.F90").write_text(source, encoding="utf-8")
    check = command([compiler, *FLAGS, "-fsyntax-only", "solver.F90"], directory, env, timeout)
    compiled = command([compiler, *FLAGS, "solver.F90", "-o", "solver.exe"], directory, env, timeout)
    (directory / "compile.log").write_text(check.stdout + check.stderr + compiled.stdout + compiled.stderr,
                                            encoding="utf-8")
    if not run:
        return ""
    result = command([str(directory / "solver.exe")], directory, env, timeout)
    output = result.stdout + result.stderr
    (directory / "run.log").write_text(output, encoding="utf-8")
    return output


def read_records(path: Path) -> list[bytes]:
    data = path.read_bytes()
    result: list[bytes] = []
    offset = 0
    while offset < len(data):
        size, = struct.unpack_from("<i", data, offset)
        if size < 0 or offset + size + 8 > len(data):
            raise AssertionError(f"Invalid sequential record: {path.name}")
        end, = struct.unpack_from("<i", data, offset + size + 4)
        if end != size:
            raise AssertionError(f"Record marker mismatch: {path.name}")
        result.append(data[offset + 4:offset + size + 4])
        offset += size + 8
    return result


def finite_records(records: list[bytes]) -> None:
    for record in records:
        values = struct.unpack(f"<{len(record)//8}d", record)
        if not all(math.isfinite(x) for x in values):
            raise AssertionError("Nonfinite checkpoint or snapshot")


def metadata(directory: Path) -> tuple[Path, dict[str, str]]:
    files = list(directory.glob("reload*-latest.meta"))
    if len(files) != 1:
        raise AssertionError("Expected exactly one latest restart metadata")
    entries = {}
    for line in files[0].read_text(encoding="utf-8").splitlines():
        pair = line.split(maxsplit=1)
        if len(pair) == 2:
            entries[pair[0]] = pair[1].strip()
    state = files[0].with_name(files[0].name.replace("latest.meta", entries["reloadFileName"] + ".bin"))
    return state, entries


def inspect_run(directory: Path, *, steady: bool) -> dict[str, object]:
    state, meta = metadata(directory)
    records = read_records(state)
    header = records.pop(0)
    if len(header) != 48 or struct.unpack("<4i", header[32:]) != (NX, NY, 9, 0):
        raise AssertionError("Checkpoint must identify one D2Q9 distribution and a scalar temperature")
    expected_sizes = [NX * NY * 9 * 8, NX * NY * 8]
    if steady:
        expected_sizes.append((3 * NX * NY + 2) * 8)
    if [len(r) for r in records] != expected_sizes:
        raise AssertionError("Incorrect f/T/steady-reference checkpoint record sizes")
    finite_records(records)
    mass = sum(struct.unpack(f"<{NX*NY*9}d", records[0]))
    mass_error = abs(mass - NX*NY) / (NX*NY)
    if mass_error > 1.0e-11:
        raise AssertionError(f"Closed-cavity mass changed by {mass_error}")
    snapshots = list(directory.glob("*Snapshot-*.bin"))
    if not snapshots:
        raise AssertionError("No post-processing snapshots")
    for snapshot in snapshots:
        fields = read_records(snapshot)
        if [len(r) for r in fields] != [NX*NY*8] * 4:
            raise AssertionError("Incorrect u/v/T/rho snapshot layout")
        finite_records(fields)
        if min(struct.unpack(f"<{NX*NY}d", fields[3])) <= 0:
            raise AssertionError("Nonpositive snapshot density")
    series = {}
    for prefix in ("Nu_VolAvg", "Re_VolAvg"):
        files = list(directory.glob(prefix + "*.dat"))
        if len(files) != 1:
            raise AssertionError(f"Missing {prefix} history")
        rows = [tuple(map(float, line.split())) for line in files[0].read_text().splitlines() if line.strip()]
        if not rows or any(len(row) != 2 or not all(map(math.isfinite, row)) for row in rows):
            raise AssertionError(f"Invalid or nonfinite {prefix}")
        if any(rows[i][0] >= rows[i+1][0] for i in range(len(rows)-1)):
            raise AssertionError(f"Repeated or unordered history times in {prefix}")
        series[prefix] = len(rows)
    saved_samples = int(meta["nureSampleCount"])
    if any(count != saved_samples for count in series.values()):
        raise AssertionError("Saved Nu/Re counter differs from actual history rows")
    return {"itc_total": int(meta["itc_total"]), "snapshot_count": len(snapshots),
            "checkpoint_records": [len(r) for r in records],
            "closed_cavity_relative_mass_error": mass_error, "finite_Nu_Re_samples": series,
            "checkpoint_saved_Nu_Re_samples": saved_samples}


def integration_checks(source: str, root: Path, compiler: str,
                       env: dict[str, str], timeout: float) -> dict[str, object]:
    results = {}
    for name, steady, rb in (("side_steady", True, False), ("rb_unsteady", False, True)):
        case = root / name
        build(small_source(source, steady=steady, rb=rb, steps=48), case, compiler, env, timeout)
        results[name] = inspect_run(case, steady=steady)
    for name, steady, rb in (("side", True, False), ("rb", False, True)):
        full = root / ("side_steady" if steady else "rb_unsteady")
        split = root / (name + "_restart")
        build(small_source(source, steady=steady, rb=rb, steps=32), split, compiler, env, timeout)
        before = inspect_run(split, steady=steady)
        histories = {p.name: p.read_bytes() for prefix in ("Nu_VolAvg", "Re_VolAvg")
                     for p in split.glob(prefix + "*.dat")}
        build(small_source(source, steady=steady, rb=rb, steps=48, restart=True), split, compiler, env, timeout)
        after = inspect_run(split, steady=steady)
        full_state, full_meta = metadata(full)
        split_state, split_meta = metadata(split)
        if int(full_meta["itc_total"]) != 48 or int(split_meta["itc_total"]) != 48:
            raise AssertionError("Absolute restart time horizon mismatch")
        if full_state.read_bytes() != split_state.read_bytes():
            raise AssertionError(f"{name}: f/T or steady reference differs after strict continuation")
        full_snapshots = {p.name: p.read_bytes() for p in full.glob("*Snapshot-*.bin")}
        split_snapshots = {p.name: p.read_bytes() for p in split.glob("*Snapshot-*.bin")}
        if (not steady and full_snapshots.keys() != split_snapshots.keys()) or any(
                split_snapshots.get(k) != v for k,v in full_snapshots.items()):
            raise AssertionError(f"{name}: strict-continuation snapshot mismatch")
        for filename, old in histories.items():
            current = (split / filename).read_bytes()
            if not current.startswith(old):
                raise AssertionError(f"{name}: previous history overwritten")
            if current.splitlines()[-1] != (full / filename).read_bytes().splitlines()[-1]:
                raise AssertionError(f"{name}: final Nu/Re differs from uninterrupted run")
        if steady and after["checkpoint_saved_Nu_Re_samples"] != before["checkpoint_saved_Nu_Re_samples"] + 1:
            raise AssertionError("Steady final Nu/Re counter not appended")
        protected = {p.name: p.read_bytes() for p in split.iterdir()
                     if p.name in histories or p.name.endswith("-latest.meta")}
        command([str(split / "solver.exe")], split, env, timeout)
        if any((split / k).read_bytes() != v for k,v in protected.items()):
            raise AssertionError("Completed-target restart changed history or metadata")
        if steady:
            # 放宽停止阈值且增加总步限：已保存的收敛状态仍应零步退出。
            converged = small_source(source, steady=True, rb=rb, steps=96, restart=True)
            converged = replace_parameter(converged, "epsU", "1.0d6")
            converged = replace_parameter(converged, "epsT", "1.0d6")
            build(converged, split, compiler, env, timeout)
            if any((split / k).read_bytes() != v for k,v in protected.items()):
                raise AssertionError("Converged-state restart changed history or metadata")
        results[name + "_strict_restart"] = {"before": before, "after": after,
                "f_T_and_steady_reference_bitwise_equal": True, "snapshot_bitwise_equal": True,
                "Nu_Re_history_preserved": True, "final_Nu_Re_row_bitwise_equal": True,
                "completed_target_restart_preserved": True}
        if steady:
            results[name + "_strict_restart"]["converged_restart_preserved"] = True
    return results


def nu_values(output: str) -> dict[str, float]:
    pattern = r"(?m)^\s*(Nu_global|Nu_hot|Nu_cold|Nu_middle)(?:\([^)]*\))?\s*=\s*([-+0-9.EeDd]+)"
    result = {label: float(value.replace("D", "E")) for label, value in re.findall(pattern, output)}
    if not {"Nu_global", "Nu_hot", "Nu_cold"}.issubset(result):
        raise AssertionError("Missing final wall or volume Nu diagnostics")
    return result


def physical_semantics_checks(source: str, root: Path, compiler: str,
                             env: dict[str, str], timeout: float) -> dict[str, object]:
    results = {}
    for label, rb in (("side", False), ("rb", True)):
        # 两个工况均采用稳态编译，以检查实际最终壁面/全场 Nu 后处理。
        reference = root / "side_steady" if not rb else root / "semantics_rb_reference"
        if rb:
            build(small_source(source, steady=True, rb=True, steps=48), reference, compiler, env, timeout)
        state, meta = metadata(reference)
        records = read_records(state)
        reference_flow = np.frombuffer(records[1], dtype="<f8")
        reference_temperature = np.frombuffer(records[2], dtype="<f8")
        reference_nu = nu_values((reference / "run.log").read_text(encoding="utf-8"))

        translated = root / ("semantics_" + label + "_temperature_shift")
        clone = small_source(source, steady=True, rb=rb, steps=48)
        physics = list(map(float, meta["physicalParameters"].replace("D", "E").split()))
        shift = 3.25
        clone = replace_parameter(clone, "Thot", f"{physics[3]+shift:.17e}".replace("e", "d"))
        clone = replace_parameter(clone, "Tcold", f"{physics[4]+shift:.17e}".replace("e", "d"))
        output = build(clone, translated, compiler, env, timeout)
        translated_state, _ = metadata(translated)
        translated_records = read_records(translated_state)
        flow_error = float(max(abs(np.frombuffer(translated_records[1], dtype="<f8") - reference_flow)))
        temperature_error = float(max(abs(np.frombuffer(translated_records[2], dtype="<f8") - reference_temperature - shift)))
        shifted_nu = nu_values(output)
        if shifted_nu.keys() != reference_nu.keys():
            raise AssertionError("Translated temperature changed diagnostic availability")
        nu_shift_error = max(abs(shifted_nu[k] - reference_nu[k]) for k in reference_nu)
        if flow_error > 2e-12 or temperature_error > 2e-12 or nu_shift_error > 1e-10:
            raise AssertionError(f"{label}: absolute temperature reference affected physical results")

        scaled = root / ("semantics_" + label + "_velocity_units")
        clone = small_source(source, steady=True, rb=rb, steps=48)
        clone = replace_parameter(clone, "velocityScaleCompare", "1.0d0")
        output = build(clone, scaled, compiler, env, timeout)
        scaled_state, _ = metadata(scaled)
        if scaled_state.read_bytes() != state.read_bytes():
            raise AssertionError("Changing display velocity units changed lattice state")
        scaled_nu = nu_values(output)
        if scaled_nu.keys() != reference_nu.keys():
            raise AssertionError("Velocity display scale changed diagnostic availability")
        nu_scale_error = max(abs(scaled_nu[k] - reference_nu[k]) for k in reference_nu)
        if nu_scale_error > 1e-12:
            raise AssertionError("Changing display velocity units changed physical Nu")
        results[label] = {"temperature_reference_shift": shift, "max_f_shift_error": flow_error,
                          "max_T_translation_error": temperature_error, "max_Nu_translation_error": nu_shift_error,
                          "velocity_display_scale_state_bitwise_equal": True,
                          "max_Nu_display_scale_error": nu_scale_error}
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "2DRBOpenaccHybridRLB4.F90")
    parser.add_argument("--compiler", default="C:/msys64/ucrt64/bin/gfortran.exe")
    parser.add_argument("--report", type=Path, default=HERE / "validation_local.json")
    parser.add_argument("--temp-root", type=Path, help="Existing writable ASCII directory")
    parser.add_argument("--timeout", type=float, default=90)
    parser.add_argument("--math-only", action="store_true")
    args = parser.parse_args()
    source_bytes = args.source.read_bytes()
    source = source_bytes.decode("utf-8").replace("\r\n", "\n")
    env = os.environ.copy()
    env["PATH"] = str(Path(args.compiler).parent) + os.pathsep + env.get("PATH", "")
    env["ACC_DEVICE_TYPE"] = "host"
    env["OMP_NUM_THREADS"] = "1"
    report = {"source": str(args.source.resolve()), "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
              "compiler": args.compiler, "backend": "gfortran OpenACC host; NVIDIA GPU not tested",
              "scope": "independent flow moments, scalar spatial/RK4 operators, finite-step smoke and restart",
              "high_Ra_accuracy_validated": False, "grid_independence_validated": False}
    start = time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix="hybrid-verify-", dir=args.temp_root) as temporary:
            if not temporary.isascii():
                raise RuntimeError("Pass --temp-root with a writable ASCII directory")
            root = Path(temporary)
            report["compiler_version"] = command([args.compiler, "--version"], root, env, args.timeout).stdout.splitlines()[0]
            report["compile_matrix"] = {}
            for steady in (True, False):
                for rb in (False, True):
                    label = ("steady" if steady else "unsteady") + ("_rb" if rb else "_side")
                    clone = small_source(source, steady=steady, rb=rb, steps=48)
                    build(clone, root / label, args.compiler, env, args.timeout, run=False)
                    report["compile_matrix"][label] = True
            report["mathematical_checks"] = {}
            for label, rb in (("side_heated", False), ("rayleigh_benard", True)):
                clone = small_source(source, steady=True, rb=rb, steps=48)
                output = build(replace_main(clone, mathematical_main()), root / ("math_" + label), args.compiler, env, args.timeout)
                metrics = {}
                for line in output.splitlines():
                    if line.startswith("CHECK "):
                        _, name, residual, limit = line.split()
                        residual, limit = float(residual.replace("D", "E")), float(limit.replace("D", "E"))
                        metrics[name] = {"max_abs_residual": max(residual, metrics.get(name, {}).get("max_abs_residual", 0.0)),
                                         "tolerance": limit}
                if "ALL_MATHEMATICAL_CHECKS_PASSED" not in output:
                    raise AssertionError("Mathematical harness did not finish")
                report["mathematical_checks"][label] = metrics
            if not args.math_only:
                report["integration_checks"] = integration_checks(source, root, args.compiler, env, args.timeout)
                report["physical_semantics_checks"] = physical_semantics_checks(source, root, args.compiler, env, args.timeout)
            report["passed"] = True
    except Exception as exc:
        report["passed"] = False
        report["failure"] = str(exc)
        raise
    finally:
        report["elapsed_seconds"] = round(time.monotonic() - start, 3)
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        summary = {"passed": report.get("passed", False), "report": str(args.report.resolve()),
                   "source_sha256": report["source_sha256"], "elapsed_seconds": report["elapsed_seconds"]}
        if "failure" in report:
            summary["failure"] = report["failure"]
        print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
