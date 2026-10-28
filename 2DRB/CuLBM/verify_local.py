"""实际 CuLBM 源码的独立局部、边界和严格续算核验。

需要 Python 3、NumPy 与 gfortran，使用 OpenACC host 后端和 ASCII 临时目录。
不修改主源码，不访问服务器，不宣称 NVIDIA GPU 或高 Ra 物理结果已验证。
热累积量载体严格要求 T>0；默认热流率下的温标背景依赖单独量化。
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, re, struct, subprocess, tempfile, time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
NX,NY=12,10
FLAGS=["-cpp","-ffree-line-length-none","-fopenacc","-O0","-g","-fcheck=all","-fbacktrace","-ffpe-trap=invalid,zero,overflow"]

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
    real(kind=8),parameter :: densities(3)=[0.83d0,1.17d0,2.73d0]
    real(kind=8),parameter :: temperatures(3)=[0.37d0,2.2d0,5.1d0]
    real(kind=8) :: pop(0:8),post(0:8),scaled(0:8),eq(0:8),raw(0:8),central(0:8)
    real(kind=8) :: a(0:8),b(0:8),before(0:8),after(0:8),target(0:8),roundtrip(0:8)
    real(kind=8) :: thermal(0:4),thermalPost(0:4),thermalScaled(0:4)
    real(kind=8) :: rr,ux,uy,forceX,forceY,temp,err,expected,meanX,meanY,varianceBefore,varianceAfter
    logical :: valid
    call acc_init(acc_device_host)
    call initial()
    call require(all(ex==cx).and.all(ey==cy),'direction_order')
    write(*,'(a,7(1x,es24.16e3))') 'CONFIG',Mach,lengthUnit,Rayleigh,thermalCs2,Prandtl,diffusivity,cuThermalFluxRate

    ! 线性变换独立 oracle：用方向矩的定义直接求和，包括非零均值。
    do n=1,16
      do k=0,8
        pop(k)=0.03d0+0.02d0*sin(dble((n+1)*(k+2)))
      enddo
      ux=0.17d0*cos(dble(n));uy=0.11d0*sin(dble(n))
      call moments(pop,0.0d0,0.0d0,a)
      call cu_pop_to_raw(pop,raw)
      call measure('raw_moment_definition',maxval(abs(raw-a)),3.0d-15)
      call cu_raw_to_pop(raw,roundtrip)
      call measure('population_raw_roundtrip',maxval(abs(roundtrip-pop)),3.0d-15)
      call moments(pop,ux,uy,a)
      call cu_raw_to_central(raw,ux,uy,central)
      call measure('central_moment_definition',maxval(abs(central-a)),3.0d-15)
      call cu_central_to_raw(central,ux,uy,b)
      call measure('central_raw_roundtrip',maxval(abs(b-raw)),3.0d-15)
    enddo

    ! 初始 f 必须满足论文的半步物理动量；g 的守恒质量是正温度。
    err=0.0d0
    do j=1,ny
      do i=1,nx
        call moments(f(i,j,:),0.0d0,0.0d0,a)
        forceY=rho(i,j)*gBeta*(T(i,j)-Tref)
        err=max(err,abs(a(0)-rho(i,j)),abs(a(1)+0.5d0*Fx(i,j)-rho(i,j)*u(i,j)))
        err=max(err,abs(a(2)+0.5d0*forceY-rho(i,j)*v(i,j)),abs(sum(g(i,j,:))-T(i,j)))
      enddo
    enddo
    call measure('initial_mass_temperature_half_force_momentum',err,4.0d-14)

    do n=1,3
      rr=densities(n);ux=0.073d0;uy=-0.041d0
      call cu_flow_equilibrium(rr,ux,uy,eq)
      call cumulants(eq,rr,a)
      target=0.0d0;target(3)=cs2;target(4)=cs2
      call measure('normalized_equilibrium_cumulants',maxval(abs(a(1:8)-target(1:8))),5.0d-15)
      pop=eq
      pop(0)=pop(0)-rr*0.004d0
      pop(1)=pop(1)+rr*0.003d0;pop(3)=pop(3)+rr*0.003d0
      pop(2)=pop(2)-rr*0.001d0;pop(4)=pop(4)-rr*0.001d0
      pop(5)=pop(5)+rr*0.0007d0;pop(7)=pop(7)+rr*0.0007d0
      pop(6)=pop(6)-rr*0.0007d0;pop(8)=pop(8)-rr*0.0007d0
      do k=0,8
        pop(k)=pop(k)+rr*0.004d0*omega(k)*(dble(cx(k)*cx(k))-cs2)*dble(cy(k)) &
                       +rr*0.006d0*omega(k)*(dble(cx(k)*cx(k))-cs2)*(dble(cy(k)*cy(k))-cs2)
      enddo
      call moments(pop,0.0d0,0.0d0,before)
      call cumulants(pop,rr,a)
      call cu_flow_collide_node(pop,rr,ux,uy,0.0d0,0.0d0,post)
      call moments(post,0.0d0,0.0d0,after)
      call measure('flow_zero_force_mass_momentum',maxval(abs(after(0:2)-before(0:2))),5.0d-14)
      call cumulants(post,rr,b)
      do k=3,5
        expected=(1.0d0-cuViscousRate)*a(k)+cuViscousRate*target(k)
        call measure('normalized_viscous_cumulant_rate',abs(b(k)-expected),5.0d-14)
      enddo
      call measure('normalized_third_cumulant_rate',maxval(abs(b(6:7)-(1.0d0-cuThirdRate)*a(6:7))),5.0d-14)
      call measure('normalized_fourth_cumulant_inverse_and_rate',abs(b(8)-(1.0d0-cuFourthRate)*a(8)),5.0d-14)
      call require(abs(a(8))>1.0d-6,'fourth_cumulant_test_nontrivial')

      forceX=1.3d-4;forceY=-2.1d-4
      meanX=before(1)/rr;meanY=before(2)/rr
      ux=meanX+0.5d0*forceX/rr;uy=meanY+0.5d0*forceY/rr
      call cu_flow_collide_node(pop,rr,ux,uy,forceX,forceY,post)
      call moments(post,0.0d0,0.0d0,after)
      call measure('forced_mass',abs(after(0)-before(0)),5.0d-14)
      call measure('forced_raw_momentum_increment',maxval(abs(after(1:2)-before(1:2)-[forceX,forceY])),5.0d-14)
      call cu_flow_collide_node(2.5d0*pop,2.5d0*rr,ux,uy,2.5d0*forceX,2.5d0*forceY,scaled)
      call measure('density_force_homogeneity',maxval(abs(scaled-2.5d0*post)),8.0d-14)
    enddo

    ! proper mass-scaled thermal cumulants 的独立方差定义，不复制碰撞展开式。
    ux=0.071d0;uy=-0.037d0
    do n=1,3
      temp=temperatures(n)
      call cu_thermal_equilibrium(temp,ux,uy,thermal)
      call moments(thermal,ux,uy,a)
      call measure('thermal_equilibrium_mass_flux',maxval(abs(a(0:2)-[temp,0.0d0,0.0d0])),5.0d-15)
      call measure('thermal_equilibrium_variances',maxval(abs(a(3:4)/temp-thermalCs2)),5.0d-15)
      thermal(1)=thermal(1)+0.002d0*temp;thermal(3)=thermal(3)-0.002d0*temp
      thermal(2)=thermal(2)-0.0013d0*temp;thermal(4)=thermal(4)+0.0013d0*temp
      call moments(thermal,ux,uy,before)
      call cu_thermal_collide_node(thermal,temp,ux,uy,thermalPost,valid)
      call require(valid,'positive_temperature_accepted')
      call moments(thermalPost,ux,uy,after)
      call measure('thermal_mass_conservation',abs(after(0)-temp),5.0d-14)
      call measure('thermal_flux_rate',maxval(abs(after(1:2)-(1.0d0-cuThermalFluxRate)*before(1:2))),5.0d-14)
      do k=3,4
        varianceBefore=before(k)/temp-(before(k-2)/temp)**2
        varianceAfter=after(k)/temp-(after(k-2)/temp)**2
        expected=varianceBefore+cuThermalSecondRate*(thermalCs2-varianceBefore)
        call measure('thermal_normalized_variance_rate',abs(varianceAfter-expected),5.0d-14)
      enddo
      call cu_thermal_collide_node(2.5d0*thermal,2.5d0*temp,ux,uy,thermalScaled,valid)
      call require(valid,'scaled_positive_temperature_accepted')
      call measure('thermal_mass_homogeneity',maxval(abs(thermalScaled-2.5d0*thermalPost)),8.0d-14)
    enddo

    ! 生产约定拒绝 T<=0：包含全零 g、零和但非零热流、整体负温度三种不同输入。
    do n=1,3
      thermal=0.0d0;temp=0.0d0
      if(n.EQ.2) then
        thermal(1)=0.001d0;thermal(3)=-0.001d0
      elseif(n.EQ.3) then
        temp=-1.0d0
        call cu_thermal_equilibrium(temp,0.0d0,0.0d0,thermal)
      endif
      call cu_thermal_collide_node(thermal,temp,0.0d0,0.0d0,thermalPost,valid)
      call require(.not.valid,'nonpositive_thermal_mass_rejected')
      call measure('invalid_input_returned_without_clipping',maxval(abs(thermalPost-thermal)),0.0d0)
    enddo
    write(*,'(a)') 'ALL_MATHEMATICAL_CHECKS_PASSED'
contains
    subroutine moments(distribution,sx,sy,result)
      real(kind=8),intent(in) :: distribution(0:),sx,sy
      real(kind=8),intent(out) :: result(0:8)
      integer :: d,l
      result=0.0d0
      do l=0,8
        do d=0,size(distribution)-1
          result(l)=result(l)+distribution(d)*(dble(cx(d))-sx)**px(l)*(dble(cy(d))-sy)**py(l)
        enddo
      enddo
    end subroutine
    subroutine cumulants(distribution,density,result)
      real(kind=8),intent(in) :: distribution(0:8),density
      real(kind=8),intent(out) :: result(0:8)
      real(kind=8) :: m(0:8),sx,sy
      call moments(distribution,0.0d0,0.0d0,m)
      sx=m(1)/density;sy=m(2)/density
      call moments(distribution,sx,sy,m)
      result=m/density
      result(8)=result(8)-result(3)*result(4)-2.0d0*result(5)**2
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

CONDUCTION_MAIN = r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,n
    real(kind=8) :: initialT(nx,ny),local(0:4),err
    logical,parameter :: exactCase=__EXACT_CASE__
    call acc_init(acc_device_host)
    call initial()
    u=0.0d0;v=0.0d0
    if(exactCase.and.abs(cuThermalFluxRate-1.0d0)>2.0d-12) error stop 'Conduction case requires heat flux rate 1'
    do j=1,ny
      do i=1,nx
#ifdef SideHeatedCell
        T(i,j)=Thot+(dble(i)-0.5d0)/dble(nx)*(Tcold-Thot)
#else
        T(i,j)=Thot+(dble(j)-0.5d0)/dble(ny)*(Tcold-Thot)
#endif
        call cu_thermal_equilibrium(T(i,j),0.0d0,0.0d0,local)
        g(i,j,:)=local
      enddo
    enddo
    initialT=T
    do n=1,600
      call collisionT()
      call streamingT()
      call bouncebackT()
      call macroT()
    enddo
    call check_fields_culbm()
    err=maxval(abs(T-initialT))
#ifdef SideHeatedCell
    call SideHeatedcalc_Nu_global()
    call SideHeatedcalc_Nu_wall_avg()
#else
    call RBcalc_Nu_global()
    call RBcalc_Nu_wall_avg()
#endif
    write(*,'(a,4(1x,es24.16e3))') 'CONDUCTION',err,Nu_global,Nu_hot,Nu_cold
    if(exactCase) then
      if(err>2.0d-11) error stop 'Linear conduction exact-case temperature mismatch'
      if(max(abs(Nu_global-1.0d0),abs(Nu_hot-1.0d0),abs(Nu_cold-1.0d0))>2.0d-9) &
          error stop 'Linear conduction exact-case Nu mismatch'
    endif
    write(*,'(a)') 'CONDUCTION_CHECK_FINISHED'
end program main
"""

INVALID_MAIN = r"""
program main
    use commondata
    use openacc
    implicit none
    call acc_init(acc_device_host)
    call initial()
    u=0.0d0;v=0.0d0
    T=__INVALID_T__
    g=0.0d0
    __INVALID_G__
    call collisionT()
    !$acc wait(1)
    if(cuCollisionStatus.NE.1) error stop 'CuLBM invalid flag was not latched'
    call check_fields_culbm()
    error stop 'CuLBM invalid state incorrectly accepted'
end program main
"""

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
    if (len(header) != 48 or header[:32].rstrip(b" ") != b"CULBM2D_D2Q9_D2Q5_V1"
            or struct.unpack("<4i", header[32:]) != (NX, NY, 9, 5)):
        raise AssertionError("Checkpoint must identify D2Q9 flow and D2Q5 thermal distributions")
    expected_sizes = [NX * NY * 9 * 8, NX * NY * 5 * 8]
    if steady:
        expected_sizes.append((3 * NX * NY + 2) * 8)
    if [len(r) for r in records] != expected_sizes:
        raise AssertionError("Incorrect f/g/steady-reference checkpoint record sizes")
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
            raise AssertionError(f"{name}: f/g or steady reference differs after strict continuation")
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
                "f_g_and_steady_reference_bitwise_equal": True, "snapshot_bitwise_equal": True,
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


def parse_metrics(output: str) -> dict[str, dict[str, float]]:
    metrics = {}
    for line in output.splitlines():
        if line.startswith("CHECK "):
            _, name, residual, limit = line.split()
            value = float(residual.replace("D", "E"))
            metrics[name] = {"max_abs_residual": max(value, metrics.get(name, {}).get("max_abs_residual", 0.0)),
                             "tolerance": float(limit.replace("D", "E"))}
    return metrics


def conduction_checks(source: str, root: Path, compiler: str, env: dict[str, str], timeout: float,
                      configurations: dict[str, list[float]]) -> dict[str, object]:
    results = {}
    for label, rb in (("side_heated", False), ("rayleigh_benard", True)):
        mach, length, rayleigh, ct2, pr, kappa, rate = configurations[label]
        # 只为独立核验边界选择 kappa=ct2/2，即热流率=1；由 Pr 调整保持 nu/kappa=Pr 一致。
        exact_kappa = ct2/2
        exact_pr = (mach*length/(exact_kappa*math.sqrt(3*rayleigh)))**2
        results[label] = {}
        for exact in (True, False):
            clone = small_source(source, steady=True, rb=rb, steps=48)
            if exact:
                clone = replace_parameter(clone, "Prandtl", f"{exact_pr:.17e}".replace("e", "d"))
            main = CONDUCTION_MAIN.replace("__EXACT_CASE__", ".true." if exact else ".false.")
            output = build(replace_main(clone, main), root / ("conduction_" + label + ("_exact" if exact else "_default")),
                           compiler, env, timeout)
            row = next(line.split()[1:] for line in output.splitlines() if line.startswith("CONDUCTION "))
            error, nu_volume, nu_hot, nu_cold = map(float, row)
            results[label]["exact_boundary_case" if exact else "default_rate_observation"] = {
                "Prandtl": exact_pr if exact else pr, "diffusivity": exact_kappa if exact else kappa,
                "thermal_flux_rate": 1.0 if exact else rate, "steps": 600,
                "max_T_linear_deviation": error, "Nu_volume": nu_volume, "Nu_hot": nu_hot, "Nu_cold": nu_cold,
                "exact_Nu_1_required": exact,
                "interpretation": "boundary identity at heat-flux rate 1" if exact else
                    "finite-grid thermal-cumulant correction; observation only, not a benchmark"}
    return results


def invalid_state_checks(source: str, root: Path, compiler: str,
                         env: dict[str, str], timeout: float) -> dict[str, object]:
    results = {}
    for name, temp, assignments in (("zero_mass_zero_population", "0.0d0", ""),
                                    ("zero_mass_nonzero_heat_flux", "0.0d0", "g(:,:,1)=0.001d0;g(:,:,3)=-0.001d0"),
                                    ("negative_thermal_mass", "-1.0d0", "g(:,:,0)=-1.0d0")):
        clone = small_source(source, steady=True, rb=False, steps=48)
        main = INVALID_MAIN.replace("__INVALID_T__", temp).replace("__INVALID_G__", assignments)
        directory = root / ("invalid_" + name)
        build(replace_main(clone, main), directory, compiler, env, timeout, run=False)
        result = subprocess.run([str(directory / "solver.exe")], cwd=directory, env=env, capture_output=True,
                                text=True, encoding="utf-8", errors="replace", timeout=timeout)
        output = result.stdout + result.stderr
        if result.returncode == 0 or "CuLBM invalid collision state" not in output or "SIGFPE" in output:
            raise AssertionError(f"{name}: nonpositive thermal mass not rejected by the persistent production guard\n{output}")
        results[name] = {"rejected": True, "persistent_status_guard": True, "division_by_zero_avoided": True}
    return results


def scalar_semantics_checks(source: str, root: Path, compiler: str,
                           env: dict[str, str], timeout: float) -> dict[str, object]:
    results = {}
    for label, rb in (("side", False), ("rb", True)):
        reference = root / "side_steady" if not rb else root / "semantics_rb_reference"
        if rb:
            build(small_source(source, steady=True, rb=True, steps=48), reference, compiler, env, timeout)
        state, meta = metadata(reference)
        records = read_records(state)
        baseline_f = np.frombuffer(records[1], dtype="<f8")
        baseline_g = np.frombuffer(records[2], dtype="<f8")
        baseline_nu = nu_values((reference / "run.log").read_text(encoding="utf-8"))
        physics = list(map(float, meta["physicalParameters"].replace("D", "E").split()))
        case_result = {}
        for mode in ("amplitude", "background", "velocity_units"):
            clone = small_source(source, steady=True, rb=rb, steps=48)
            if mode == "amplitude":
                clone = replace_parameter(clone, "Thot", f"{2*physics[3]:.17e}".replace("e", "d"))
                clone = replace_parameter(clone, "Tcold", f"{2*physics[4]:.17e}".replace("e", "d"))
            elif mode == "background":
                clone = replace_parameter(clone, "Thot", f"{physics[3]+3.25:.17e}".replace("e", "d"))
                clone = replace_parameter(clone, "Tcold", f"{physics[4]+3.25:.17e}".replace("e", "d"))
            else:
                clone = replace_parameter(clone, "velocityScaleCompare", "1.0d0")
            directory = root / ("semantics_" + label + "_" + mode)
            output = build(clone, directory, compiler, env, timeout)
            current, _ = metadata(directory)
            current_records = read_records(current)
            f_values = np.frombuffer(current_records[1], dtype="<f8")
            g_values = np.frombuffer(current_records[2], dtype="<f8")
            nu = nu_values(output)
            f_error = float(max(abs(f_values-baseline_f)))
            nu_error = max(abs(nu[k]-baseline_nu[k]) for k in baseline_nu)
            if mode == "amplitude":
                g_error = float(max(abs(g_values-2*baseline_g)))
                if f_error > 2e-12 or g_error > 2e-12 or nu_error > 1e-10:
                    raise AssertionError("Positive thermal mass scaling changed dimensionless results")
                case_result[mode] = {"factor": 2.0, "max_f_error": f_error, "max_scaled_g_error": g_error,
                                     "max_Nu_error": nu_error, "required_invariance": True}
            elif mode == "background":
                # /T 热流平方项使有限网格上的加法平移不严格不变；保留数据而不伪装为通过条件。
                cells = NX*NY
                shift_T_error = float(max(abs(g_values.reshape(5, cells).sum(axis=0)
                                             -baseline_g.reshape(5, cells).sum(axis=0)-3.25)))
                case_result[mode] = {"shift": 3.25, "max_f_change": f_error,
                                     "max_T_translation_deviation": shift_T_error, "max_Nu_change": nu_error,
                                     "required_invariance": False,
                                     "interpretation": "finite-grid background dependence of proper thermal cumulants"}
            else:
                if current.read_bytes() != state.read_bytes() or nu_error > 1e-12:
                    raise AssertionError("Display velocity scale changed CuLBM state or physical Nu")
                case_result[mode] = {"state_bitwise_equal": True, "max_Nu_change": nu_error}
        results[label] = case_result
    return results


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,default=HERE/"2DRBOpenaccCuLBM.F90")
    parser.add_argument("--compiler",default="C:/msys64/ucrt64/bin/gfortran.exe")
    parser.add_argument("--report",type=Path,default=HERE/"validation_local.json")
    parser.add_argument("--temp-root",type=Path,help="Existing writable ASCII directory")
    parser.add_argument("--timeout",type=float,default=90)
    parser.add_argument("--math-only",action="store_true")
    args=parser.parse_args()
    source_bytes=args.source.read_bytes()
    source=source_bytes.decode("utf-8").replace("\r\n","\n")
    env=os.environ.copy()
    env["PATH"]=str(Path(args.compiler).parent)+os.pathsep+env.get("PATH","")
    env["ACC_DEVICE_TYPE"]="host";env["OMP_NUM_THREADS"]="1"
    report={"source":str(args.source.resolve()),"source_sha256":hashlib.sha256(source_bytes).hexdigest(),
            "compiler":args.compiler,"backend":"gfortran OpenACC host; NVIDIA GPU not tested",
            "scope":"normalized cumulants and inverse, positive thermal carrier, boundary identities, finite-step strict restart",
            "high_Ra_accuracy_validated":False,"grid_independence_validated":False,
            "thermal_carrier_domain":"strictly positive; zero and negative mass rejected",
            "temperature_origin_invariance_claimed":False}
    start=time.monotonic()
    try:
      with tempfile.TemporaryDirectory(prefix="culbm-verify-",dir=args.temp_root) as temporary:
        if not temporary.isascii():
            raise RuntimeError("Pass --temp-root with a writable ASCII directory")
        root=Path(temporary)
        report["compiler_version"]=command([args.compiler,"--version"],root,env,args.timeout).stdout.splitlines()[0]
        report["compile_matrix"]={}
        for steady in (True,False):
          for rb in (False,True):
            label=("steady" if steady else "unsteady")+("_rb" if rb else "_side")
            build(small_source(source,steady=steady,rb=rb,steps=48),root/label,args.compiler,env,args.timeout,run=False)
            report["compile_matrix"][label]=True
        report["mathematical_checks"]={}
        configurations={}
        for label,rb in (("side_heated",False),("rayleigh_benard",True)):
            clone=small_source(source,steady=True,rb=rb,steps=48)
            output=build(replace_main(clone,MATH_MAIN),root/("math_"+label),args.compiler,env,args.timeout)
            if "ALL_MATHEMATICAL_CHECKS_PASSED" not in output:
                raise AssertionError("Mathematical harness did not finish")
            report["mathematical_checks"][label]=parse_metrics(output)
            configurations[label]=list(map(float,next(line.split()[1:] for line in output.splitlines() if line.startswith("CONFIG "))))
        report["conduction_checks"]=conduction_checks(source,root,args.compiler,env,args.timeout,configurations)
        report["invalid_state_checks"]=invalid_state_checks(source,root,args.compiler,env,args.timeout)
        if not args.math_only:
            report["integration_checks"]=integration_checks(source,root,args.compiler,env,args.timeout)
            report["scalar_semantics_checks"]=scalar_semantics_checks(source,root,args.compiler,env,args.timeout)
        report["passed"]=True
    except Exception as exc:
        report["passed"]=False;report["failure"]=str(exc)
        raise
    finally:
        report["elapsed_seconds"]=round(time.monotonic()-start,3)
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        summary={"passed":report.get("passed",False),"report":str(args.report.resolve()),
                 "source_sha256":report["source_sha256"],"elapsed_seconds":report["elapsed_seconds"]}
        if "failure" in report:summary["failure"]=report["failure"]
        print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
