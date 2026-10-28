"""在临时目录核验实际二维 FPC 源码；不改参数原件，不访问服务器。

需要 Python 3 和 gfortran。默认使用本机 MSYS2 gfortran 的 OpenACC host
后端，以执行实际数据区和异步等待路径；这不验证 NVIDIA GPU 运行。
所有算例均为小网格、有限步的开发核验，不构成高 Ra 或网格独立性验证。
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


HERE = Path(__file__).resolve().parent


def replace_parameter(source: str, name: str, value: str) -> str:
    """只改临时源文件的单个标量参数；找不到时明确失败。"""
    pattern = rf"(?im)(\b{re.escape(name)}\s*=\s*)[^,!\n]+"
    changed, count = re.subn(pattern, lambda m: m[1] + value, source, count=1)
    if count != 1:
        raise AssertionError(f"Parameter {name} not found")
    return changed


def select_macros(source: str, *, steady: bool, rb: bool) -> str:
    wanted = {
        "steadyFlow": steady,
        "unsteadyFlow": not steady,
        "SideHeatedCell": not rb,
        "RayleighBenardCell": rb,
    }
    for macro, enabled in wanted.items():
        pattern = rf"(?m)^\s*!?\s*#define\s+{macro}\s*(?:![^\n]*)?$"
        replacement = ("#define " if enabled else "!#define ") + macro
        source, count = re.subn(pattern, replacement, source)
        if count == 0 and enabled:
            source = replacement + "\n" + source
    return source


def small_source(source: str, *, steady: bool, rb: bool, steps: int,
                 restart: bool = False) -> str:
    source = select_macros(source, steady=steady, rb=rb)
    source = replace_parameter(source, "nx", "12")
    source = replace_parameter(source, "ny", "10")
    source = replace_parameter(source, "Rayleigh", "1.0d4")
    source = replace_parameter(source, "loadInitField", "1" if restart else "0")
    # 输出时刻直接指定格子步，保留原有 t_ff <-> itc 换算逻辑。
    for key, value in (("outputSnapshotInterval", "7.0d0/timeUnit"),
                       ("reloadFileInterval", "16.0d0/timeUnit"),
                       ("outputPltFileInterval", "16.0d0/timeUnit")):
        # baseline 在稳态和非稳态条件块中分别定义同名参数。
        pattern = rf"(?im)^(\s*real\(kind=8\),\s*parameter\s*::\s*{key}\s*=\s*)[^!\n]+"
        source, count = re.subn(pattern, lambda m: m[1] + value + "   ", source)
        if not count:
            raise AssertionError(f"Output parameter {key} not found")
    source = replace_parameter(source, "unsteadyRunDuration", f"{steps}.0d0/timeUnit")
    # 稳态只检查有限步运行和诊断路径，收敛另在服务器物理验证中判断。
    if steady:
        source = replace_parameter(source, "itc_max", str(steps))
    return source


def replace_main(source: str, main: str) -> str:
    pattern = r"(?ims)^\s*program\s+main\s*$.*?^\s*end\s+program\s+main\s*$"
    source, count = re.subn(pattern, lambda _: "\n" + main + "\n", source, count=1)
    if count != 1:
        raise AssertionError("The actual source must contain one program main")
    return source


MATH_MAIN = r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,k,n
    integer, parameter :: px(0:8)=[0,1,0,2,0,1,2,1,2]
    integer, parameter :: py(0:8)=[0,0,1,0,2,1,1,2,2]
    integer, parameter :: cx(0:8)=[0,1,0,-1,0,1,-1,-1,1]
    integer, parameter :: cy(0:8)=[0,0,1,0,-1,1,1,-1,-1]
    real(kind=8) :: p(0:8),q(0:8),r(0:8),c(0:8),z(0:8),a(0:8),b(0:8)
    real(kind=8) :: ux,uy,rr,ff,err,scale,temp,expected,base(0:8),stress0(0:8)
    real(kind=8) :: initialMass, finalMass

    call acc_init(acc_device_host)
    call initial()
    call require(maxval(abs(ex-cx))+maxval(abs(ey-cy))==0,'direction_order')

    ! 独立 oracle：按离散速度和矩定义求和，不重写程序的展开式。
    do n=1,20
        do k=0,8
            p(k)=0.03d0+0.02d0*sin(dble((n+2)*(k+1)))
        enddo
        ux=0.17d0*cos(dble(n));uy=0.11d0*sin(dble(n))
        call moment_definition(p,0.0d0,0.0d0,a)
        call fpc_pop_to_raw(p,r)
        call measure('raw_definition',maxval(abs(r-a)),2.0d-15)
        call fpc_raw_to_pop(r,q)
        call measure('raw_roundtrip',maxval(abs(q-p)),2.0d-15)
        call fpc_raw_to_central(r,ux,uy,c)
        call moment_definition(p,ux,uy,a)
        call measure('central_definition',maxval(abs(c-a)),2.0d-15)
        call fpc_central_to_raw(c,ux,uy,z)
        call measure('central_roundtrip',maxval(abs(z-r)),2.0d-15)
    enddo

    ! 实际初值平衡构造器在非零速度及带符号温标下满足 Maxwellian 中心矩定义。
    ux=0.073d0;uy=-0.041d0
    do n=-1,1
        temp=0.3d0*dble(n)
        call fpc_equilibrium(temp,ux,uy,p)
        call moment_definition(p,ux,uy,a)
        b=0.0d0
        b(0)=temp;b(3)=cs2*temp;b(4)=cs2*temp;b(8)=cs2*cs2*temp
        call measure('equilibrium_central_moment_definition',maxval(abs(a-b)),3.0d-15)
    enddo

    ! 初值由实际 initial() 生成：包括 f 中的半步浮力修正，不能仅检查 sum(f)。
    err=0.0d0
    do j=1,ny
        do i=1,nx
            p=f(i,j,:);q=g(i,j,:)
            call moment_definition(p,0.0d0,0.0d0,a)
            call moment_definition(q,0.0d0,0.0d0,b)
            ff=rho(i,j)*gBeta*(T(i,j)-Tref)
            err=max(err,abs(a(0)-rho(i,j)),abs(b(0)-T(i,j)))
            err=max(err,abs(a(1)+0.5d0*Fx(i,j)-rho(i,j)*u(i,j)))
            err=max(err,abs(a(2)+0.5d0*ff-rho(i,j)*v(i,j)))
            err=max(err,abs(b(1)-T(i,j)*u(i,j)),abs(b(2)-T(i,j)*v(i,j)))
        enddo
    enddo
    call measure('initial_mass_temperature_physical_momentum',err,2.0d-14)

    ! 零力下任意非平衡输入的守恒及二阶应力松弛。
    rr=1.17d0
    base=rr*[4.0d0/9,1.0d0/9,1.0d0/9,1.0d0/9,1.0d0/9, &
             1.0d0/36,1.0d0/36,1.0d0/36,1.0d0/36]
    p=base
    p(0)=p(0)-0.004d0
    p(1)=p(1)+0.003d0;p(3)=p(3)+0.003d0
    p(2)=p(2)-0.001d0;p(4)=p(4)-0.001d0
    p(5)=p(5)+0.0007d0;p(7)=p(7)+0.0007d0
    p(6)=p(6)-0.0007d0;p(8)=p(8)-0.0007d0
    call moment_definition(p,0.0d0,0.0d0,stress0)
    do k=0,8
        f(:,:,k)=p(k)
    enddo
    rho=rr;u=0.0d0;v=0.0d0;T=Tref
    call collision()
    !$acc wait(1)
    q=f_post(1,1,:)
    call moment_definition(q,0.0d0,0.0d0,a)
    call measure('flow_zero_force_mass',abs(a(0)-stress0(0)),3.0d-14)
    call measure('flow_zero_force_momentum',maxval(abs(a(1:2)-stress0(1:2))),3.0d-14)
    expected=(1.0d0-fpcShearRate)*(stress0(3)-stress0(4))
    call measure('flow_normal_stress_rate',abs((a(3)-a(4))-expected),3.0d-14)
    expected=(1.0d0-fpcShearRate)*stress0(5)
    call measure('flow_shear_stress_rate',abs(a(5)-expected),3.0d-14)
    expected=stress0(3)+stress0(4)+fpcBulkRate*(2.0d0*cs2*rr-stress0(3)-stress0(4))
    call measure('flow_bulk_rate',abs(a(3)+a(4)-expected),3.0d-14)

    ! 非零力的独立检查：先按物理动量定义设速度，再检查分布一阶矩增量恰为 F。
    T=Tref+0.2d0
    ff=rr*gBeta*0.2d0
    u=stress0(1)/rr;v=(stress0(2)+0.5d0*ff)/rr
    call collision()
    !$acc wait(1)
    q=f_post(1,1,:)
    call moment_definition(q,0.0d0,0.0d0,a)
    call measure('flow_forced_mass',abs(a(0)-stress0(0)),3.0d-14)
    call measure('flow_forced_x_momentum',abs(a(1)-stress0(1)),3.0d-14)
    call measure('flow_forced_y_momentum_increment',abs(a(2)-stress0(2)-ff),3.0d-14)
    call require(abs(ff)>1.0d-9,'force_test_nontrivial')

    ! 含正温、零温、负温；独立检查唯一守恒量和两个一阶热流中心矩。
    ux=0.073d0;uy=-0.041d0
    u=ux;v=uy
    do n=-1,1
        temp=0.3d0*dble(n)
        p=base*(temp/rr)
        p(1)=p(1)+0.0012d0;p(3)=p(3)-0.0012d0
        p(2)=p(2)-0.0008d0;p(4)=p(4)+0.0008d0
        do k=0,8
            g(:,:,k)=p(k)
        enddo
        T=temp
        call moment_definition(p,ux,uy,b)
        call collisionT()
        !$acc wait(1)
        q=g_post(1,1,:)
        call moment_definition(q,ux,uy,a)
        call measure('thermal_temperature_conservation',abs(a(0)-b(0)),3.0d-14)
        expected=(1.0d0-fpcThermalFluxRate)*b(1)
        call measure('thermal_x_flux_rate',abs(a(1)-expected),3.0d-14)
        expected=(1.0d0-fpcThermalFluxRate)*b(2)
        call measure('thermal_y_flux_rate',abs(a(2)-expected),3.0d-14)
    enddo

    ! 纯导热核验：只推进实际温度碰撞/迁移/边界，u=v=0，避免人为改浮力定义。
    u=0.0d0;v=0.0d0
    do j=1,ny
        do i=1,nx
#ifdef SideHeatedCell
            T(i,j)=Thot+(dble(i)-0.5d0)/dble(nx)*(Tcold-Thot)
#else
            T(i,j)=Thot+(dble(j)-0.5d0)/dble(ny)*(Tcold-Thot)
#endif
            g(i,j,:)=omegaT(:)*T(i,j)
        enddo
    enddo
    initialMass=sum(T)
    do n=1,600
        call collisionT()
        call streamingT()
        call bouncebackT()
        call macroT()
    enddo
    !$acc wait(1)
    err=0.0d0
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
    call measure('linear_conduction_temperature',err,1.0d-10)
    finalMass=sum(T)
    call measure('linear_conduction_heat_balance',abs(finalMass-initialMass),1.0d-10)
#ifdef SideHeatedCell
    call SideHeatedcalc_Nu_global()
    call SideHeatedcalc_Nu_wall_avg()
#else
    call RBcalc_Nu_global()
    call RBcalc_Nu_wall_avg()
#endif
    call measure('linear_conduction_Nu_volume',abs(Nu_global-1.0d0),1.0d-9)
    call measure('linear_conduction_Nu_hot',abs(Nu_hot-1.0d0),1.0d-9)
    call measure('linear_conduction_Nu_cold',abs(Nu_cold-1.0d0),1.0d-9)
    write(*,'(a)') 'ALL_MATHEMATICAL_CHECKS_PASSED'
contains
    subroutine moment_definition(pop,sx,sy,mom)
        real(kind=8),intent(in) :: pop(0:8),sx,sy
        real(kind=8),intent(out) :: mom(0:8)
        integer :: d,l
        mom=0.0d0
        do l=0,8
            do d=0,8
                mom(l)=mom(l)+pop(d)*(dble(cx(d))-sx)**px(l)*(dble(cy(d))-sy)**py(l)
            enddo
        enddo
    end subroutine
    subroutine measure(label,residual,limit)
        character(len=*),intent(in) :: label
        real(kind=8),intent(in) :: residual,limit
        write(*,'(a,1x,a,2(1x,es24.16e3))') 'CHECK',label,residual,limit
        if(.not.(residual<=limit)) then
            write(*,*) 'FAILED ',label
            error stop 41
        endif
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


def run_command(command: list[str], directory: Path, env: dict[str, str],
                timeout: float) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, cwd=directory, env=env, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=timeout)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(command)}\n"
                           + (result.stdout + result.stderr)[-12000:])
    return result


def build_run(source: str, directory: Path, compiler: str,
              env: dict[str, str], timeout: float) -> str:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "solver.F90").write_text(source, encoding="utf-8")
    flags = ["-cpp", "-ffree-line-length-none", "-fopenacc", "-O0", "-g",
             "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
    syntax = run_command([compiler, *flags, "-fsyntax-only", "solver.F90"], directory, env, timeout)
    compiled = run_command([compiler, *flags, "solver.F90", "-o", "solver.exe"], directory, env, timeout)
    (directory / "compile.log").write_text(syntax.stdout + syntax.stderr + compiled.stdout + compiled.stderr,
                                            encoding="utf-8")
    result = run_command([str(directory / "solver.exe")], directory, env, timeout)
    output = result.stdout + result.stderr
    (directory / "run.log").write_text(output, encoding="utf-8")
    return output


def read_records(path: Path) -> list[bytes]:
    """读取 gfortran 默认四字节 record marker 的 sequential-unformatted 输出。"""
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
        values = struct.unpack(f"<{len(record) // 8}d", record)
        if not all(math.isfinite(v) for v in values):
            raise AssertionError("Nonfinite distribution or snapshot data")


def metadata(directory: Path) -> tuple[Path, dict[str, str]]:
    files = list(directory.glob("reload*-latest.meta"))
    if len(files) != 1:
        raise AssertionError(f"Expected one latest restart metadata, found {files}")
    entries: dict[str, str] = {}
    for line in files[0].read_text(encoding="utf-8").splitlines():
        pair = line.split(maxsplit=1)
        if len(pair) == 2:
            entries[pair[0]] = pair[1].strip()
    state = files[0].with_name(files[0].name.replace("latest.meta", entries["reloadFileName"] + ".bin"))
    return state, entries


def inspect_run(directory: Path) -> dict[str, object]:
    state, meta = metadata(directory)
    distribution_records = read_records(state)
    header = distribution_records.pop(0)
    if len(header) != 48 or not header[:32].rstrip(b" ").startswith(b"FPC2D_D2Q9_D2Q9_V1"):
        raise AssertionError("Checkpoint FPC identifier missing")
    if struct.unpack("<4i", header[32:]) != (12, 10, 9, 9):
        raise AssertionError("Checkpoint FPC lattice dimensions mismatch")
    if [len(r) for r in distribution_records] != [12 * 10 * 9 * 8] * 2:
        raise AssertionError("Checkpoint must contain two full D2Q9 distributions")
    finite_records(distribution_records)
    flow = struct.unpack(f"<{12 * 10 * 9}d", distribution_records[0])
    relative_mass_error = abs(sum(flow) - 12 * 10) / (12 * 10)
    if relative_mass_error > 1.0e-11:
        raise AssertionError(f"Closed cavity mass changed: relative error {relative_mass_error}")
    snapshots = sorted(directory.glob("*Snapshot-*.bin"))
    if not snapshots:
        raise AssertionError("No post-processing snapshot written")
    for snapshot in snapshots:
        records = read_records(snapshot)
        if [len(r) for r in records] != [12 * 10 * 8] * 4:
            raise AssertionError("Snapshot must contain u/v/T/rho as four full fields")
        finite_records(records)
        if min(struct.unpack(f"<{12 * 10}d", records[3])) <= 0:
            raise AssertionError("Nonpositive density in snapshot")
    series: dict[str, int] = {}
    for prefix in ("Nu_VolAvg", "Re_VolAvg"):
        files = list(directory.glob(prefix + "*.dat"))
        if len(files) != 1:
            raise AssertionError(f"Expected one {prefix} time series")
        rows = [tuple(map(float, line.split())) for line in files[0].read_text().splitlines()
                if line.strip() and not line.lstrip().startswith("#")]
        if not rows or any(len(row) != 2 or not all(math.isfinite(x) for x in row) for row in rows):
            raise AssertionError(f"Invalid or nonfinite {prefix} time series")
        if any(rows[i][0] >= rows[i+1][0] for i in range(len(rows)-1)):
            raise AssertionError(f"Duplicate or unordered {prefix} time stamps")
        series[prefix] = len(rows)
    saved_samples = int(meta["nureSampleCount"])
    if any(count != saved_samples for count in series.values()):
        raise AssertionError("Checkpoint Nu/Re count differs from preserved history rows")
    return {"itc_total": int(meta["itc_total"]), "snapshot_count": len(snapshots),
            "checkpoint_records": [len(r) for r in distribution_records],
            "closed_cavity_relative_mass_error": relative_mass_error,
            "finite_Nu_Re_samples": series,
            "checkpoint_saved_Nu_Re_samples": saved_samples}


def integration_checks(source: str, root: Path, compiler: str,
                       env: dict[str, str], timeout: float) -> dict[str, object]:
    results: dict[str, object] = {}
    for name, steady, rb in (("side_steady", True, False),
                            ("rb_unsteady", False, True)):
        directory = root / name
        clone = small_source(source, steady=steady, rb=rb, steps=48)
        build_run(clone, directory, compiler, env, timeout)
        results[name] = inspect_run(directory)

    # 不断电 48 步 vs 32 步 + 严格续算到 48 步；检查 f/g 的二进制逐位一致。
    for name, steady, rb in (("side", True, False), ("rb", False, True)):
        full = root / ("side_steady" if steady else "rb_unsteady")
        split = root / (name + "_restart")
        first = small_source(source, steady=steady, rb=rb, steps=32)
        build_run(first, split, compiler, env, timeout)
        before = inspect_run(split)
        histories_before = {p.name: p.read_bytes() for prefix in ("Nu_VolAvg", "Re_VolAvg")
                            for p in split.glob(prefix + "*.dat")}
        resumed = small_source(source, steady=steady, rb=rb, steps=48, restart=True)
        build_run(resumed, split, compiler, env, timeout)
        after = inspect_run(split)
        if steady and after["checkpoint_saved_Nu_Re_samples"] != before["checkpoint_saved_Nu_Re_samples"] + 1:
            raise AssertionError("Steady restart must append and save one additional final Nu/Re row")
        full_state, full_meta = metadata(full)
        split_state, split_meta = metadata(split)
        if int(full_meta["itc_total"]) != 48 or int(split_meta["itc_total"]) != 48:
            raise AssertionError(f"{name}: restart or absolute step horizon mismatch")
        if full_state.read_bytes() != split_state.read_bytes():
            a = read_records(full_state)[1:]
            b = read_records(split_state)[1:]
            maxdiff = max(max(abs(x-y) for x,y in zip(struct.unpack(f'<{len(xr)//8}d',xr),
                                                     struct.unpack(f'<{len(yr)//8}d',yr)))
                          for xr,yr in zip(a,b))
            raise AssertionError(f"{name}: strict restart differs, max distribution error={maxdiff}")
        # 同一物理时刻的无量纲后处理快照也应逐位一致。
        full_snapshots = {p.name:p.read_bytes() for p in full.glob("*Snapshot-*.bin")}
        split_snapshots = {p.name:p.read_bytes() for p in split.glob("*Snapshot-*.bin")}
        # 稳态第一段结束时会额外写出 step=32 的停止快照，这个额外文件可保留。
        if (not steady and full_snapshots.keys() != split_snapshots.keys()) or any(
            split_snapshots.get(name) != content for name, content in full_snapshots.items()
        ):
            raise AssertionError(f"{name}: restart snapshot timeline or fields differ")
        for filename, old_history in histories_before.items():
            current = (split / filename).read_bytes()
            if not current.startswith(old_history):
                raise AssertionError(f"{name}: restart overwrote old history {filename}")
            if current.splitlines()[-1] != (full / filename).read_bytes().splitlines()[-1]:
                raise AssertionError(f"{name}: final Nu/Re history row differs {filename}")
        # 同一总目标已完成时再次启动，应直接结束，保留历史与检查点账本。
        protected = {p.name: p.read_bytes() for p in split.iterdir()
                     if p.name in histories_before or p.name.endswith("-latest.meta")}
        run_command([str(split / "solver.exe")], split, env, timeout)
        if any((split / filename).read_bytes() != content for filename, content in protected.items()):
            raise AssertionError(f"{name}: restarting an already completed target modified history or metadata")
        results[name + "_strict_restart"] = {"before": before, "after": after,
                                             "fg_bitwise_equal": True,
                                             "snapshot_bitwise_equal": True,
                                             "Nu_Re_history_preserved": True,
                                             "final_Nu_Re_row_bitwise_equal": True,
                                             "completed_target_restart_preserved": True}
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "2DRBOpenaccFPC.F90")
    parser.add_argument("--compiler", default="C:/msys64/ucrt64/bin/gfortran.exe")
    parser.add_argument("--report", type=Path, default=HERE / "validation_local.json")
    parser.add_argument("--timeout", type=float, default=90)
    parser.add_argument("--temp-root", type=Path,
                        help="An existing writable ASCII directory for temporary builds")
    parser.add_argument("--math-only", action="store_true")
    args = parser.parse_args()
    source_bytes = args.source.read_bytes()
    source = source_bytes.decode("utf-8").replace("\r\n", "\n")
    env = os.environ.copy()
    env["PATH"] = str(Path(args.compiler).parent) + os.pathsep + env.get("PATH", "")
    env["ACC_DEVICE_TYPE"] = "host"
    env["OMP_NUM_THREADS"] = "1"
    report: dict[str, object] = {
        "source": str(args.source.resolve()),
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "compiler": args.compiler,
        "backend": "gfortran OpenACC host; NVIDIA GPU not tested",
        "scope": "independent moment identities, pure conduction, finite-step small-grid smoke/restart",
        "high_Ra_accuracy_validated": False,
        "grid_independence_validated": False,
    }
    start = time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix="fpc-local-", dir=args.temp_root) as temporary:
            if not temporary.isascii():
                raise RuntimeError("An ASCII temporary directory is required; set TMP/TEMP to an ASCII path")
            root = Path(temporary)
            report["compiler_version"] = run_command([args.compiler, "--version"], root, env, args.timeout).stdout.splitlines()[0]
            report["mathematical_checks"] = {}
            for case, rb in (("side_heated", False), ("rayleigh_benard", True)):
                clone = small_source(source, steady=True, rb=rb, steps=48)
                output = build_run(replace_main(clone, MATH_MAIN), root / ("mathematics_" + case), args.compiler, env, args.timeout)
                metrics: dict[str, dict[str, float]] = {}
                for line in output.splitlines():
                    if line.startswith("CHECK "):
                        _, name, residual, limit = line.split()
                        value = float(residual.replace("D", "E"))
                        metrics[name] = {"max_abs_residual": max(value, metrics.get(name, {}).get("max_abs_residual", 0.0)),
                                         "tolerance": float(limit.replace("D", "E"))}
                if "ALL_MATHEMATICAL_CHECKS_PASSED" not in output:
                    raise AssertionError("Mathematical harness did not finish")
                report["mathematical_checks"][case] = metrics
            if not args.math_only:
                report["integration_checks"] = integration_checks(source, root, args.compiler, env, args.timeout)
            report["passed"] = True
    except Exception as exc:
        report["passed"] = False
        report["failure"] = str(exc)
        raise
    finally:
        report["elapsed_seconds"] = round(time.monotonic() - start, 3)
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        summary = {"passed": report.get("passed", False),
                   "report": str(args.report.resolve()),
                   "source_sha256": report["source_sha256"],
                   "elapsed_seconds": report["elapsed_seconds"]}
        if "failure" in report:
            summary["failure"] = report["failure"]
        print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
