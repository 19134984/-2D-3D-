"""核验实际 Luo-MRT/Feng-ADE Fortran：局部矩、正支、导热、输出与严格续算。

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
ENTROPY_HISTORY = "EntropyHistory2DOpenaccEntropicADE.txt"


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
    integer :: i,j,k,n,status,fallbackCases
    integer, parameter :: cx(0:8)=[0,1,0,-1,0,1,-1,-1,1]
    integer, parameter :: cy(0:8)=[0,0,1,0,-1,1,1,-1,-1]
    real(kind=8), parameter :: wi(0:8)=[4d0/9,1d0/9,1d0/9,1d0/9,1d0/9, &
                                       1d0/36,1d0/36,1d0/36,1d0/36]
    real(kind=8) :: p(0:8),q(0:8),eq(0:8),neq(0:8),ss(0:8),hh(0:8),sp(0:8),om(0:8)
    real(kind=8) :: ux,uy,phi,cu,jx,jy,num,den,lambdaOracle,lv,av,err,forceY
    real(kind=8) :: beforeMass,beforeJx,beforeJy,expected,initialHeat
    call acc_init(acc_device_host)
    call initial()
    call require(all(ex==cx).and.all(ey==cy),'D2Q9_direction_order')
    call require(maxval(abs(omegaT-wi))<1d-15,'D2Q9_weight_definition')

    ! Eq.7 的矩定义及严格平衡态退化：调用生产 equilibrium/collision。
    do n=1,8
      ux=0.11d0*cos(dble(n));uy=0.08d0*sin(dble(n));phi=0.7d0+0.2d0*dble(n)
      call feng_thermal_equilibrium(phi,ux,uy,p)
      call measure('equilibrium_mass',abs(sum(p)-phi),5d-14)
      call measure('equilibrium_flux_x',abs(sum(p*dble(cx))-phi*ux),5d-14)
      call measure('equilibrium_flux_y',abs(sum(p*dble(cy))-phi*uy),5d-14)
      call feng_thermal_collision(p,ux,uy,q,lv,av,status)
      call require(status==0,'equilibrium_positive_branch')
      call measure('equilibrium_fixed_point',maxval(abs(q-p)),5d-14)
      call measure('equilibrium_fallback_unneeded',abs(av-1d0),1d-14)
    enddo

    ! 独立逐方向求和 oracle：Eq.7、12、34、38、39，不复用生产分解例程。
    fallbackCases=0
    do n=1,400
      ux=0.08d0*cos(dble(n));uy=0.06d0*sin(dble(n))
      p=wi*(1d0+0.88d0*sin(dble(n)*0.719d0+[(dble(k)*1.713d0,k=0,8)]))
      p=2d0*p/sum(p);phi=sum(p)
      do k=0,8
        cu=dble(cx(k))*ux+dble(cy(k))*uy
        eq(k)=phi*wi(k)*(1d0+3d0*cu+4.5d0*cu**2-1.5d0*(ux**2+uy**2))
      enddo
      neq=p-eq;jx=sum(dble(cx)*neq);jy=sum(dble(cy)*neq)
      do k=0,8
        ss(k)=3d0*wi(k)*(dble(cx(k))*jx+dble(cy(k))*jy)
      enddo
      hh=neq-ss;sp=(1d0-fengFluxRate)*ss
      den=sum(hh*hh/eq)
      num=sum(hh*sp/eq)+sum(hh*log(eq/wi))
      lambdaOracle=1d0+num/den
      om=-fengFluxRate*ss-lambdaOracle*hh
      call feng_thermal_collision(p,ux,uy,q,lv,av,status)
      call require(status==0,'nonequilibrium_positive_branch')
      call measure('ADE_lambda_with_J',abs(lv-lambdaOracle),2d-11)
      call measure('thermal_collision_scalar_conservation',abs(sum(q)-phi),5d-13)
      call measure('full_increment_alpha',maxval(abs(q-p-av*om)),8d-13)
      call measure('post_x_flux_including_fallback', &
          abs(sum(dble(cx)*(q-eq))-(1d0-av*fengFluxRate)*jx),8d-13)
      call measure('post_y_flux_including_fallback', &
          abs(sum(dble(cy)*(q-eq))-(1d0-av*fengFluxRate)*jy),8d-13)
      call require(minval(q)>0d0,'post_collision_population_positive')
      call require(av>0d0.and.av<=1d0,'fallback_alpha_range')
      if(av<1d0-1d-12) fallbackCases=fallbackCases+1
    enddo
    call require(fallbackCases>0,'fallback_test_is_nontrivial')
    write(*,'(a,i0)') 'FALLBACK_CASES ',fallbackCases
    p=wi;p(1)=-0.001d0
    call feng_thermal_collision(p,0d0,0d0,q,lv,av,status)
    call require(status/=0,'negative_population_is_rejected')
    p=-wi
    call feng_thermal_collision(p,0d0,0d0,q,lv,av,status)
    call require(status/=0,'negative_temperature_is_rejected')

    ! Luo 流场力项：只用质量/动量矩定义作 oracle。
    p=wi*(1d0+0.05d0*[(sin(dble(k+1)),k=0,8)])
    p=p/sum(p)
    beforeMass=sum(p);beforeJx=sum(p*dble(cx));beforeJy=sum(p*dble(cy))
    rho=beforeMass;T=Tref+0.2d0
    forceY=beforeMass*gBeta*0.2d0
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

    ! 两个方向的正温标纯导热：实际碰撞、迁移、半程壁面、Nu 后处理。
    u=0d0;v=0d0
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
        expected=Thot+(dble(i)-0.5d0)/dble(nx)*(Tcold-Thot)
#else
        expected=Thot+(dble(j)-0.5d0)/dble(ny)*(Tcold-Thot)
#endif
        err=max(err,abs(T(i,j)-expected))
      enddo
    enddo
    call measure('linear_conduction_temperature',err,2d-10)
    call measure('linear_conduction_heat_balance',abs(sum(T)-initialHeat),2d-9)
    call require(minval(g)>0d0,'linear_conduction_positive_populations')
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


def inspect_entropy_history(directory: Path, meta: dict[str, str]) -> int:
    lines = (directory / ENTROPY_HISTORY).read_text().splitlines()
    if not lines or not lines[0].startswith("# itc t_ff fallback_count "):
        raise AssertionError("Entropy history header is missing")
    rows = [line.split() for line in lines[1:] if line.strip()]
    if len(rows) != int(meta["nureSampleCount"]):
        raise AssertionError("Entropy history count differs from Nu/Re/checkpoint")
    total_step = int(meta["itc_total"])
    time_unit = total_step / float(meta["time_tf"])
    previous_step, previous_time, previous_count = -1, -1.0, -1
    for row in rows:
        if len(row) != 10:
            raise AssertionError("Invalid entropy history column count")
        step, count = int(row[0]), int(row[2])
        values = [float(value) for value in row[1:2] + row[3:]]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("Nonfinite entropy history")
        time_tf, fraction, alpha_min, *_ = values
        if not (previous_step < step <= total_step and time_tf > previous_time):
            raise AssertionError("Entropy history steps/times are not ordered within checkpoint")
        if abs(time_tf - step / time_unit) > 1e-12 * max(1.0, abs(time_tf)):
            raise AssertionError("Entropy history step and time disagree")
        if not (0 <= count <= NX * NY * step and count >= previous_count and 0.0 < alpha_min <= 1.0):
            raise AssertionError("Invalid cumulative entropy diagnostics")
        if abs(fraction - count / (NX * NY * max(1, step))) > 1e-13:
            raise AssertionError("Entropy fallback fraction differs from cumulative count")
        previous_step, previous_time, previous_count = step, time_tf, count
    return len(rows)


def rejected_entropy_restarts(seed: Path, executable: Path, root: Path,
                              env: dict[str, str], timeout: float) -> dict:
    """破坏真实32步检查点配套历史，验证生产读取器拒绝而非继续拼接。"""
    entropy_bytes = (seed / ENTROPY_HISTORY).read_bytes()
    lines = entropy_bytes.decode().splitlines(keepends=True)
    row_indices = [index for index, line in enumerate(lines)
                   if line.strip() and not line.lstrip().startswith("#")]
    if not row_indices:
        raise AssertionError("Restart rejection seed must contain entropy samples")
    _, meta = metadata(seed)
    total_step = int(meta["itc_total"])
    time_unit = total_step / float(meta["time_tf"])
    cases = {}
    for label in ("missing", "short", "ahead", "duplicate_tail", "wrong_time", "wrong_step"):
        directory = root / label
        shutil.copytree(seed, directory)
        path = directory / ENTROPY_HISTORY
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
                     or p.name == ENTROPY_HISTORY
                     or p.name.startswith(("Nu_VolAvg", "Re_VolAvg"))}
        result = run_command([str(executable)], directory, env, timeout, expect_success=False)
        output = result.stdout + result.stderr
        (directory / "rejected_restart.log").write_text(output, encoding="utf-8")
        if result.returncode == 0 or "EntropicADE restart entropy" not in output:
            raise AssertionError(f"Entropy {label} restart was not refused by the production validator: "
                                 + output[-4000:])
        if any((directory / filename).read_bytes() != content
               for filename, content in protected.items()):
            raise AssertionError(f"Rejected entropy {label} restart changed saved state/history")
        if label == "missing" and path.exists():
            raise AssertionError("Rejected missing entropy history was silently recreated")
        cases[label] = {"rejected_by_entropy_validator": True,
                        "checkpoint_and_histories_preserved": True}
    return cases


def inspect_run(directory: Path, expected_magic: str) -> dict:
    state, meta = metadata(directory)
    records = read_records(state)
    if len(records[0]) != 48 or records[0][:32].rstrip().decode() != expected_magic:
        raise AssertionError("Wrong scheme identifier in checkpoint")
    if struct.unpack("<4i", records[0][32:]) != (NX, NY, 9, 9):
        raise AssertionError("Wrong checkpoint mesh or velocity counts")
    if any(len(r) != NX * NY * 9 * 8 for r in records[1:3]):
        raise AssertionError("Checkpoint does not contain complete f/g")
    flow, thermal = [struct.unpack(f"<{NX*NY*9}d", r) for r in records[1:3]]
    if not all(math.isfinite(v) for v in (*flow, *thermal)) or min(thermal) <= 0:
        raise AssertionError("Nonfinite f/g or nonpositive thermal population")
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
            "entropy_history_samples": inspect_entropy_history(directory, meta),
            "checkpoint_record_bytes": [len(r) for r in records]}


def integration_checks(source: str, root: Path, compiler: Path,
                       env: dict[str, str], timeout: float) -> dict:
    match = re.search(r"reloadMagic\s*=\s*['\"]([^'\"]+)['\"]", source, re.I)
    if not match:
        raise AssertionError("Scheme-specific reloadMagic missing")
    magic = match[1]
    results = {}
    for name, steady, rb in (("side_steady", True, False), ("rb_unsteady", False, True)):
        full = root / name
        build_run(small_source(source, steady=steady, rb=rb, steps=48), full, compiler, env, timeout)
        results[name] = inspect_run(full, magic)
        split = root / (name + "_restart")
        build_run(small_source(source, steady=steady, rb=rb, steps=32), split, compiler, env, timeout)
        seed = root / (name + "_restart_seed32")
        shutil.copytree(split, seed)
        histories = {p.name: p.read_bytes() for prefix in ("Nu_VolAvg", "Re_VolAvg")
                     for p in split.glob(prefix + "*.dat")}
        histories[ENTROPY_HISTORY] = (split / ENTROPY_HISTORY).read_bytes()
        build_run(small_source(source, steady=steady, rb=rb, steps=48, restart=True),
                  split, compiler, env, timeout)
        after = inspect_run(split, magic)
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
        if not steady and (split / ENTROPY_HISTORY).read_bytes() != (full / ENTROPY_HISTORY).read_bytes():
            raise AssertionError("RB continuous/split entropy history bytes differ")
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
            "entropy_history_preserved": True, "entropy_final_row_bitwise_equal": True,
            "entropy_full_history_bitwise_equal": not steady,
            "entropy_history_comparison_note": ("Steady stop at32 adds a final sample; compare old prefix and final48 row"
                                                if steady else "Every RB sample and header are bitwise equal"),
            "completed_target_restart_preserved": True, "after": after}
        results[name + "_invalid_entropy_restarts"] = rejected_entropy_restarts(
            seed, split / "solver.exe", root / (name + "_invalid_entropy"), env, timeout)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "2DRBOpenaccEntropicADE.F90")
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
    root = Path(tempfile.mkdtemp(prefix="feng-ade-local-", dir=args.temp_root)).resolve()
    if not str(root).isascii():
        parser.error("Use an ASCII --temp-root for portable Fortran builds")
    env = compiler_environment(compiler)
    report = {"passed": False, "verified_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "source": str(args.source.resolve()), "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
              "compiler": str(compiler), "backend": "GNU Fortran OpenACC host; NVIDIA GPU not tested",
              "build_directory": str(root), "high_Ra_accuracy_validated": False,
              "grid_independence_validated": False,
              "scope": "Actual production kernels; moment/positive-branch/conduction and finite-step output/restart checks"}
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
        for label, rb in (("side_heated", False), ("rayleigh_benard", True)):
            clone = small_source(source, steady=True, rb=rb, steps=48)
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
