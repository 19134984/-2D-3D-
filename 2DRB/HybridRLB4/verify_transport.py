"""以周期剪切波检验实际 RLB-4 流场核恢复的运动黏度，不修改生产源码。

测试副本仅替换参数和主程序，调用实际 collision / streaming / macro。
周期缺失链由测试 driver 补齐，不调用生产壁面边界、不增加生产周期宏。
使用 gfortran OpenACC host；这不是 NVIDIA GPU、热方程或高 Ra 的验证。
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
import subprocess
import tempfile
import time


HERE = Path(__file__).resolve().parent
NX, NY = 64, 6
VISCOSITY = 0.05
AMPLITUDE = 1.0e-5


DRIVER = r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,alpha,step
    real(kind=8) :: eq(0:8),waveNumber,x,amplitudeSin,amplitudeCos
    real(kind=8) :: maxDensityError,maxLongitudinalVelocity,maxRowDifference,totalMass
    real(kind=8),parameter :: initialAmplitude=1.0d-5

    if(nx.NE.64.OR.ny.NE.6) error stop 'Unexpected transport-test grid'
    if(abs(viscosity-0.05d0).GT.1.0d-14) error stop 'Unexpected transport-test viscosity'
    if(gBeta.NE.0.0d0) error stop 'Transport-test force must vanish'
    call acc_init(acc_device_host)
    allocate(f(nx,ny,0:8),f_post(0:nx+1,0:ny+1,0:8))
    allocate(rho(nx,ny),u(nx,ny),v(nx,ny),T(nx,ny),Fx(nx,ny),Fy(nx,ny))
    waveNumber=2.0d0*acos(-1.0d0)/dble(nx)
    rho=1.0d0;u=0.0d0;T=Tref;Fx=0.0d0;Fy=0.0d0
    do j=1,ny
      do i=1,nx
        x=dble(i)-0.5d0
        v(i,j)=initialAmplitude*sin(waveNumber*x)
        call rlb_equilibrium(1.0d0,u(i,j),v(i,j),eq)
        do alpha=0,8
          f(i,j,alpha)=eq(alpha)
        enddo
      enddo
    enddo
    !$acc enter data copyin(f,rho,u,v,T,Fx,Fy,ex,ey) create(f_post)
    call sample_mode(0)
    do step=1,__STEPS__
      call collision()
      call streaming()
      call periodic_missing_links_driver()
      call macro()
      if(mod(step,20).EQ.0.OR.step.EQ.__STEPS__) call sample_mode(step)
    enddo
    !$acc wait(1)
    !$acc exit data delete(f,f_post,rho,u,v,T,Fx,Fy,ex,ey)
    deallocate(f,f_post,rho,u,v,T,Fx,Fy)
    call acc_shutdown(acc_device_host)

contains
    subroutine sample_mode(stepNumber)
      integer,intent(in) :: stepNumber
      !$acc wait(1)
      !$acc update self(rho,u,v)
      amplitudeSin=0.0d0;amplitudeCos=0.0d0
      do j=1,ny
        do i=1,nx
          x=dble(i)-0.5d0
          amplitudeSin=amplitudeSin+v(i,j)*sin(waveNumber*x)
          amplitudeCos=amplitudeCos+v(i,j)*cos(waveNumber*x)
        enddo
      enddo
      amplitudeSin=2.0d0*amplitudeSin/dble(nx*ny)
      amplitudeCos=2.0d0*amplitudeCos/dble(nx*ny)
      maxDensityError=maxval(abs(rho-1.0d0))
      maxLongitudinalVelocity=maxval(abs(u))
      maxRowDifference=0.0d0
      do j=2,ny
        maxRowDifference=max(maxRowDifference,maxval(abs(v(:,j)-v(:,1))))
      enddo
      totalMass=sum(rho)
      write(*,'(I7,1x,7(ES25.16E3,1x))') stepNumber,amplitudeSin,amplitudeCos, &
          maxDensityError,maxLongitudinalVelocity,maxRowDifference,totalMass,minval(rho)
    end subroutine sample_mode

    subroutine periodic_missing_links_driver()
      integer :: ii,jj,aa,ip,jp,wrappedI,wrappedJ
      !$acc parallel loop gang vector collapse(2) default(none) &
      !$acc& present(f,f_post,ex,ey) async(1) private(aa,ip,jp,wrappedI,wrappedJ)
      do jj=1,ny
        do ii=1,nx
          do aa=0,8
            ip=ii-ex(aa);jp=jj-ey(aa)
            if(ip.LT.1.OR.ip.GT.nx.OR.jp.LT.1.OR.jp.GT.ny) then
              wrappedI=modulo(ip-1,nx)+1
              wrappedJ=modulo(jp-1,ny)+1
              f(ii,jj,aa)=f_post(wrappedI,wrappedJ,aa)
            endif
          enddo
        enddo
      enddo
    end subroutine periodic_missing_links_driver
end program main
"""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_parameter(source: str, name: str, value: str) -> str:
    pattern = rf"(?im)(\b{re.escape(name)}\s*=\s*)[^,!\n]+"
    result, count = re.subn(pattern, lambda match: match[1] + value, source, count=1)
    if count != 1:
        raise AssertionError(f"Expected parameter {name}")
    return result


def make_driver_source(source: str, steps: int) -> str:
    for name, value in (("nx", str(NX)), ("ny", str(NY)),
                        ("tauf", "0.5d0+3.0d0*0.05d0"),
                        ("gBeta", "0.0d0"), ("timeUnit", "1.0d0"),
                        ("velocityUnit", "1.0d0")):
        source = replace_parameter(source, name, value)
    source, count = re.subn(
        r"(?ims)^\s*program\s+main\s*$.*?^\s*end\s+program\s+main\s*$",
        lambda _: "\n" + DRIVER.replace("__STEPS__", str(steps)) + "\n",
        source, count=1)
    if count != 1:
        raise AssertionError("Expected exactly one production program main")
    return source


def compiler_environment(compiler: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = str(compiler.parent) + os.pathsep + env.get("PATH", "")
    env["ACC_DEVICE_TYPE"] = "host"
    env["LANG"] = "C"
    env["LC_ALL"] = "C"
    return env


def run_logged(command: list[str], cwd: Path, env: dict[str, str],
               log_name: str, timeout: float = 120.0) -> str:
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=timeout)
    (cwd / log_name).write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}); inspect {cwd / log_name}")
    return result.stdout


def analyze(output: str, steps: int, fit_start: int, tolerance: float) -> dict:
    rows = []
    for line in output.splitlines():
        values = line.split()
        if len(values) != 8:
            raise AssertionError(f"Unexpected driver output: {line}")
        row = [float(value.replace("D", "E")) for value in values]
        if not all(math.isfinite(value) for value in row):
            raise AssertionError("Nonfinite shear-wave sample")
        rows.append(row)
    if not rows or rows[0][0] != 0 or rows[-1][0] != steps:
        raise AssertionError("Shear-wave driver did not reach its full horizon")
    fit = [row for row in rows if row[0] >= fit_start]
    if len(fit) < 5 or any(row[1] <= 0.0 for row in rows):
        raise AssertionError("Invalid positive-amplitude decay fit")
    mean_t = sum(row[0] for row in fit) / len(fit)
    mean_log = sum(math.log(row[1]) for row in fit) / len(fit)
    denominator = sum((row[0] - mean_t) ** 2 for row in fit)
    slope = sum((row[0] - mean_t) * (math.log(row[1]) - mean_log)
                for row in fit) / denominator
    wave_number = 2.0 * math.pi / NX
    fitted_viscosity = -slope / wave_number ** 2
    relative_error = abs(fitted_viscosity / VISCOSITY - 1.0)
    predicted_final = AMPLITUDE * math.exp(-VISCOSITY * wave_number ** 2 * steps)
    endpoint_error = abs(rows[-1][1] / predicted_final - 1.0)
    max_mass_drift = max(abs(row[6] - NX * NY) / (NX * NY) for row in rows)
    max_phase = max(abs(row[2]) for row in rows) / AMPLITUDE
    max_density = max(row[3] for row in rows)
    max_longitudinal = max(row[4] for row in rows)
    max_row_difference = max(row[5] for row in rows)
    min_density = min(row[7] for row in rows)
    checks = {
        "full_horizon_finite": True,
        "viscosity_error_below_tolerance": relative_error < tolerance,
        "analytical_endpoint_error_below_tolerance": endpoint_error < tolerance,
        "mass_conservation": max_mass_drift < 1.0e-11,
        "density_uniformity": max_density < 1.0e-10,
        "no_longitudinal_mode": max_longitudinal < 1.0e-11,
        "no_transverse_row_variation": max_row_difference < 1.0e-12,
        "no_phase_shift": max_phase < 1.0e-8,
        "positive_density": min_density > 0.0,
    }
    return {
        "passed": all(checks.values()), "checks": checks,
        "grid": [NX, NY], "target_viscosity": VISCOSITY,
        "initial_amplitude": AMPLITUDE, "wave_number": wave_number,
        "steps": steps, "fit_start_step": fit_start, "fit_samples": len(fit),
        "fitted_viscosity": fitted_viscosity, "relative_viscosity_error": relative_error,
        "relative_error_tolerance": tolerance, "initial_measured_amplitude": rows[0][1],
        "final_measured_amplitude": rows[-1][1], "final_analytical_amplitude": predicted_final,
        "relative_endpoint_error": endpoint_error, "max_relative_mass_drift": max_mass_drift,
        "max_density_error": max_density, "max_longitudinal_velocity": max_longitudinal,
        "max_row_difference": max_row_difference, "max_cosine_amplitude_relative": max_phase,
        "minimum_density": min_density,
        "samples": [{"step": int(row[0]), "sin_amplitude": row[1], "cos_amplitude": row[2]}
                    for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "2DRBOpenaccHybridRLB4.F90")
    parser.add_argument("--compiler", type=Path)
    parser.add_argument("--temp-root", type=Path,
                        help="Use an ASCII path for gfortran builds on Windows")
    parser.add_argument("--output", type=Path, default=HERE / "transport_validation_local.json")
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--fit-start", type=int, default=200)
    parser.add_argument("--relative-tolerance", type=float, default=0.01)
    args = parser.parse_args()
    if args.steps < 100 or not 0 <= args.fit_start <= args.steps - 80:
        parser.error("Need at least five fit samples before the full horizon")
    if not 0.0 < args.relative_tolerance < 1.0:
        parser.error("relative-tolerance must lie between zero and one")
    compiler = args.compiler
    if compiler is None:
        discovered = shutil.which("gfortran")
        compiler = Path(discovered) if discovered else Path("C:/msys64/ucrt64/bin/gfortran.exe")
    if not compiler.is_file():
        parser.error(f"gfortran compiler not found: {compiler}")
    source_path = args.source.resolve()
    source_sha = sha256(source_path)
    source = source_path.read_text(encoding="utf-8-sig")
    test_source = make_driver_source(source, args.steps)
    if args.temp_root is not None:
        args.temp_root.mkdir(parents=True, exist_ok=True)
    build_dir = Path(tempfile.mkdtemp(prefix="hybrid-transport-", dir=args.temp_root)).resolve()
    if not str(build_dir).isascii():
        parser.error("The build directory must be ASCII; specify --temp-root")
    (build_dir / "shear_wave.F90").write_text(test_source, encoding="utf-8")
    env = compiler_environment(compiler)
    executable = build_dir / ("shear_wave.exe" if os.name == "nt" else "shear_wave")
    report = {
        "passed": False, "verified_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source": str(source_path), "source_sha256": source_sha,
        "test_source_sha256": hashlib.sha256(test_source.encode("utf-8")).hexdigest(),
        "compiler": str(compiler.resolve()), "backend": "gfortran OpenACC host",
        "build_directory": str(build_dir),
        "production_subroutines": ["collision", "streaming", "macro", "rlb_collide_node", "rlb_equilibrium"],
        "test_only_changes": "64x6 grid; nu=0.05; zero force; shear-wave main; periodic missing-link driver",
        "limitations": "Tests zero-force shear viscosity only; not NVIDIA GPU, thermal convection, high-Ra stability or paper benchmark agreement",
    }
    started = time.perf_counter()
    try:
        version = run_logged([str(compiler), "--version"], build_dir, env, "compiler_version.log")
        report["compiler_version"] = version.splitlines()[0]
        flags = ["-cpp", "-fopenacc", "-O2", "-ffree-line-length-none", "-fcheck=all",
                 "-ffpe-trap=invalid,zero,overflow", "-fbacktrace"]
        report["compiler_flags"] = flags
        run_logged([str(compiler), *flags, "shear_wave.F90", "-o", str(executable)],
                   build_dir, env, "build.log")
        output = run_logged([str(executable)], build_dir, env, "shear_wave_samples.log")
        report.update(analyze(output, args.steps, args.fit_start, args.relative_tolerance))
    except Exception as error:
        report["error"] = str(error)
    finally:
        report["wall_seconds"] = time.perf_counter() - started
        report["source_unchanged"] = sha256(source_path) == source_sha
        report["passed"] = bool(report["passed"] and report["source_unchanged"])
        args.output.resolve().write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                                         encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key not in ("samples",)},
                     indent=2, ensure_ascii=False))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
