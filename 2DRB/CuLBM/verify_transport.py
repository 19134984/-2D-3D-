"""用实际 CuLBM 双分布推进核检验运动黏度和热扩散率，不改生产源码。

64x6 周期测试由 driver 补齐缺失链，两例均调用实际 collision/collisionT、
streaming/streamingT、macroT/macro。测试为正温度载体上的小幅线性波；
gfortran OpenACC host 通过不代表 NVIDIA GPU 或高 Ra 热对流验证。
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
COEFFICIENT = 0.05
AMPLITUDE = 1.0e-5
BACKGROUND = 1.5


DRIVER = r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,alpha,step
    integer,parameter :: transportCase=__CASE__
    real(kind=8) :: eq(0:8),eqT(0:4),waveNumber,x
    real(kind=8) :: amplitudeSin,amplitudeCos,fieldValue,firstRowValue
    real(kind=8) :: maxDensityError,maxUnwantedVelocity,maxRowDifference,maxOtherFieldError
    real(kind=8),parameter :: initialAmplitude=1.0d-5

    if(nx.NE.64.OR.ny.NE.6) error stop 'Unexpected transport-test grid'
    if(abs(viscosity-0.05d0).GT.1.0d-14) error stop 'Unexpected transport-test viscosity'
    if(abs(diffusivity-0.05d0).GT.1.0d-14) error stop 'Unexpected transport-test diffusivity'
    if(abs(taug-0.75d0).GT.1.0d-14) error stop 'Unexpected transport-test thermal rate'
    if(abs(Tref-1.5d0).GT.1.0d-14) error stop 'Unexpected positive temperature background'
    if(gBeta.NE.0.0d0) error stop 'Transport-test force must vanish'
    call acc_init(acc_device_host)
    allocate(f(nx,ny,0:8),f_post(0:nx+1,0:ny+1,0:8))
    allocate(g(nx,ny,0:4),g_post(0:nx+1,0:ny+1,0:4))
    allocate(rho(nx,ny),u(nx,ny),v(nx,ny),T(nx,ny),Fx(nx,ny),Fy(nx,ny))
    waveNumber=2.0d0*acos(-1.0d0)/dble(nx)
    rho=1.0d0;u=0.0d0;v=0.0d0;T=Tref;Fx=0.0d0;Fy=0.0d0
    cuCollisionStatus=0
    do j=1,ny
      do i=1,nx
        x=dble(i)-0.5d0
        if(transportCase.EQ.1) then
          v(i,j)=initialAmplitude*sin(waveNumber*x)
        else
          T(i,j)=Tref+initialAmplitude*sin(waveNumber*x)
        endif
        call cu_initial_population(1.0d0,u(i,j),v(i,j),0.0d0,0.0d0,eq)
        call cu_thermal_equilibrium(T(i,j),u(i,j),v(i,j),eqT)
        do alpha=0,8
          f(i,j,alpha)=eq(alpha)
        enddo
        do alpha=0,4
          g(i,j,alpha)=eqT(alpha)
        enddo
      enddo
    enddo
    !$acc enter data copyin(f,g,rho,u,v,T,Fx,Fy,ex,ey,cuCollisionStatus) create(f_post,g_post)
    call sample_mode(0)
    do step=1,__STEPS__
      call collision()
      call collisionT()
      call streaming()
      call streamingT()
      call periodic_missing_links_driver()
      call macroT()
      call macro()
      if(mod(step,20).EQ.0.OR.step.EQ.__STEPS__) call sample_mode(step)
    enddo
    !$acc wait(1)
    !$acc exit data delete(f,f_post,g,g_post,rho,u,v,T,Fx,Fy,ex,ey,cuCollisionStatus)
    deallocate(f,f_post,g,g_post,rho,u,v,T,Fx,Fy)
    call acc_shutdown(acc_device_host)

contains
    subroutine sample_mode(stepNumber)
      integer,intent(in) :: stepNumber
      !$acc wait(1)
      !$acc update self(rho,u,v,T,cuCollisionStatus)
      if(cuCollisionStatus.NE.0) error stop 'CuLBM collision rejected the test state'
      amplitudeSin=0.0d0;amplitudeCos=0.0d0;maxRowDifference=0.0d0
      do j=1,ny
        do i=1,nx
          x=dble(i)-0.5d0
          if(transportCase.EQ.1) then
            fieldValue=v(i,j);firstRowValue=v(i,1)
          else
            fieldValue=T(i,j)-Tref;firstRowValue=T(i,1)-Tref
          endif
          amplitudeSin=amplitudeSin+fieldValue*sin(waveNumber*x)
          amplitudeCos=amplitudeCos+fieldValue*cos(waveNumber*x)
          maxRowDifference=max(maxRowDifference,abs(fieldValue-firstRowValue))
        enddo
      enddo
      amplitudeSin=2.0d0*amplitudeSin/dble(nx*ny)
      amplitudeCos=2.0d0*amplitudeCos/dble(nx*ny)
      maxDensityError=maxval(abs(rho-1.0d0))
      if(transportCase.EQ.1) then
        maxUnwantedVelocity=maxval(abs(u))
        maxOtherFieldError=maxval(abs(T-Tref))
      else
        maxUnwantedVelocity=max(maxval(abs(u)),maxval(abs(v)))
        maxOtherFieldError=maxUnwantedVelocity
      endif
      write(*,'(I7,1x,10(ES25.16E3,1x))') stepNumber,amplitudeSin,amplitudeCos, &
          maxDensityError,maxUnwantedVelocity,maxRowDifference,sum(rho),minval(rho), &
          sum(T),minval(T),maxOtherFieldError
    end subroutine sample_mode

    subroutine periodic_missing_links_driver()
      integer :: ii,jj,aa,ip,jp,wrappedI,wrappedJ
      !$acc parallel loop gang vector collapse(2) default(none) &
      !$acc& present(f,f_post,g,g_post,ex,ey) async(1) private(aa,ip,jp,wrappedI,wrappedJ)
      do jj=1,ny
        do ii=1,nx
          do aa=0,8
            ip=ii-ex(aa);jp=jj-ey(aa)
            if(ip.LT.1.OR.ip.GT.nx.OR.jp.LT.1.OR.jp.GT.ny) then
              wrappedI=modulo(ip-1,nx)+1
              wrappedJ=modulo(jp-1,ny)+1
              f(ii,jj,aa)=f_post(wrappedI,wrappedJ,aa)
              if(aa.LE.4) g(ii,jj,aa)=g_post(wrappedI,wrappedJ,aa)
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
    changed, count = re.subn(pattern, lambda m: m[1] + value, source, count=1)
    if count != 1:
        raise AssertionError(f"Expected parameter {name}")
    return changed


def make_driver_source(source: str, steps: int, case: int) -> str:
    for name, value in (("nx", str(NX)), ("ny", str(NY)), ("Prandtl", "1.0d0"),
                        ("tauf", "0.5d0+3.0d0*0.05d0"), ("gBeta", "0.0d0"),
                        ("timeUnit", "1.0d0"), ("velocityUnit", "1.0d0")):
        source = replace_parameter(source, name, value)
    driver = DRIVER.replace("__STEPS__", str(steps)).replace("__CASE__", str(case))
    source, count = re.subn(
        r"(?ims)^\s*program\s+main\s*$.*?^\s*end\s+program\s+main\s*$",
        lambda _: "\n" + driver + "\n", source, count=1)
    if count != 1:
        raise AssertionError("Expected exactly one production program main")
    return source


def run_logged(command: list[str], cwd: Path, env: dict[str, str], log_name: str) -> str:
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=120)
    (cwd / log_name).write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}); inspect {cwd / log_name}")
    return result.stdout


def analyze(output: str, case_name: str, steps: int, fit_start: int, tolerance: float) -> dict:
    rows = []
    for line in output.splitlines():
        values = line.split()
        if len(values) != 11:
            raise AssertionError(f"Unexpected test output: {line}")
        row = [float(value.replace("D", "E")) for value in values]
        if not all(math.isfinite(value) for value in row):
            raise AssertionError("Nonfinite transport sample")
        rows.append(row)
    if not rows or rows[0][0] != 0 or rows[-1][0] != steps:
        raise AssertionError("Transport driver did not reach its full horizon")
    fit = [row for row in rows if row[0] >= fit_start]
    if len(fit) < 5 or any(row[1] <= 0.0 for row in rows):
        raise AssertionError("Invalid positive-amplitude decay fit")
    mean_t = sum(row[0] for row in fit) / len(fit)
    mean_log = sum(math.log(row[1]) for row in fit) / len(fit)
    slope = sum((row[0] - mean_t) * (math.log(row[1]) - mean_log) for row in fit)
    slope /= sum((row[0] - mean_t) ** 2 for row in fit)
    wave_number = 2.0 * math.pi / NX
    fitted = -slope / wave_number ** 2
    relative_error = abs(fitted / COEFFICIENT - 1.0)
    predicted_final = AMPLITUDE * math.exp(-COEFFICIENT * wave_number ** 2 * steps)
    endpoint_error = abs(rows[-1][1] / predicted_final - 1.0)
    max_mass_drift = max(abs(row[6] - NX * NY) / (NX * NY) for row in rows)
    max_heat_drift = max(abs(row[8] - BACKGROUND * NX * NY) / (BACKGROUND * NX * NY) for row in rows)
    max_phase = max(abs(row[2]) for row in rows) / AMPLITUDE
    max_density = max(row[3] for row in rows)
    max_unwanted_velocity = max(row[4] for row in rows)
    max_row_difference = max(row[5] for row in rows)
    max_other_error = max(row[10] for row in rows)
    checks = {
        "full_horizon_finite_and_collision_accepted": True,
        "transport_coefficient_error_below_tolerance": relative_error < tolerance,
        "analytical_endpoint_error_below_tolerance": endpoint_error < tolerance,
        "mass_conservation": max_mass_drift < 1.0e-11,
        "heat_conservation": max_heat_drift < 1.0e-11,
        "density_uniformity": max_density < 1.0e-10,
        "no_unwanted_velocity": max_unwanted_velocity < 1.0e-10,
        "no_transverse_row_variation": max_row_difference < 1.0e-12,
        "no_phase_shift": max_phase < 1.0e-7,
        "other_field_preserved": max_other_error < 1.0e-10,
        "positive_density": min(row[7] for row in rows) > 0.0,
        "positive_temperature_carrier": min(row[9] for row in rows) > 0.0,
    }
    return {
        "case": case_name, "passed": all(checks.values()), "checks": checks,
        "grid": [NX, NY], "target_coefficient": COEFFICIENT,
        "coefficient": "kinematic_viscosity" if case_name == "shear_wave" else "thermal_diffusivity",
        "temperature_background": BACKGROUND, "initial_amplitude": AMPLITUDE,
        "steps": steps, "wave_number": wave_number, "fit_start_step": fit_start,
        "fit_samples": len(fit), "fitted_coefficient": fitted,
        "relative_coefficient_error": relative_error, "relative_error_tolerance": tolerance,
        "initial_measured_amplitude": rows[0][1], "final_measured_amplitude": rows[-1][1],
        "final_analytical_amplitude": predicted_final, "relative_endpoint_error": endpoint_error,
        "max_relative_mass_drift": max_mass_drift, "max_relative_heat_drift": max_heat_drift,
        "max_density_error": max_density, "max_unwanted_velocity": max_unwanted_velocity,
        "max_row_difference": max_row_difference, "max_cosine_amplitude_relative": max_phase,
        "max_other_field_error": max_other_error, "minimum_density": min(row[7] for row in rows),
        "minimum_temperature": min(row[9] for row in rows),
        "samples": [{"step": int(row[0]), "sin_amplitude": row[1], "cos_amplitude": row[2]}
                    for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "2DRBOpenaccCuLBM.F90")
    parser.add_argument("--compiler", type=Path)
    parser.add_argument("--temp-root", type=Path, help="Use an ASCII directory for Windows gfortran")
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
        found = shutil.which("gfortran")
        compiler = Path(found) if found else Path("C:/msys64/ucrt64/bin/gfortran.exe")
    if not compiler.is_file():
        parser.error(f"gfortran compiler not found: {compiler}")
    compiler = compiler.resolve()
    env = os.environ.copy()
    env["PATH"] = str(compiler.parent) + os.pathsep + env.get("PATH", "")
    env["ACC_DEVICE_TYPE"] = "host"; env["LANG"] = "C"; env["LC_ALL"] = "C"
    source_path = args.source.resolve()
    source_sha = sha256(source_path)
    source = source_path.read_text(encoding="utf-8-sig")
    if args.temp_root is not None:
        args.temp_root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="culbm-transport-", dir=args.temp_root)).resolve()
    if not str(work).isascii():
        parser.error("The build directory must be ASCII; specify --temp-root")
    report = {
        "passed": False, "verified_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source": str(source_path), "source_sha256": source_sha, "compiler": str(compiler),
        "backend": "gfortran OpenACC host", "build_directory": str(work),
        "production_subroutines": ["collision", "collisionT", "streaming", "streamingT", "macroT", "macro"],
        "test_only_changes": "64x6; Pr=1; nu=kappa=0.05; zero force; main replaced; periodic missing-link driver",
        "limitations": "Small-amplitude zero-force linear transport on positive T only; not NVIDIA GPU, high-Ra stability, nonlinear thermal accuracy or paper benchmark agreement",
        "cases": [],
    }
    started = time.perf_counter()
    try:
        version = run_logged([str(compiler), "--version"], work, env, "compiler_version.log")
        report["compiler_version"] = version.splitlines()[0]
        flags = ["-cpp", "-fopenacc", "-O2", "-ffree-line-length-none", "-fcheck=all",
                 "-ffpe-trap=invalid,zero,overflow", "-fbacktrace"]
        report["compiler_flags"] = flags
        for case, name in ((1, "shear_wave"), (2, "thermal_diffusion_wave")):
            case_dir = work / name; case_dir.mkdir()
            test_source = make_driver_source(source, args.steps, case)
            (case_dir / "transport.F90").write_text(test_source, encoding="utf-8")
            executable = case_dir / ("transport.exe" if os.name == "nt" else "transport")
            run_logged([str(compiler), *flags, "transport.F90", "-o", str(executable)],
                       case_dir, env, "build.log")
            output = run_logged([str(executable)], case_dir, env, "samples.log")
            result = analyze(output, name, args.steps, args.fit_start, args.relative_tolerance)
            result["test_source_sha256"] = hashlib.sha256(test_source.encode("utf-8")).hexdigest()
            report["cases"].append(result)
        report["passed"] = all(case["passed"] for case in report["cases"])
    except Exception as error:
        report["error"] = str(error)
    finally:
        report["wall_seconds"] = time.perf_counter() - started
        report["source_unchanged"] = sha256(source_path) == source_sha
        report["passed"] = bool(report["passed"] and report["source_unchanged"])
        args.output.resolve().write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    compact = {key: value for key, value in report.items() if key != "cases"}
    compact["cases"] = [{key: value for key, value in case.items() if key != "samples"}
                        for case in report["cases"]]
    print(json.dumps(compact, indent=2, ensure_ascii=False))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
