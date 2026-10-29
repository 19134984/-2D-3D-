"""用生产碰撞核做周期剪切波和温度波的独立输运系数检查。

测试副本只改网格/物性、零浮力并替换主程序；周期补链由测试驱动完成。
不改生产求解器，不将这些小振幅周期结果称作热对流或高 Ra 验证。
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
import tempfile
import time

from verify_local import compiler_environment, replace_main, replace_parameter, run_command

HERE = Path(__file__).resolve().parent
NX, NY = 64, 6
COEFFICIENT = 0.05

DRIVER = r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,k,step
    integer,parameter :: mode=__MODE__
    real(kind=8),parameter :: advection=__ADVECTION__,amp=__AMPLITUDE__
    real(kind=8) :: eq(0:8),waveNumber,x,as,ac,totalScalar,quantity
    call acc_init(acc_device_host)
    call initial()
    if(nx/=64.or.ny/=6) error stop 'Wrong transport grid'
    if(abs(viscosity-0.05d0)>1d-14.or.abs(diffusivity-0.05d0)>1d-14) &
        error stop 'Wrong transport coefficients'
    if(gBeta/=0d0) error stop 'Transport test must have zero force'
    waveNumber=2d0*pi/dble(nx)
    rho=1d0;u=advection;v=0d0;T=Tref;Fx=0d0;Fy=0d0
    do j=1,ny
      do i=1,nx
        x=dble(i)-0.5d0
        if(mode==1) v(i,j)=amp*sin(waveNumber*x)
        if(mode==2) T(i,j)=Tref+amp*sin(waveNumber*x)
        call luo_flow_equilibrium(rho(i,j),u(i,j),v(i,j),eq)
        f(i,j,:)=eq
        call zhang_thermal_equilibrium(T(i,j),u(i,j),v(i,j),eq)
        g(i,j,:)=eq
      enddo
    enddo
    zhangFluxPreviousX=T*u;zhangFluxPreviousY=T*v
    call enter_data_2d_openacc()
    call sample_mode(0)
    do step=1,__STEPS__
      if(mode==1) then
        call collision()
        call streaming()
        call periodic_missing_links_driver(.true.)
        call macro()
      else
        call periodic_scalar_collision_driver()
        call streamingT()
        call periodic_missing_links_driver(.false.)
        call macroT()
      endif
      if(mod(step,20)==0.or.step==__STEPS__) call sample_mode(step)
    enddo
    call exit_data_2d_openacc()
contains
    subroutine sample_mode(stepNumber)
      integer,intent(in) :: stepNumber
      !$acc wait(1)
      !$acc update self(rho,u,v,T,g)
      as=0d0;ac=0d0
      do j=1,ny
        do i=1,nx
          x=dble(i)-0.5d0
          quantity=v(i,j)
          if(mode==2) quantity=T(i,j)-Tref
          as=as+quantity*sin(waveNumber*x)
          ac=ac+quantity*cos(waveNumber*x)
        enddo
      enddo
      as=2d0*as/dble(nx*ny);ac=2d0*ac/dble(nx*ny)
      totalScalar=sum(rho)
      if(mode==2) totalScalar=sum(T)
      write(*,'(a,1x,i7,7(1x,es25.16e3))') 'SAMPLE',stepNumber,as,ac,totalScalar, &
          sum(rho),minval(T),minval(g),minval(rho)
    end subroutine sample_mode
    subroutine periodic_scalar_collision_driver()
      integer :: ii,jj,aa,il,ir,jd,ju
      real(kind=8) :: gx,gy,dx,dy,pop(0:8),post(0:8)
      !$acc parallel loop collapse(2) default(none) &
      !$acc& present(T,u,v,g,g_post,zhangFluxPreviousX,zhangFluxPreviousY) async(1) &
      !$acc& private(aa,il,ir,jd,ju,gx,gy,dx,dy,pop,post)
      do jj=1,ny
        do ii=1,nx
          il=modulo(ii-2,nx)+1;ir=modulo(ii,nx)+1
          jd=modulo(jj-2,ny)+1;ju=modulo(jj,ny)+1
          gx=0.5d0*(T(ir,jj)-T(il,jj));gy=0.5d0*(T(ii,ju)-T(ii,jd))
          dx=T(ii,jj)*u(ii,jj)-zhangFluxPreviousX(ii,jj)
          dy=T(ii,jj)*v(ii,jj)-zhangFluxPreviousY(ii,jj)
          do aa=0,8
            pop(aa)=g(ii,jj,aa)
          enddo
          call zhang_thermal_collision(pop,T(ii,jj),u(ii,jj),v(ii,jj),gx,gy,dx,dy,post)
          do aa=0,8
            g_post(ii,jj,aa)=post(aa)
          enddo
          zhangFluxPreviousX(ii,jj)=T(ii,jj)*u(ii,jj)
          zhangFluxPreviousY(ii,jj)=T(ii,jj)*v(ii,jj)
        enddo
      enddo
    end subroutine periodic_scalar_collision_driver
    subroutine periodic_missing_links_driver(flow)
      logical,intent(in) :: flow
      integer :: ii,jj,aa,ip,jp,wi,wj
      if(flow) then
        !$acc parallel loop collapse(2) default(none) present(f,f_post,ex,ey) async(1) &
        !$acc& private(aa,ip,jp,wi,wj)
        do jj=1,ny
          do ii=1,nx
            do aa=0,8
              ip=ii-ex(aa);jp=jj-ey(aa)
              if(ip<1.or.ip>nx.or.jp<1.or.jp>ny) then
                wi=modulo(ip-1,nx)+1;wj=modulo(jp-1,ny)+1
                f(ii,jj,aa)=f_post(wi,wj,aa)
              endif
            enddo
          enddo
        enddo
      else
        !$acc parallel loop collapse(2) default(none) present(g,g_post,ex,ey) async(1) &
        !$acc& private(aa,ip,jp,wi,wj)
        do jj=1,ny
          do ii=1,nx
            do aa=0,8
              ip=ii-ex(aa);jp=jj-ey(aa)
              if(ip<1.or.ip>nx.or.jp<1.or.jp>ny) then
                wi=modulo(ip-1,nx)+1;wj=modulo(jp-1,ny)+1
                g(ii,jj,aa)=g_post(wi,wj,aa)
              endif
            enddo
          enddo
        enddo
      endif
    end subroutine periodic_missing_links_driver
end program main
"""

GAUSSIAN_MAIN = r"""
program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,k,n
    real(kind=8) :: eq(0:8),r2,width2,analytical,err2,norm2,initialHeat,relativeMass
    call acc_init(acc_device_host)
    call initial()
    if(nx/=400.or.ny/=400.or.zhangEnableFilter/=0) error stop 'Wrong paper Gaussian configuration'
    rho=1d0;u=0.1d0;v=0.1d0;Fx=0d0;Fy=0d0
    do j=1,ny
      do i=1,nx
        r2=(dble(i)-200d0)**2+(dble(j)-200d0)**2
        T(i,j)=exp(-r2/200d0)
        call zhang_thermal_equilibrium(T(i,j),u(i,j),v(i,j),eq)
        g(i,j,:)=eq
      enddo
    enddo
    initialHeat=sum(T)
    zhangFluxPreviousX=T*u;zhangFluxPreviousY=T*v
    call enter_data_2d_openacc()
    do n=1,200
      call periodic_scalar_collision_driver()
      call streamingT()
      call periodic_missing_links_driver(.false.)
      call macroT()
    enddo
    !$acc wait(1)
    !$acc update self(T,g)
    if(.not.all(abs(T)<=huge(1d0)).or..not.all(abs(g)<=huge(1d0))) error stop 'Nonfinite Gaussian state'
    width2=100d0+2d0*diffusivity*200d0
    err2=0d0;norm2=0d0
    do j=1,ny
      do i=1,nx
        r2=(dble(i)-220d0)**2+(dble(j)-220d0)**2
        analytical=100d0/width2*exp(-r2/(2d0*width2))
        err2=err2+(T(i,j)-analytical)**2
        norm2=norm2+analytical**2
      enddo
    enddo
    relativeMass=abs(sum(T)-initialHeat)/initialHeat
    write(*,'(a,5(1x,es25.16e3))') 'GAUSSIAN',sqrt(err2/norm2),relativeMass,minval(T),maxval(T),diffusivity
    call exit_data_2d_openacc()
contains
__PERIODIC_HELPERS__
end program main
"""


def driver_source(source: str, mode: int, advection: float, steps: int) -> str:
    for name, value in (("nx", str(NX)), ("ny", str(NY)), ("Prandtl", "1d0"),
                        ("viscosity", "0.05d0"), ("diffusivity", "0.05d0"),
                        ("tauf", "0.65d0"), ("taug", "0.65d0"), ("zhangEnableFilter", "0"),
                        ("gBeta", "0d0"), ("timeUnit", "1d0"), ("velocityUnit", "1d0")):
        source = replace_parameter(source, name, value)
    amplitude = "1d-5" if mode == 1 else "1d-3"
    driver = DRIVER.replace("__MODE__", str(mode)).replace("__STEPS__", str(steps))
    driver = driver.replace("__ADVECTION__", f"{advection:.8f}d0").replace("__AMPLITUDE__", amplitude)
    return replace_main(source, driver)


def gaussian_source(source: str, diffusivity: float, model: int, beta: float) -> str:
    for name, value in (("nx", "400"), ("ny", "400"), ("Prandtl", "1d0"),
                        ("viscosity", "0.05d0"), ("diffusivity", f"{diffusivity:.12f}d0"),
                        ("tauf", "0.65d0"), ("taug", f"{0.5+3*diffusivity:.12f}d0"),
                        ("gBeta", "0d0"), ("timeUnit", "1d0"), ("velocityUnit", "1d0"),
                        ("zhangEnableFilter", "0"), ("zhangThermalModel", str(model)),
                        ("zhangBeta", f"{beta:.8f}d0")):
        source = replace_parameter(source, name, value)
    helpers = DRIVER[DRIVER.index("    subroutine periodic_scalar_collision_driver()"):DRIVER.rindex("end program main")]
    return replace_main(source, GAUSSIAN_MAIN.replace("__PERIODIC_HELPERS__", helpers))


def analyze_gaussian(output: str, *, pe: int, label: str, reference: float,
                     reference_tolerance: float) -> dict:
    lines = [line for line in output.splitlines() if line.startswith("GAUSSIAN ")]
    if len(lines) != 1:
        raise AssertionError("Gaussian driver did not emit its final analytical error")
    values = [float(value.replace("D", "E")) for value in lines[0].split()[1:]]
    if len(values) != 5 or not all(math.isfinite(value) for value in values):
        raise AssertionError("Invalid Gaussian final metrics")
    l2_error, mass_drift, t_min, t_max, diffusivity = values
    reference_error = abs(l2_error/reference - 1.0)
    checks = {"analytical_relative_L2_below_0p03": l2_error < 0.03,
              "scalar_mass_conservation": mass_drift < 1e-10,
              "Table1_reference_within_requested_tolerance": reference_error < reference_tolerance,
              "paper_diffusivity": abs(diffusivity - 1.0/pe) < 1e-12}
    return {"passed": all(checks.values()), "checks": checks, "case": label,
            "Pe": pe, "grid": [400, 400], "steps": 200, "velocity": [0.1, 0.1],
            "initial_center": [200, 200], "initial_sigma": 10, "filter_enabled": False,
            "initial_previous_Cu": "C(x,0)*u; the first-step dt(Cu) is exactly zero",
            "gradient": "periodic second-order centered difference; dt=dx=dy=1",
            "relative_L2_error": l2_error, "Table1_relative_L2_reference": reference,
            "relative_deviation_from_Table1": reference_error,
            "Table1_relative_tolerance": reference_tolerance,
            "relative_scalar_mass_drift": mass_drift, "T_min": t_min, "T_max": t_max,
            "qualification": "Periodic scalar production kernel, not coupled RB; paper leaves first-step flux history unspecified"}


def analyze(output: str, *, mode: int, advection: float, steps: int,
            fit_start: int, tolerance: float) -> dict:
    rows = []
    for line in output.splitlines():
        if line.startswith("SAMPLE "):
            row = [float(x.replace("D", "E")) for x in line.split()[1:]]
            if len(row) != 8 or not all(math.isfinite(x) for x in row):
                raise AssertionError("Invalid periodic mode sample")
            rows.append(row)
    if not rows or rows[0][0] != 0 or rows[-1][0] != steps:
        raise AssertionError("Periodic driver did not reach the full horizon")
    fit = [r for r in rows if r[0] >= fit_start]
    if len(fit) < 5:
        raise AssertionError("At least five decay fit samples are required")
    amplitudes = [math.hypot(r[1], r[2]) for r in fit]
    if min(amplitudes) <= 0:
        raise AssertionError("Zero mode amplitude cannot be fitted")
    mt = sum(r[0] for r in fit) / len(fit)
    ml = sum(math.log(a) for a in amplitudes) / len(fit)
    slope = sum((r[0]-mt)*(math.log(a)-ml) for r, a in zip(fit, amplitudes)) / sum((r[0]-mt)**2 for r in fit)
    wave_number = 2*math.pi/NX
    fitted = -slope / wave_number**2
    relative_error = abs(fitted/COEFFICIENT - 1)
    initial_amplitude = 1e-5 if mode == 1 else 1e-3
    endpoint = initial_amplitude*math.exp(-COEFFICIENT*wave_number**2*steps)
    endpoint_error = abs(math.hypot(rows[-1][1], rows[-1][2])/endpoint - 1)
    # Signed zero-mean thermal waves have an almost-zero total scalar; normalize
    # conservation by the known initial field scale rather than round-off in its sum.
    scalar_scale = max(abs(rows[0][3]), NX*NY*initial_amplitude)
    scalar_drift = max(abs(r[3]-rows[0][3])/scalar_scale for r in rows)
    mass_drift = max(abs(r[4]-NX*NY)/(NX*NY) for r in rows)
    phase_errors = [abs(math.remainder(math.atan2(r[2], r[1])+wave_number*advection*r[0], 2*math.pi)) for r in rows]
    checks = {"full_horizon_finite": True, "transport_coefficient_error": relative_error < tolerance,
              "analytical_endpoint_error": endpoint_error < tolerance,
              "periodic_scalar_conservation": scalar_drift < 1e-10,
              "flow_mass_conservation": mass_drift < 1e-10,
              "phase_error_below_0p03_radian": max(phase_errors) < 0.03,
              "positive_flow_density": min(r[7] for r in rows) > 0}
    return {"passed": all(checks.values()), "checks": checks,
            "mode": "shear" if mode == 1 else "thermal", "prescribed_advection": advection,
            "grid": [NX, NY], "steps": steps, "fit_start_step": fit_start,
            "target_coefficient": COEFFICIENT, "fitted_coefficient": fitted,
            "relative_coefficient_error": relative_error, "relative_tolerance": tolerance,
            "relative_analytical_endpoint_error": endpoint_error,
            "maximum_relative_scalar_drift": scalar_drift, "maximum_relative_mass_drift": mass_drift,
            "scalar_conservation_normalization": scalar_scale,
            "maximum_phase_error_radians": max(phase_errors),
            "samples": [{"step": int(r[0]), "sin_amplitude": r[1], "cos_amplitude": r[2]} for r in rows]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "2DRBOpenaccZhangHR.F90")
    parser.add_argument("--compiler", type=Path)
    parser.add_argument("--temp-root", type=Path)
    parser.add_argument("--output", type=Path, default=HERE / "transport_validation_local.json")
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--fit-start", type=int, default=200)
    parser.add_argument("--relative-tolerance", type=float, default=0.01)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--skip-gaussian", action="store_true", help="Only run small-amplitude periodic modes")
    parser.add_argument("--paper-reference-tolerance", type=float, default=0.01,
                        help="Allowed relative deviation from the rounded Gaussian errors in Table1")
    args = parser.parse_args()
    if args.steps < args.fit_start+100 or args.fit_start < 0:
        parser.error("Need at least five samples after fit-start")
    compiler = args.compiler or Path(shutil.which("gfortran") or "C:/msys64/ucrt64/bin/gfortran.exe")
    if not compiler.is_file():
        parser.error(f"Compiler not found: {compiler}")
    source_bytes = args.source.read_bytes()
    source = source_bytes.decode("utf-8-sig").replace("\r\n", "\n")
    if args.temp_root:
        args.temp_root.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="zhang-hr-transport-", dir=args.temp_root)).resolve()
    env = compiler_environment(compiler)
    report = {"passed": False, "verified_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "source": str(args.source.resolve()), "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
              "compiler": str(compiler), "backend": "GNU Fortran OpenACC host; not NVIDIA GPU",
              "build_directory": str(root), "cases": [],
              "test_only_changes": "Periodic scalar drivers with zero buoyancy: 64x6 modes at nu=kappa=.05; 400x400 Gaussian at kappa=1/Pe",
              "limitations": "Production local collision with driver periodic gradients; not walls, thermal convection, high-Ra or grid independence"}
    started = time.monotonic()
    flags = ["-cpp", "-fopenacc", "-O2", "-ffree-line-length-none", "-fcheck=all",
             "-ffpe-trap=invalid,zero,overflow", "-fbacktrace"]
    try:
        report["compiler_version"] = run_command([str(compiler), "--version"], root, env, args.timeout).stdout.splitlines()[0]
        report["compiler_flags"] = flags
        for label, mode, advection in (("shear", 1, 0.0), ("thermal_rest", 2, 0.0), ("thermal_advected", 2, 0.03)):
            directory = root / label
            directory.mkdir()
            test_source = driver_source(source, mode, advection, args.steps)
            (directory / "transport.F90").write_text(test_source, encoding="utf-8")
            syntax = run_command([str(compiler), *flags, "-fsyntax-only", "transport.F90"], directory, env, args.timeout)
            compiled = run_command([str(compiler), *flags, "transport.F90", "-o", "transport.exe"], directory, env, args.timeout)
            (directory / "compile.log").write_text(syntax.stdout+syntax.stderr+compiled.stdout+compiled.stderr)
            output = run_command([str(directory / "transport.exe")], directory, env, args.timeout)
            (directory / "samples.log").write_text(output.stdout+output.stderr)
            result = analyze(output.stdout, mode=mode, advection=advection, steps=args.steps,
                             fit_start=args.fit_start, tolerance=args.relative_tolerance)
            result["case"] = label
            report["cases"].append(result)
            if not result["passed"]:
                raise AssertionError(f"Transport case {label} failed acceptance checks")
        report["paper_Gaussian_cases"] = []
        if not args.skip_gaussian:
            for pe, references in ((1, (2.2573e-2, 1.0054e-2, 3.4286e-3, 1.7807e-3)),
                                   (1000, (1.1564e-3, 1.0842e-3, 3.1075e-3, 8.3536e-3))):
                for (label, model, beta), reference in zip(
                        (("BGK", 0, .95), ("PR1", 1, .95),
                         ("HR_beta0p5", 2, .5), ("HR_beta0", 2, 0.0)), references):
                    case = f"Gaussian_Pe{pe}_{label}"
                    directory = root / case
                    directory.mkdir()
                    clone = gaussian_source(source, 1.0/pe, model, beta)
                    (directory / "gaussian.F90").write_text(clone, encoding="utf-8")
                    syntax = run_command([str(compiler), *flags, "-fsyntax-only", "gaussian.F90"], directory, env, args.timeout)
                    compiled = run_command([str(compiler), *flags, "gaussian.F90", "-o", "gaussian.exe"], directory, env, args.timeout)
                    (directory / "compile.log").write_text(syntax.stdout+syntax.stderr+compiled.stdout+compiled.stderr)
                    output = run_command([str(directory / "gaussian.exe")], directory, env, args.timeout)
                    (directory / "samples.log").write_text(output.stdout+output.stderr)
                    result = analyze_gaussian(output.stdout, pe=pe, label=case, reference=reference,
                                              reference_tolerance=args.paper_reference_tolerance)
                    report["paper_Gaussian_cases"].append(result)
                    if not result["passed"]:
                        raise AssertionError(f"Paper Gaussian {case} failed acceptance checks")
        report["passed"] = True
    except Exception as error:
        report["failure"] = str(error)
    finally:
        report["source_unchanged"] = args.source.read_bytes() == source_bytes
        report["passed"] = bool(report["passed"] and report["source_unchanged"])
        report["elapsed_seconds"] = round(time.monotonic()-started, 3)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        print(json.dumps({k: report[k] for k in ("passed", "source_sha256", "build_directory", "elapsed_seconds", "failure")
                          if k in report}, ensure_ascii=False, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
