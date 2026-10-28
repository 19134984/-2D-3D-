"""Independent periodic-mode transport checks using the actual Entropic kernels.

Creates isolated test drivers under _local_checks/transport. Adapted from the independent review driver.
This transport test deliberately overrides the production geometry and forcing. The solver source is
read verbatim; only commondata test parameters and the periodic driver differ.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import os
import sys

ROOT = Path(__file__).resolve().parent
HERE = ROOT / "_local_checks" / "transport"
SOURCE = ROOT / "2DRBOpenaccEntropic.F90"
COMPILER = Path(sys.argv[1] if len(sys.argv) > 1 else r"C:/msys64/ucrt64/bin/gfortran.exe")
ENV = os.environ.copy()
ENV["PATH"] = str(COMPILER.parent) + os.pathsep + ENV.get("PATH", "")

DRIVER = r"""
program decay_check
  use commondata
  implicit none
  integer :: mode,i,j,a,n,ip,jp
  real(8) :: feq(0:8),geq(0:8),phase,amp0,amp1,wavenumber,expected,effective,amp
  allocate(f(nx,ny,0:8),g(nx,ny,0:8),f_post(0:nx+1,0:ny+1,0:8),g_post(0:nx+1,0:ny+1,0:8))
  allocate(rho(nx,ny),u(nx,ny),v(nx,ny),T(nx,ny),Fx(nx,ny),Fy(nx,ny),collisionStatus(nx,ny))
  omega=(/4d0/9d0,1d0/9d0,1d0/9d0,1d0/9d0,1d0/9d0,1d0/36d0,1d0/36d0,1d0/36d0,1d0/36d0/)
  omegaT=omega
  wavenumber=2d0*pi/dble(nx)
  do mode=1,2
    collisionStatus=0
    do j=1,ny
      do i=1,nx
        phase=wavenumber*(dble(i)-0.5d0)
        rho(i,j)=1d0; u(i,j)=0d0; v(i,j)=0d0; T(i,j)=Tref
        if(mode==1) v(i,j)=1d-5*sin(phase)
        if(mode==2) T(i,j)=Tref+1d-6*sin(phase)
        call entropic_equilibria(rho(i,j),u(i,j),v(i,j),T(i,j),feq,geq)
        f(i,j,:)=feq; g(i,j,:)=geq
      enddo
    enddo
    amp0=0d0; amp1=0d0
    do n=1,5000
      call collision()
      call collisionT()
      do j=1,ny
        do i=1,nx
          do a=0,8
            ip=modulo(i-1-ex(a),nx)+1
            jp=modulo(j-1-ey(a),ny)+1
            f(i,j,a)=f_post(ip,jp,a)
            g(i,j,a)=g_post(ip,jp,a)
          enddo
        enddo
      enddo
      call macro()
      if(any(collisionStatus/=0)) error stop 'Entropic collision rejected periodic mode'
      if(abs(sum(rho)-dble(nx*ny))>1d-8) error stop 'Periodic mode lost mass'
      if(abs(sum(rho*v))>1d-8) error stop 'Periodic mode lost momentum'
      if(mode==1.and.maxval(abs(v))>2d-5) error stop 'Equilibrium shear mode spuriously amplified'
      if(n==100.or.n==5000) then
        amp=0d0
        do j=1,ny
          do i=1,nx
            phase=wavenumber*(dble(i)-0.5d0)
            if(mode==1) amp=amp+v(i,j)*sin(phase)
            if(mode==2) amp=amp+(T(i,j)-Tref)*sin(phase)
          enddo
        enddo
        amp=2d0*amp/dble(nx*ny)
        if(n==100) amp0=amp
        if(n==5000) amp1=amp
      endif
    enddo
    expected=viscosity
    if(mode==2) expected=diffusivity
    effective=-log(amp1/amp0)/(wavenumber*wavenumber*4900d0)
    write(*,'(I2,5ES25.16E3)') mode,expected,effective,(effective/expected-1d0),amp0,amp1
  enddo
end program decay_check
"""


def replace_parameter(module, name, value):
    pattern = rf"(?im)^\s*real\(kind=8\),\s*parameter\s*::\s*{name}\s*=.*$"
    updated, count = re.subn(pattern, f"        real(kind=8), parameter :: {name}={value}", module)
    if count != 1:
        raise RuntimeError(f"Parameter {name}: expected one match, found {count}")
    return updated


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    text = SOURCE.read_text(encoding="utf-8-sig")
    module_end = re.search(r"(?im)^\s*end module commondata\s*$", text).end()
    module = text[:module_end]
    module = re.sub(r"(?im)^\s*integer\(kind=4\),\s*parameter\s*::\s*nx=\d+,\s*ny=\d+.*$",
                    "        integer(kind=4), parameter :: nx=64, ny=1", module)
    for name, value in {"tauf": "0.65d0", "gBeta": "0d0", "timeUnit": "1d0", "velocityUnit": "1d0"}.items():
        module = replace_parameter(module, name, value)
    names = ["collision", "collisionT", "macro", "entropic_pop_to_raw", "entropic_raw_to_pop",
             "entropic_moment_population", "entropic_equilibria", "entropic_relax"]
    kernels = []
    for name in names:
        match = re.search(rf"(?ims)^\s*subroutine {name}\(.*?^\s*end subroutine {name}\s*$", text)
        if not match:
            raise RuntimeError(f"Missing source kernel: {name}")
        kernels.append(match.group())
    results = []
    for pr in [0.7, 1.0, 7.0]:
        test_module = replace_parameter(module, "Prandtl", f"{pr}d0")
        outdir = HERE / ("entropic_decay_pr" + str(pr).replace(".", "_"))
        outdir.mkdir(exist_ok=True)
        source = outdir / "periodic_decay.F90"
        source.write_text(test_module + "\n" + DRIVER + "\n" + "\n".join(kernels), encoding="utf-8")
        exe = outdir / "periodic_decay.exe"
        compile_result = subprocess.run([str(COMPILER), "-cpp", "-ffree-line-length-none", "-O2", "-fcheck=all",
                                         "-fbacktrace", source.name, "-o", exe.name], cwd=outdir, env=ENV, capture_output=True, text=True)
        if compile_result.returncode:
            raise RuntimeError(compile_result.stdout + compile_result.stderr)
        run = subprocess.run([str(exe)], cwd=outdir, env=ENV, capture_output=True, text=True, timeout=120)
        if run.returncode:
            raise RuntimeError(run.stdout + run.stderr)
        rows = []
        for line in run.stdout.splitlines():
            fields = line.split()
            if len(fields) != 6:
                continue
            mode = int(fields[0])
            expected, effective, relative_error, amp_start, amp_end = map(float, fields[1:])
            rows.append(dict(mode="shear" if mode == 1 else "thermal", expected=expected,
                             effective=effective, relative_error=relative_error,
                             amplitude_at_step_100=amp_start, amplitude_at_step_5000=amp_end))
        if len(rows) != 2:
            raise RuntimeError(run.stdout + run.stderr)
        if any(abs(row['relative_error']) > 0.005 for row in rows):
            raise RuntimeError('Transport coefficient error exceeds 0.5 percent: ' + repr(rows))
        if any(row['amplitude_at_step_100'] <= 0 or row['amplitude_at_step_5000'] <= 0 for row in rows):
            raise RuntimeError('Fourier mode changed sign: ' + repr(rows))
        results.append(dict(Prandtl=pr, nx=64, ny=1, steps=5000, results=rows))
    output = ROOT / "transport_validation_local.json"
    payload = dict(status='passed', source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(), normalized_source_sha256=hashlib.sha256(text.encode('utf-8')).hexdigest(), method='actual Fortran kernels, GNU host, periodic Fourier decay; test overrides nx=64 ny=1, nu=0.05 and zero buoyancy; no wall or GPU validation', cases=results)
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
