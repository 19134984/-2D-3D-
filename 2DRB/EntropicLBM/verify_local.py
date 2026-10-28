"""Compile and check the actual Fortran source locally; never submits a server job.

GNU -fopenacc runs on ACC_DEVICE_TYPE=host. This checks directive compilation and
host execution, not NVHPC code generation or GPU execution. Generated variants
and outputs stay in _local_checks; the production source is never edited.
"""
from pathlib import Path
import json
import hashlib
import os
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SOURCE = (ROOT / "2DRBOpenaccEntropic.F90").read_text(encoding="utf-8-sig")
SOURCE_SHA256 = hashlib.sha256((ROOT / "2DRBOpenaccEntropic.F90").read_bytes()).hexdigest()
COMPILER = Path(sys.argv[1] if len(sys.argv) > 1 else r"C:/msys64/ucrt64/bin/gfortran.exe")
BUILD = ROOT / "_local_checks"
ENV = os.environ.copy()
ENV["PATH"] = str(COMPILER.parent) + os.pathsep + ENV.get("PATH", "")
ENV["ACC_DEVICE_TYPE"] = "host"
RESULTS = []


def run(args, cwd, expected=0):
    result = subprocess.run(args, cwd=cwd, env=ENV, text=True,
                            encoding="utf-8", errors="replace", capture_output=True)
    if result.returncode != expected:
        raise RuntimeError(f"{args[0]} returned {result.returncode}\n{result.stdout}\n{result.stderr}")
    return result.stdout + result.stderr


def compile_source(text, folder, name="solver", acc=True):
    folder.mkdir(parents=True, exist_ok=True)
    filename = f"{name}.F90"
    (folder / filename).write_text(text, encoding="utf-8")
    flags = [str(COMPILER), "-cpp", "-ffree-line-length-none", "-O1", "-g", "-fcheck=all", "-fbacktrace"]
    if acc:
        flags += ["-fopenacc"]
    run(flags + [filename, "-o", f"{name}.exe"], folder)
    return folder / f"{name}.exe"


def variant(case, *, flow=None, pr=0.71, duration=0.1, restart=False, steady_steps=160):
    text = SOURCE.replace("nx=256, ny=256", "nx=16, ny=16")
    text = text.replace("Rayleigh=1.0d6", "Rayleigh=1.0d4")
    text = text.replace("Prandtl=0.71d0", f"Prandtl={pr:.16f}d0")
    text = text.replace("steadyCheckIntervalItc=2000", "steadyCheckIntervalItc=20")
    text = text.replace("itc_max=20000000", f"itc_max={steady_steps}")
    text = text.replace("loadInitField=0", "loadInitField=1" if restart else "loadInitField=0")
    if (flow or ("unsteady" if case == "rb" else "steady")) == "unsteady":
        text = text.replace("#define steadyFlow\n!#define unsteadyFlow", "!#define steadyFlow\n#define unsteadyFlow")
    if case == "rb":
        text = text.replace("#define SideHeatedCell\n!#define RayleighBenardCell", "!#define SideHeatedCell\n#define RayleighBenardCell")
    text = text.replace("unsteadyRunDuration=1000.0d0", f"unsteadyRunDuration={duration:.16f}d0")
    text = text.replace("outputSnapshotInterval=0.5d0", "outputSnapshotInterval=0.02d0")
    text = text.replace("reloadFileInterval=100.0d0", "reloadFileInterval=0.04d0")
    text = text.replace("outputPltFileInterval=100.0d0", "outputPltFileInterval=0.05d0")
    return text


DRIVER = r"""
program contracts
  use commondata
  implicit none
  real(kind=8) :: p(0:8),r(0:8),back(0:8),feq(0:8),geq(0:8),eqraw(0:8)
  real(kind=8) :: dsraw(0:8),ds(0:8),dh(0:8),fstar(0:8),post(0:8),delta(0:8)
  real(kind=8) :: rho0,ux,uy,temp0,energy0,ay,expected,gradient,hpost,htest
  integer :: k,i,j,statusLocal
  rho0=1.01d0; ux=0.015d0; uy=-0.008d0; temp0=Tcold+0.007d0
  omega=(/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0, &
      1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
  omegaT=omega
  call entropic_equilibria(rho0,ux,uy,temp0,feq,geq)
  call entropic_pop_to_raw(feq,eqraw)
  if(maxval(abs(eqraw(0:2)-(/rho0,rho0*ux,rho0*uy/))).gt.2.0d-15) error stop 'f moments'
  call entropic_pop_to_raw(geq,r)
  energy0=rho0*(2.0d0*temp0+ux*ux+uy*uy)
  if(abs(r(0)-energy0).gt.2.0d-15) error stop 'g total energy'
  if(abs(r(1)-(energy0+2.0d0*rho0*cs2)*ux).gt.2.0d-15) error stop 'g heat flux'
  if(abs(r(3)-energy0*(cs2+ux*ux)-2.0d0*rho0*cs2*(cs2+2.0d0*ux*ux)).gt.2.0d-15) &
      error stop 'g contracted fourth moment'
  r=eqraw
  r(3)=r(3)+0.0002d0; r(4)=r(4)-0.0001d0; r(5)=r(5)+0.00015d0
  r(6)=r(6)+0.00012d0; r(7)=r(7)-0.00009d0; r(8)=r(8)+0.00008d0
  call entropic_raw_to_pop(r,p)
  call entropic_pop_to_raw(p,back)
  if(maxval(abs(r-back)).gt.2.0d-15) error stop 'natural moment inverse'
  dsraw=0.0d0; dsraw(3:5)=r(3:5)-eqraw(3:5)
  call entropic_raw_to_pop(dsraw,ds)
  dsraw=0.0d0; dsraw(6:8)=r(6:8)-eqraw(6:8)
  call entropic_raw_to_pop(dsraw,dh)
  call entropic_pop_to_raw(ds,back)
  if(maxval(abs(back(6:8))).gt.2.0d-15) error stop 's is not natural second order'
  if(Prandtl.le.1.0d0) then
    fstar=feq
  else
    call entropic_moment_population(rho0,rho0*ux,rho0*uy,r(3),r(4),r(5),fstar)
  endif
  call entropic_relax(p,feq,ds,dh,fstar,post,statusLocal)
  if(statusLocal.ne.0) error stop 'admissible entropy state rejected'
  call entropic_pop_to_raw(post,back)
  if(maxval(abs(back(0:2)-eqraw(0:2))).gt.2.0d-15) error stop 'collision conservation'
  if(Prandtl.le.1.0d0) then
    expected=r(3)+2.0d0*beta*(eqraw(3)-r(3))
  else
    expected=r(3)+2.0d0*beta1*(eqraw(3)-r(3))
  endif
  if(abs(back(3)-expected).gt.2.0d-15) error stop 'Pr-dependent viscosity'
  gradient=sum((-beta*dh)*(log(post/omega)+1.0d0))
  if(abs(gradient).gt.2.0d-13) error stop 'entropy gradient'
  hpost=sum(post*log(post/omega))
  do k=-5,5
    delta=post+0.1d0*dble(k)*dh
    if(any(delta.le.0.0d0)) cycle
    htest=sum(delta*log(delta/omega))
    if(htest.lt.hpost-2.0d-15) error stop 'entropy minimum'
  enddo
  allocate(f(nx,ny,0:8),g(nx,ny,0:8),f_post(0:nx+1,0:ny+1,0:8),g_post(0:nx+1,0:ny+1,0:8))
  allocate(rho(nx,ny),u(nx,ny),v(nx,ny),T(nx,ny),Fx(nx,ny),Fy(nx,ny),collisionStatus(nx,ny))
  collisionStatus=0; rho=rho0; T=temp0; u=ux; v=uy
  do k=0,8
    f(:,:,k)=feq(k); g(:,:,k)=geq(k)
  enddo
  call collision()
  call collisionT()
  ay=gBeta*(temp0-Tref)
  call entropic_pop_to_raw(f_post(1,1,0:8),back)
  if(abs(back(2)-rho0*(uy+ay)).gt.2.0d-15) error stop 'EDM momentum'
  if(abs(sum(g_post(1,1,0:8))-energy0).gt.2.0d-15) error stop 'g collision conservation'
  f=f_post(1:nx,1:ny,0:8); g=g_post(1:nx,1:ny,0:8)
  call macro()
  expected=temp0-(2.0d0*uy*ay+ay*ay)/2.0d0
  if(abs(T(1,1)-expected).gt.2.0d-15) error stop 'f-only EDM total-energy convention'
  if(abs(v(1,1)-(uy+ay+0.5d0*gBeta*(expected-Tref))).gt.2.0d-15) error stop 'diagnostic velocity'
  if(any(collisionStatus.ne.0)) error stop 'collision failure'
  ! Near-equilibrium regression: numerical low-order residues must not be amplified by gamma.
  T=Tref
  do j=1,ny
    do i=1,nx
      uy=1.0d-5*sin(2.0d0*pi*dble(i-1)/dble(nx))
      call entropic_equilibria(rho0,0.0d0,uy,Tref,feq,geq)
      f(i,j,0:8)=feq; g(i,j,0:8)=geq
    enddo
  enddo
  do k=1,20
    call collision()
    call collisionT()
    f=f_post(1:nx,1:ny,0:8); g=g_post(1:nx,1:ny,0:8)
    call macro()
    if(abs(sum(f)-rho0*dble(nx*ny)).gt.2.0d-10) error stop 'near-equilibrium mass amplification'
    if(maxval(abs(v)).gt.1.0001d-5) error stop 'near-equilibrium velocity amplification'
  enddo
  print *, 'PASS contracts Pr=',Prandtl,'entropy gradient=',gradient
end program contracts
"""


def records(path):
    """GNU sequential unformatted records; compare f/g, not output counters."""
    import struct
    data = path.read_bytes()
    pos = 0
    chunks = []
    while pos < len(data):
        count, = struct.unpack_from("<i", data, pos)
        pos += 4
        chunks.append(data[pos:pos + count])
        pos += count
        tail, = struct.unpack_from("<i", data, pos)
        assert tail == count
        pos += 4
    return chunks


def latest_checkpoint(folder):
    meta = (folder / "reloadFile2DOpenaccEntropic-latest.meta").read_text()
    suffix = re.search(r"reloadFileName\s+(\d+)", meta).group(1)
    return folder / f"reloadFile2DOpenaccEntropic-{suffix}.bin"


def main():
    BUILD.mkdir(exist_ok=True)
    for pr in (0.71, 1.0, 7.0):
        text = variant("side", pr=pr)
        kernels = re.sub(r"(?ims)^\s*program main\b.*?^\s*end program main\s*$", "", text)
        folder = BUILD / f"contracts_pr{pr}"
        exe = compile_source(kernels + "\n" + DRIVER, folder, acc=False)
        output = run([str(exe)], folder)
        RESULTS.append({"check": f"moments_entropy_total_energy_Pr{pr}", "result": output.strip()})
    for case in ("side", "rb"):
        folder = BUILD / f"{case}_smoke"
        exe = compile_source(variant(case), folder)
        output = run([str(exe)], folder)
        (folder / "stdout.txt").write_text(output)
        RESULTS.append({"check": f"GNU_OpenACC_host_{case}_160_steps", "result": "PASS"})
    for case, flow in (("side", "unsteady"), ("rb", "steady")):
        compile_source(variant(case, flow=flow), BUILD / f"{case}_{flow}_compile")
        RESULTS.append({"check": f"GNU_OpenACC_{case}_{flow}_compile", "result": "PASS"})
    for case in ("side", "rb"):
        split = BUILD / f"{case}_smoke"
        full = BUILD / f"{case}_continuous"
        resumed = variant(case, duration=0.2, restart=True, steady_steps=320)
        exe = compile_source(resumed, split)
        output = run([str(exe)], split)
        (split / "restart_stdout.txt").write_text(output)
        exe = compile_source(variant(case, duration=0.2, steady_steps=320), full)
        run([str(exe)], full)
        assert records(latest_checkpoint(split))[1:] == records(latest_checkpoint(full))[1:]
        if case == "rb":
            for history in ("Nu_VolAvg_2DOpenaccEntropic.dat", "Re_VolAvg_2DOpenaccEntropic.dat"):
                assert (split / history).read_bytes() == (full / history).read_bytes()
        RESULTS.append({"check": f"{case}_restart_equals_continuous_320_steps", "result": "PASS bitwise f/g"})
    split = BUILD / "side_offphase_restart"
    exe = compile_source(variant("side", steady_steps=157), split)
    run([str(exe)], split)
    exe = compile_source(variant("side", restart=True, steady_steps=320), split)
    run([str(exe)], split)
    full = BUILD / "side_continuous"
    assert records(latest_checkpoint(split))[1:] == records(latest_checkpoint(full))[1:]
    history = "convergence2DOpenaccEntropic.plt"
    assert (split / history).read_bytes() == (full / history).read_bytes()
    RESULTS.append({"check": "steady_offphase_checkpoint_157_to_320", "result": "PASS bitwise fields/reference fields/errors and full convergence history"})
    partial = BUILD / "rb_incomplete_statistics"
    text = variant("rb", duration=0.2).replace(
        "unsteadyItcRemaining = max(0, itc_max - restartItcOffset)",
        "unsteadyItcRemaining = min(160, max(0, itc_max - restartItcOffset))")
    exe = compile_source(text, partial)
    run([str(exe)], partial)
    summary = (partial / "NuRe_TimeAverage_2DOpenaccEntropic.txt").read_text()
    assert "INCOMPLETE" in summary and "full-window means not emitted" in summary
    RESULTS.append({"check": "unsteady_partial_window", "result": "PASS INCOMPLETE label; no full-window mean"})
    assert SOURCE_SHA256 == hashlib.sha256((ROOT / "2DRBOpenaccEntropic.F90").read_bytes()).hexdigest()
    payload = {"source_sha256": SOURCE_SHA256,
               "compiler": run([str(COMPILER), "--version"], ROOT).splitlines()[0],
               "scope": "GNU OpenACC host; no NVHPC, GPU, convection benchmark or grid convergence validation",
               "checks": RESULTS}
    report = json.dumps(payload, indent=2, ensure_ascii=False)
    (BUILD / "results.json").write_text(report, encoding="utf-8")
    (ROOT / "validation_local.json").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
