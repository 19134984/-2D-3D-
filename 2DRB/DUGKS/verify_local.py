"""Local host-only verification. The production source is never modified.

Run with Python 3 and gfortran available. All generated sources, binaries, and
outputs go under the explicitly supplied ASCII build directory.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True)
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    root = Path(args.build_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    original = Path(__file__).with_name("2DRBOpenaccDUGKS.F90").read_text(encoding="utf-8")
    env = dict(os.environ)
    compiler = str(Path(args.compiler).resolve())
    env["PATH"] = str(Path(compiler).parent) + os.pathsep + env.get("PATH", "")
    results = {}

    def make_source(case="side", unsteady=False, periodic=False, conduction=False, steps=80, restart=False):
        s = original.replace("nx=256, ny=256", "nx=16, ny=16")
        s = s.replace("steadyCheckIntervalItc=2000", "steadyCheckIntervalItc=20")
        s = s.replace("outputSnapshotFile=1, outputPltFile=1", "outputSnapshotFile=0, outputPltFile=0")
        s = s.replace("outputSnapshotInterval=10.0d0", "outputSnapshotInterval=0.001d0")
        s = s.replace("outputSnapshotInterval=0.5d0", "outputSnapshotInterval=0.001d0")
        s = s.replace("itc_max=20000000", f"itc_max={steps}")
        s = s.replace("itc_max=ceiling(unsteadyRunDuration*timeUnit/dt)", f"itc_max={steps}")
        if case == "rb":
            s = s.replace("#define SideHeatedCell\n!#define RayleighBenardCell", "!#define SideHeatedCell\n#define RayleighBenardCell")
        if unsteady:
            s = s.replace("#define steadyFlow\n!#define unsteadyFlow", "!#define steadyFlow\n#define unsteadyFlow")
        if periodic:
            s = s.replace("!#define VerticalWallsPeriodical", "#define VerticalWallsPeriodical")
        if conduction:
            s = s.replace("gBeta=velocityUnit**2/(lengthUnit*temperatureDifference)", "gBeta=0.0d0")
            s = s.replace("rbInitialPerturbation=1.0d-3", "rbInitialPerturbation=0.0d0")
            # Exact transformed steady conductive state, g~=geq-(tau+dt/2)c.grad(T)w.
            gradient = "cx/(dble(nx)*dx)" if case == "side" else "cy/(dble(ny)*dy)"
            old = "g(i,j,k)=equilibrium(omega(k),cx,cy,T(i,j),0d0,0d0)"
            s = s.replace(old, old + f"+(tauc+0.5d0*dt)*omega(k)*temperatureDifference*({gradient})")
        if restart:
            s = s.replace("loadInitField=0", "loadInitField=1")
        return s

    def compile_source(directory, source, acc=False, acc_runtime=False):
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "case.F90").write_text(source, encoding="utf-8")
        flags = ["-cpp", "-ffree-line-length-none", "-O1", "-fcheck=all", "-fbacktrace"]
        if acc:
            flags += ["-fopenacc", "-fsyntax-only"]
        else:
            if acc_runtime:
                flags += ["-fopenacc"]
            flags += ["-o", "case.exe"]
        p = subprocess.run([compiler, *flags, "case.F90"], cwd=directory, env=env, text=True, capture_output=True)
        (directory / "compile.log").write_text(p.stdout + p.stderr, encoding="utf-8")
        if p.returncode:
            raise RuntimeError(f"Compile failed: {directory}\n{p.stderr}")

    def run_case(directory, expected_success=True):
        p = subprocess.run([str(directory / "case.exe")], cwd=directory, env=env, text=True, capture_output=True, timeout=120)
        (directory / "run.log").write_text(p.stdout + p.stderr, encoding="utf-8")
        if expected_success != (p.returncode == 0):
            raise RuntimeError(f"Unexpected return code: {directory}\n{p.stdout}\n{p.stderr}")
        return p

    def diagnostic(directory):
        lines = (directory / "FinalDiagnostics_2DOpenaccDUGKS.dat").read_text().splitlines()
        n = next(i for i, x in enumerate(lines) if "NuVolAvg ReVolAvg" in x)
        return [float(x) for x in lines[n + 1].split()]

    # Both directions and both flow modes; periodic RB is checked separately.
    for case in ("side", "rb"):
        for mode in (False, True):
            name = f"{case}_{'unsteady' if mode else 'steady'}"
            directory = root / name
            source = make_source(case=case, unsteady=mode)
            compile_source(directory, source)
            run_case(directory)
            values = diagnostic(directory)
            assert all(abs(x) < 1e8 for x in values), (name, values)
            results[name] = dict(NuVolAvg=values[0], ReVolAvg=values[1])
    directory = root / "rb_periodic"
    compile_source(directory, make_source(case="rb", unsteady=True, periodic=True))
    run_case(directory)
    results["rb_periodic"] = dict(zip(("NuVolAvg", "ReVolAvg"), diagnostic(directory)[:2]))

    for case in ("side", "rb"):
        directory = root / f"conduction_{case}"
        compile_source(directory, make_source(case=case, conduction=True))
        run_case(directory)
        values = diagnostic(directory)
        defect = max(abs(values[k] - 1.0) for k in (0, 2, 3, 4, 5))
        # Paper Eq.(44) reverses tangential nonequilibrium populations at an
        # adiabatic wall. Corner-adjacent fluxes therefore have a small defect
        # even for an initially linear conductive profile; retain that boundary.
        assert defect < 5e-5, (case, values)
        assert abs(values[1]) < 1e-8
        results[f"conduction_{case}"] = dict(max_Nu_error=defect, Re=values[1])

    # Stronger low-Mach buoyancy case: ensure wall forces do not leak mass.
    directory = root / "buoyant_mass"
    source = make_source(steps=600).replace("Mach=0.01d0", "Mach=0.05d0")
    compile_source(directory, source)
    run_case(directory)
    history = [list(map(float, x.split())) for x in (directory / "NuRe_2DOpenaccDUGKS.dat").read_text().splitlines() if not x.startswith("#")]
    mass_error = max(abs(x[8] / history[0][8] - 1.0) for x in history)
    assert mass_error < 1e-11, mass_error
    results["buoyant_mass"] = dict(max_relative_mass_error=mass_error, steps=600)

    # Exact same numerical state after an uninterrupted and a checkpointed run.
    directory = root / "restart"
    compile_source(directory, make_source(steps=37))
    run_case(directory)
    with (directory / "Convergence_2DOpenaccDUGKS.dat").open("a") as handle:
        handle.write("100 1.0 1e-5 1e-5\n")  # Uncommitted output tail from a later failed segment.
    compile_source(directory, make_source(steps=80, restart=True))
    run_case(directory)
    meta = (directory / "reloadFile2DOpenaccDUGKS-latest.meta").read_text().splitlines()
    restarted = (directory / meta[1].strip()).read_bytes()
    fresh_dir = root / "side_steady"
    meta = (fresh_dir / "reloadFile2DOpenaccDUGKS-latest.meta").read_text().splitlines()
    fresh = (fresh_dir / meta[1].strip()).read_bytes()
    # Header is 32-byte magic + five int32 + twenty float64; state is seven int32 and two float64.
    offset = 32 + 5 * 4 + 20 * 8 + 7 * 4 + 2 * 8
    length = 2 * 16 * 16 * 9 * 8
    assert restarted[offset:offset + length] == fresh[offset:offset + length], "Restart transformed populations differ"
    assert (directory / "NuRe_2DOpenaccDUGKS.dat").read_bytes() == (fresh_dir / "NuRe_2DOpenaccDUGKS.dat").read_bytes()
    assert (directory / "Convergence_2DOpenaccDUGKS.dat").read_bytes() == (fresh_dir / "Convergence_2DOpenaccDUGKS.dat").read_bytes()
    results["restart_37_plus_43"] = dict(transformed_populations_bitwise_equal=True, history_bitwise_equal=True, convergence_tail_trimmed=True)
    # Reject changed physical parameters rather than silently restarting.
    compile_source(directory, make_source(steps=100, restart=True).replace("Mach=0.01d0", "Mach=0.02d0"))
    p = run_case(directory, expected_success=False)
    assert "Incompatible checkpoint physical" in p.stderr
    results["restart_parameter_guard"] = "passed"

    for case, periodic in (("side", False), ("rb", False), ("rb", True)):
        for mode in (False, True):
            directory = root / f"acc_syntax_{case}_{int(periodic)}_{int(mode)}"
            compile_source(directory, make_source(case=case, periodic=periodic, unsteady=mode), acc=True)
    results["gnu_openacc_directive_syntax"] = "passed; device execution remains unverified"
    # Execute OpenACC data regions, asynchronous queues, and routine calls on
    # the GNU host device. This still does not verify an NVIDIA GPU runtime.
    env["ACC_DEVICE_TYPE"] = "host"
    for case, mode in (("side", False), ("rb", True)):
        directory = root / f"acc_host_{case}"
        compile_source(directory, make_source(case=case, unsteady=mode), acc_runtime=True)
        run_case(directory)
        reference = diagnostic(root / f"{case}_{'unsteady' if mode else 'steady'}")
        actual = diagnostic(directory)
        delta = max(abs(a - b) for a, b in zip(actual, reference))
        assert delta < 1e-10, (case, delta)
        results[f"gnu_openacc_host_{case}"] = dict(max_diagnostic_difference_from_serial=delta)
    version = subprocess.run([compiler, "--version"], capture_output=True, text=True, env=env, check=True).stdout.splitlines()[0]
    source_path = Path(__file__).with_name("2DRBOpenaccDUGKS.F90")
    baseline = source_path.parent.parent / "均匀网格" / "2DRBOpenmp.F90"
    evidence = dict(
        date="2026-10-28",
        production_source=source_path.name,
        production_source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        primary_baseline="均匀网格/2DRBOpenmp.F90",
        primary_baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
        compiler=version,
        compiler_path=compiler,
        host_flags=["-cpp", "-ffree-line-length-none", "-O1", "-fcheck=all", "-fbacktrace"],
        openacc_flags=["-fopenacc"],
        openacc_runtime_environment={"ACC_DEVICE_TYPE": "host"},
        mesh=[16,16],
        normal_smoke_steps=80,
        tests=results,
        limitations=[
            "No nvfortran compilation or NVIDIA GPU execution was performed.",
            "No server submission was made.",
            "No paper benchmark agreement, grid independence, or statistical convergence is claimed.",
            "Zero-force conduction retains the paper full bounce-back adiabatic wall: mean numerical flux Nu defect is about 3.6e-6 on this test.",
            "Only the explicitly listed generated small-grid variants were run; the production 256x256 default was not run."
        ],
    )
    report = root / "verification.json"
    report.write_text(json.dumps(evidence, indent=2), encoding="utf-8")
    print(report)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
