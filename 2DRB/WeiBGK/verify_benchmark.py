"""小规模侧热方腔物理核验；仅修改隔离副本，不代表高Ra或网格无关验证。"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile

from verify_local import compiler_environment, replace_parameter, run_command, select_macros

HERE = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", type=Path)
    parser.add_argument("--source", type=Path, default=HERE / "2DRBOpenaccWeiBGK.F90")
    parser.add_argument("--grid", type=int, default=32)
    parser.add_argument("--steps", type=int, default=100000)
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--report", type=Path, default=HERE / "convection_validation_local.json")
    args = parser.parse_args()
    compiler = args.compiler or Path(shutil.which("gfortran") or "gfortran")
    if not compiler.is_file():
        parser.error("Supply an existing --compiler")
    original = args.source.read_bytes()
    source = select_macros(original.decode("utf-8"), steady=True, rb=False)
    for name, value in (("nx", str(args.grid)), ("ny", str(args.grid)), ("Rayleigh", "1.0d3"),
                        ("Prandtl", "0.71d0"), ("Mach", "0.05d0"), ("itc_max", str(args.steps)),
                        ("steadyCheckIntervalItc", "500"), ("loadInitField", "0"),
                        ("outputSnapshotFile", "0"), ("outputPltFile", "0"), ("outputReloadFile", "0")):
        source = replace_parameter(source, name, value)
    work = Path(tempfile.mkdtemp(prefix="wei-bgk-convection-")).resolve()
    (work / "solver.F90").write_text(source)
    env = compiler_environment(compiler)
    report = {"passed": False, "source_sha256": hashlib.sha256(original).hexdigest(),
              "verified_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "case": "2D closed side-heated square cavity, steady", "Ra": 1000, "Pr": 0.71,
              "Ma": 0.05, "grid": [args.grid, args.grid], "maximum_steps": args.steps,
              "backend": "GNU Fortran OpenACC host; not NVIDIA GPU", "build_directory": str(work),
              "reference": {"paper": "Wang et al. (2013), DOI 10.1016/j.camwa.2012.07.001, Table 2/9",
                            "Nu": 1.1178},
              "limitations": "One low-Ra grid; no grid independence, high-Ra RB or GPU validation"}
    try:
        flags = ["-cpp", "-fopenacc", "-O2", "-ffree-line-length-none", "-fcheck=all",
                 "-ffpe-trap=invalid,zero,overflow", "-fbacktrace"]
        run_command([str(compiler), *flags, "-fsyntax-only", "solver.F90"], work, env, args.timeout)
        run_command([str(compiler), *flags, "solver.F90", "-o", "solver.exe"], work, env, args.timeout)
        output = run_command([str(work / "solver.exe")], work, env, args.timeout)
        (work / "run.log").write_text(output.stdout + output.stderr)
        settings = (work / "SimulationSettings2DOpenaccWeiBGK.txt").read_text()
        metrics = {}
        for key in ("Nu_global", "Nu_hot", "Nu_cold"):
            match = re.findall(rf"{key}\s*=\s*([+\-0-9.EeDd]+)", settings)
            if not match:
                raise AssertionError(f"Missing {key}")
            metrics[key] = float(match[-1].replace("D", "E"))
        lines = [line.split() for line in output.stdout.splitlines()
                 if re.fullmatch(r"\s*\d+\s+[+\-0-9.EeDd]+\s+[+\-0-9.EeDd]+\s*", line)]
        if not lines:
            raise AssertionError("Missing steady error history")
        last = lines[-1]
        metrics.update(final_check_step=int(last[0]), errorU=float(last[1]), errorT=float(last[2]))
        metrics["Nu_relative_reference_error"] = abs(metrics["Nu_global"] / 1.1178 - 1)
        metrics["maximum_wall_volume_Nu_difference"] = max(abs(metrics[k]-metrics["Nu_global"])
                                                            for k in ("Nu_hot", "Nu_cold"))
        report["metrics"] = metrics
        report["checks"] = {"steady_thresholds_met": metrics["errorU"] <= 1e-7 and metrics["errorT"] <= 1e-7,
                            "Nu_within_two_percent": metrics["Nu_relative_reference_error"] < 0.02,
                            "wall_volume_Nu_difference_below_0p02": metrics["maximum_wall_volume_Nu_difference"] < 0.02}
        report["passed"] = all(report["checks"].values())
    except Exception as error:
        report["failure"] = str(error)
    report["source_unchanged"] = args.source.read_bytes() == original
    report["passed"] = report["passed"] and report["source_unchanged"]
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
