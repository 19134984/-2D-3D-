"""核对多级求解器与原版文本后处理、配套快照后处理的兼容性。"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import re
import subprocess

import numpy as np

REPO = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ["MULTILEVEL_VALIDATION_ROOT"])
REPORT_PATH = REPO / "运行脚本/run_cases/Multilevel_implementation_20261027/postprocess_validation.json"
OLD_PATH = REPO / "多块网格/2DRBOpenaccMultiblockPostprocess.F90"
NEW_PATH = REPO / "多块网格/2DRBOpenaccMultilevelPostprocess.F90"
old = OLD_PATH.read_text(encoding="utf-8")
new = NEW_PATH.read_text(encoding="utf-8")
COMPILER = Path("C:/msys64/ucrt64/bin/gfortran.exe")
ENV = os.environ.copy()
ENV["PATH"] = str(COMPILER.parent) + os.pathsep + ENV["PATH"]
ENV["ACC_DEVICE_TYPE"] = "host"
ENV["OMP_NUM_THREADS"] = "1"
REPORT = {"old_postprocess_sha256": hashlib.sha256(OLD_PATH.read_bytes()).hexdigest(),
          "new_postprocess_sha256": hashlib.sha256(NEW_PATH.read_bytes()).hexdigest(), "tests": {}}


def save():
    REPORT_PATH.write_text(json.dumps(REPORT, ensure_ascii=False, indent=2), encoding="utf-8")


def post_source(text, steady, start, end, snapshots):
    text = text.replace("processSnapshots=.false.", "processSnapshots=" + (".true." if snapshots else ".false."))
    text = re.sub(r"unsteadyAverageStartTf=500.0d0",
                  "unsteadyAverageStartTf=" + format(start, ".17e").replace("e", "d"), text)
    text = re.sub(r"unsteadyAverageEndTf=1000.0d0",
                  "unsteadyAverageEndTf=" + format(end, ".17e").replace("e", "d"), text)
    if steady:
        text = text.replace("!#define steadyFlow", "#define steadyFlow")
        text = text.replace("#define unsteadyFlow", "!#define unsteadyFlow")
    return text


def run_post(folder, stem, source):
    (folder / (stem + ".F90")).write_text(source, encoding="utf-8")
    flags = [str(COMPILER), "-cpp", "-ffree-line-length-none", "-O0", "-g", "-fcheck=all",
             "-ffpe-trap=invalid,zero,overflow", stem + ".F90"]
    syntax = subprocess.run(flags + ["-fsyntax-only"], cwd=folder, env=ENV, capture_output=True)
    assert syntax.returncode == 0, syntax.stderr.decode(errors="replace")
    built = subprocess.run(flags + ["-o", stem + ".exe"], cwd=folder, env=ENV, capture_output=True)
    assert built.returncode == 0, built.stderr.decode(errors="replace")
    result = subprocess.run([str(folder / (stem + ".exe"))], cwd=folder, env=ENV,
                            capture_output=True, timeout=60)
    (folder / (stem + ".log")).write_bytes(result.stdout + result.stderr)
    assert result.returncode == 0, str(folder) + "\n" + (result.stdout + result.stderr).decode(errors="replace")


def digest(folder):
    files = [folder / "NuRe_2DOpenaccMultiblock.dat", folder / "OutputFiles_2DOpenaccMultiblock.dat"]
    for row in files[1].read_text().splitlines():
        if row.strip() and not row.startswith("#"):
            name = row.split()[1]
            if name.endswith(".bin"):
                files.append(folder / name)
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def main():
    # 添加一个终点不在采样时刻的稳态算例，核对最终定量结果来自真实终场。
    spec = importlib.util.spec_from_file_location("multilevel_verifier", REPO / "运行脚本/tools/verify_multilevel.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ROOT = ROOT
    folder = ROOT / "steady_final36"
    module.compile_case(folder, module.variant(steady=True, side=True, steps=36))
    module.execute(folder)
    final = (folder / "FinalNuRe_2DOpenaccMultiblock.dat").read_text()
    assert "# itc_total = 36" in final
    assert "# error_check_itc = 32" in final
    assert "# converged = " in final
    assert "# errorU errorT = " in final
    cases = ["levels1", "levels2", "levels3", "levels4", "periodic3", "steady3", "side3",
             "conduction3", "conduction4", "old1", "old2", "steady_final36"]
    for label in cases:
        folder = ROOT / label
        steady = label.startswith("steady")
        history = np.atleast_2d(np.loadtxt(folder / "NuRe_2DOpenaccMultiblock.dat"))
        assert history.shape[1] == 14
        start, end = float(history[0, 2]), float(history[-1, 2])
        before = digest(folder)
        # 原后处理的文件名、读取参数和历史解析保持原样；仅在验证副本中缩短平均窗口。
        run_post(folder, "old_post", post_source(old, steady, start, end, False))
        old_curve = (folder / "NuReCurve_2DOpenaccMultiblock.dat").read_bytes()
        old_stats = None
        if steady:
            old_convergence = (folder / "ConvergenceCurve_2DOpenaccMultiblock.dat").read_bytes()
        else:
            old_stats = (folder / "NuReStatistics_2DOpenaccMultiblock.dat").read_bytes()
            assert b"INCOMPLETE" not in old_stats
        # 配套副本读取原始块快照，重新独立计算全部空间指标。
        run_post(folder, "new_post", post_source(new, steady, start, end, True))
        assert old_curve == (folder / "NuReCurve_2DOpenaccMultiblock.dat").read_bytes()
        if steady:
            assert old_convergence == (folder / "ConvergenceCurve_2DOpenaccMultiblock.dat").read_bytes()
        else:
            assert old_stats == (folder / "NuReStatistics_2DOpenaccMultiblock.dat").read_bytes()
            means = re.search(r"(?m)^whole:\s*(.*)$", old_stats.decode()).group(1).split()
            expected = np.trapezoid(history[:, 3:8].copy() * np.array([1, 0, 1, 1, 1]),
                                    x=history[:, 2], axis=0) / (end - start)
            expected[1] = np.sqrt(np.trapezoid(history[:, 4] ** 2, x=history[:, 2]) / (end - start))
            np.testing.assert_allclose(np.array(means, dtype=float), expected, rtol=2e-13, atol=2e-13)
        snapshots = np.atleast_2d(np.loadtxt(folder / "SnapshotNuRe_2DOpenaccMultiblock.dat"))
        assert snapshots.shape[1] == 13
        matched = snapshots[np.isin(snapshots[:, 0], history[:, 1])]
        assert matched.shape[0] == history.shape[0]
        np.testing.assert_allclose(matched, history[:, 1:], rtol=2e-12, atol=2e-12)
        if label == "steady_final36":
            last = snapshots[snapshots[:, 0] == 36]
            assert last.shape[0] == 1
            np.testing.assert_allclose(last[0, 1:], np.loadtxt(folder / "FinalNuRe_2DOpenaccMultiblock.dat"),
                                       rtol=2e-12, atol=2e-12)
        assert before == digest(folder), "Postprocessor modified solver data"
        tecplot = sorted(folder.glob("buoyancyCavity2DOpenaccMultiblockTecplot-*.dat"))[-1]
        with tecplot.open() as file:
            assert file.readline().strip() == 'VARIABLES="x/L","y/L","u","v","T","rho"'
            assert "SOLUTIONTIME=" in file.readline()
            assert len(file.readline().split()) == 6
        REPORT["tests"][label] = {
            "original_text_postprocessor": "passed without filename/format changes",
            "snapshot_recomputation_max_abs_difference": float(np.max(np.abs(matched - history[:, 1:]))),
            "statistics_and_curves_equal": True,
            "solver_data_unchanged": True,
            "tecplot_six_columns": True,
        }
        save()
        print(label, "POSTPROCESS PASS", flush=True)
    assert hashlib.sha256(OLD_PATH.read_bytes()).hexdigest() == "3c0c5b6b598cfc5f5f5e44dc9d85c198dd139babd4cc3a4f6df11e2c307d3a86"
    REPORT["status"] = "passed"
    save()
    print("ALL POSTPROCESS CHECKS PASSED", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        REPORT["status"] = "failed"
        REPORT["failure"] = str(exc)
        save()
        raise
