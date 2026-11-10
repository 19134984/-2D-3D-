"""独立验证混合粗细比后处理的旧格式兼容和新快照空间积分；不改算例输入文件。"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / "多块网格/2DRBOpenaccMultilevelMixedRatioPostprocess.F90"
OLD = REPO / "多块网格/2DRBOpenaccMultilevelPostprocess.F90"
COMPILER = Path("C:/msys64/ucrt64/bin/gfortran.exe")


def digest(folder):
    paths = [folder / "NuRe_2DOpenaccMultiblock.dat", folder / "OutputFiles_2DOpenaccMultiblock.dat"]
    for row in paths[1].read_text().splitlines():
        if row.strip() and not row.startswith("#"):
            path = folder / row.split()[1]
            if path.suffix == ".bin" and "Snapshot" in path.name:
                paths.append(path)
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, required=True, help="ASCII 可写验证目录")
    parser.add_argument("--legacy-root", type=Path, help="原 verify_multilevel.py 的验证目录")
    parser.add_argument("--case", type=Path, action="append", default=[], help="额外真实混合比例算例目录，可重复")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    args.build_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PATH"] = str(COMPILER.parent) + os.pathsep + env["PATH"]
    source = SOURCE.read_text(encoding="utf-8").replace("processSnapshots=.false.", "processSnapshots=.true.")
    report = {"source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(), "tests": {}}
    cases = []
    if args.legacy_root:
        cases += [args.legacy_root / label for label in ["levels1", "levels2", "levels3", "levels4", "old1", "old2",
                                                        "steady3", "steady_final36", "periodic3", "side3", "conduction3", "conduction4"]]
    cases += args.case
    executables = {}
    for mode in ["steady", "unsteady"]:
        folder = args.build_root / ("build_" + mode)
        folder.mkdir(exist_ok=True)
        text = source
        if mode == "steady":
            text = text.replace("!#define steadyFlow", "#define steadyFlow").replace("#define unsteadyFlow", "!#define unsteadyFlow")
        (folder / "post.F90").write_text(text, encoding="utf-8")
        base = [str(COMPILER), "-cpp", "-ffree-line-length-none", "-O0", "-g", "-fcheck=all",
                "-ffpe-trap=invalid,zero,overflow", "post.F90"]
        for extra in [["-fsyntax-only"], ["-o", "post.exe"]]:
            result = subprocess.run(base + extra, cwd=folder, env=env, capture_output=True)
            assert result.returncode == 0, result.stderr.decode(errors="replace")
        executables[mode] = folder / "post.exe"

    for original in cases:
        before = digest(original)
        solver_text = (original / "solver.F90").read_text(encoding="utf-8") if (original / "solver.F90").exists() else ""
        steady = bool(re.search(r"(?m)^#define steadyFlow\s*$", solver_text)) or original.name.startswith("steady")
        mode = "steady" if steady else "unsteady"
        for convert in [False, True] if args.legacy_root and original.parent == args.legacy_root and original.name not in ["old1", "old2"] else [False]:
            label = original.name + ("_v2" if convert else "")
            folder = args.build_root / label
            folder.mkdir(exist_ok=True)
            # 仅复制索引、历史与小网格快照；输入不被后处理或格式转换覆盖。
            for name in ["NuRe_2DOpenaccMultiblock.dat", "OutputFiles_2DOpenaccMultiblock.dat", "Convergence_2DOpenaccMultiblock.dat"]:
                if (original / name).exists():
                    shutil.copy2(original / name, folder / name)
            for name in before:
                if not name.endswith(".bin"):
                    continue
                data = (original / name).read_bytes()
                if convert:
                    assert data[:16].rstrip() == b"ML2DSNAPSHOT0001"
                    levels = struct.unpack_from("<i", data, 16 + 2 * 4)[0]
                    header = 16 + 16 * 4 + 12 * 8
                    ratios = [1] if levels == 1 else [2] * (levels - 1)
                    data = b"ML2DSNAPSHOT0002" + data[16:header] + struct.pack("<" + "i" * len(ratios), *ratios) + data[header:]
                (folder / name).write_bytes(data)
            history = np.atleast_2d(np.loadtxt(folder / "NuRe_2DOpenaccMultiblock.dat"))
            executable = executables[mode]
            if not steady:
                # 验证副本覆盖完整短测窗口，不改生产后处理默认的500--1000 TFF。
                text = re.sub(r"unsteadyAverageStartTf=500.0d0", "unsteadyAverageStartTf=" + format(float(history[0, 2]), ".17e").replace("e", "d"), source)
                text = re.sub(r"unsteadyAverageEndTf=1000.0d0", "unsteadyAverageEndTf=" + format(float(history[-1, 2]), ".17e").replace("e", "d"), text)
                (folder / "window_post.F90").write_text(text, encoding="utf-8")
                built = subprocess.run([str(COMPILER), "-cpp", "-ffree-line-length-none", "-O0", "-fcheck=all", "window_post.F90", "-o", "window_post.exe"],
                                       cwd=folder, env=env, capture_output=True)
                assert built.returncode == 0, built.stderr.decode(errors="replace")
                executable = folder / "window_post.exe"
            result = subprocess.run([str(executable)], cwd=folder, env=env, capture_output=True, timeout=60)
            (folder / "post.log").write_bytes(result.stdout + result.stderr)
            assert result.returncode == 0, label + "\n" + (result.stdout + result.stderr).decode(errors="replace")
            snapshots = np.atleast_2d(np.loadtxt(folder / "SnapshotNuRe_2DOpenaccMultiblock.dat"))
            matched = snapshots[np.isin(snapshots[:, 0], history[:, 1])]
            assert matched.shape[0] == history.shape[0], "Snapshot/history sample mismatch: " + label
            np.testing.assert_allclose(matched, history[:, 1:], rtol=2e-12, atol=2e-12)
            curve = np.atleast_2d(np.loadtxt(folder / "NuReCurve_2DOpenaccMultiblock.dat"))
            np.testing.assert_array_equal(curve, history)
            stats_error = None
            if not steady:
                stats = (folder / "NuReStatistics_2DOpenaccMultiblock.dat").read_text()
                assert "INCOMPLETE" not in stats
                means = np.array(re.search(r"(?m)^whole:\s*(.*)$", stats).group(1).split(), dtype=float)
                values = history[:, 3:8].copy()
                values[:, 1] **= 2
                expected = np.trapezoid(values, x=history[:, 2], axis=0) / (history[-1, 2] - history[0, 2])
                expected[1] = np.sqrt(expected[1])
                np.testing.assert_allclose(means, expected, rtol=2e-13, atol=2e-13)
                stats_error = float(np.max(np.abs(means - expected)))
            assert before == digest(original), "Source case modified: " + str(original)
            report["tests"][label] = {"snapshot_recomputation_max_abs_difference": float(np.max(np.abs(matched - history[:, 1:]))),
                                       "time_statistics_max_abs_difference": stats_error,
                                       "text_curve_equal": True,
                                       "source_case_unchanged": True}
            print(label, "PASS", flush=True)
    report["status"] = "passed"
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("ALL MIXED POSTPROCESS CHECKS PASSED", flush=True)


if __name__ == "__main__":
    main()
