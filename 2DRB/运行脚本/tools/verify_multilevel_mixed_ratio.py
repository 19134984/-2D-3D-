"""混合逐接口比例的独立本地回归；仅改写 ASCII 验证目录中的算例副本。"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import re
import shutil
import struct
import tempfile

import numpy as np

REPO = Path(__file__).resolve().parents[2]
BASE_PATH = REPO / "多块网格/2DRBOpenaccMultilevel.F90"
SOURCE_PATH = REPO / "多块网格/2DRBOpenaccMultilevelMixedRatio.F90"
REPORT_PATH = REPO / "运行脚本/run_cases/Multilevel_mixed_ratio_implementation_20261030/validation.json"
ROOT = Path(os.environ.get("MULTILEVEL_MIXED_VALIDATION_ROOT",
                          str(Path(tempfile.gettempdir()) / "2drb-mixed-ratio-checks")))
BASE = BASE_PATH.read_text(encoding="utf-8")
SOURCE = SOURCE_PATH.read_text(encoding="utf-8")
BASE_HASH = hashlib.sha256(BASE_PATH.read_bytes()).hexdigest()

# 复用已有 GCC OpenACC host 编译/运行基础设施，不修改原验证脚本。
spec = importlib.util.spec_from_file_location("multilevel_validation_helpers", Path(__file__).with_name("verify_multilevel.py"))
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)
REPORT = {
    "compiler": str(h.COMPILER),
    "backend": "GCC OpenACC host; not NVIDIA GPU validation",
    "build_root": str(ROOT),
    "baseline_path": str(BASE_PATH),
    "baseline_sha256": BASE_HASH,
    "source_sha256": hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest(),
    "tests": {},
}
h.ROOT = ROOT
h.REPORT_PATH = REPORT_PATH
h.REPORT = REPORT
REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)


def save():
    REPORT_PATH.write_text(json.dumps(REPORT, ensure_ascii=False, indent=2), encoding="utf-8")


def variant(ratios=(4, 2), n=256, fine=32, inner=16, steps=32,
            restart=False, conduction=False, baseline=False, periodic=False, steady=False):
    h.SOURCE = BASE if baseline else SOURCE
    text = h.variant(levels=len(ratios) + 1, n=n, fine=fine, inner=inner,
                     steps=steps, restart=restart, conduction=conduction,
                     periodic=periodic, steady=steady, side=steady)
    # 两个主文件的正式默认Ra不同；回归比较必须使用完全相同的物理参数。
    text, count = re.subn(r"\bRayleigh=[0-9.]+d[0-9]+", "Rayleigh=1.0d7", text, count=1)
    assert count == 1
    if not baseline:
        # NVHPC 24.3 对零长 parameter 数组的 product 会崩溃；单级使用[1]占位。
        literal = "[" + ",".join(map(str, ratios)) + "]" if ratios else "[1]"
        text, count = re.subn(r"(interfaceRefineRatio\(max\(1,nGridLevels-1\)\)\s*=)\[[^\n]*?\]",
                              lambda m: m[1] + literal, text, count=1)
        assert count == 1, "Could not replace interfaceRefineRatio parameter"
    return text


def checkpoint(folder):
    candidates = sorted(folder.glob("reloadFile2DOpenaccMultilevel*-latest.meta"))
    assert len(candidates) == 1, candidates
    return folder / candidates[0].read_text().splitlines()[1]


def inspect_checkpoint(path, ratios):
    raw = path.read_bytes()
    offset = 0

    def read(fmt):
        nonlocal offset
        answer = struct.unpack_from("<" + fmt, raw, offset)
        offset += struct.calcsize("<" + fmt)
        return answer

    assert read("16s")[0].rstrip() == b"ML2DRESTART0002"
    config = read("16i")
    physics = read("12d")
    assert config[2] == len(ratios) + 1
    stored_ratios = tuple(ratios) if ratios else (1,)
    assert read(f"{max(1, len(ratios))}i") == stored_ratios
    total, samples = read("2i")
    read("2d")
    nb = read("i")[0]
    assert nb == 4 * len(ratios) + 1
    cumulative = [1] + np.cumprod(ratios, dtype=int).tolist()
    shapes = []
    for _ in range(nb):
        level, ni, nj, gi, gj = read("5i")
        dx, xoff, yoff = read("3d")
        assert dx == cumulative[level - 1], (level, dx, cumulative)
        count = ni * nj * (23 if config[-1] == 1 else 20)
        data = np.frombuffer(raw, dtype="<f8", count=count, offset=offset)
        assert np.isfinite(data).all()
        offset += count * 8
        shapes.append([level, ni, nj, gi, gj, dx, xoff, yoff])
    history = []
    for index in range(len(ratios)):
        nf, cycles = read("2i")
        assert cycles == total // cumulative[index + 1]
        data = np.frombuffer(raw, dtype="<f8", count=20 * nf * 3, offset=offset)
        assert np.isfinite(data).all()
        offset += data.size * 8
        history.append({"ratio": ratios[index], "fine_level": index + 1,
                        "receivers": nf, "coarse_cycles": cycles})
    assert offset == len(raw), (offset, len(raw))
    return {"steps": total, "samples": samples, "ratios": list(ratios),
            "stored_ratio_array": list(stored_ratios),
            "cumulative_grid_steps": cumulative, "blocks": shapes, "interfaces": history,
            "finest_Luo_Qk": physics[8], "finest_Luo_Qnu": physics[9]}


def valid_case(label, ratios=(4, 2), **opts):
    folder = h.run_case(label, variant(ratios=ratios, **opts))
    REPORT["tests"][label]["checkpoint"] = inspect_checkpoint(checkpoint(folder), ratios)
    field = np.fromfile(folder / "raw-state.bin", dtype=np.float64)
    assert np.isfinite(field).all()
    history = np.atleast_2d(np.loadtxt(folder / "NuRe_2DOpenaccMultiblock.dat"))
    assert np.isfinite(history).all()
    save()
    return folder


def restart_case(label, ratios=(4, 2), reference=None, **opts):
    folder = valid_case(label, ratios=ratios, steps=16, **opts)
    h.compile_case(folder, variant(ratios=ratios, steps=32, restart=True, **opts))
    result = h.execute(folder)
    assert checkpoint(folder).read_bytes() == checkpoint(reference).read_bytes(), label
    assert (folder / "raw-state.bin").read_bytes() == (reference / "raw-state.bin").read_bytes()
    np.testing.assert_array_equal(np.loadtxt(folder / "NuRe_2DOpenaccMultiblock.dat"),
                                  np.loadtxt(reference / "NuRe_2DOpenaccMultiblock.dat"))
    REPORT["tests"][label]["restart"] = "bitwise identical checkpoint, full field, and history"
    REPORT["tests"][label]["continuation"] = result
    save()
    print(label, "RESTART PASS", flush=True)
    return folder


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    # 块内算法和原始逐矩非平衡转换不得因混合比例而改变。
    numerical_routines = ["collision", "streaming", "bounceback", "macro", "collisionT", "streamingT",
                          "bouncebackT", "macroT", "pack_interface_values", "unpack_interface_values",
                          "coarse_source_stencil", "lagrange_weights"]
    for name in numerical_routines:
        pattern = r"(?ims)^  subroutine " + name + r"\b.*?^  end subroutine " + name + r"\b[^\n]*"

        def numerical(body):
            return "\n".join(line for line in body.splitlines()
                             if not line.lstrip().startswith("!") or line.lstrip().startswith("!$acc"))

        assert numerical(re.search(pattern, SOURCE).group()) == numerical(re.search(pattern, BASE).group()), name
    REPORT["unchanged_numerical_routines"] = numerical_routines

    baseline = h.run_case("baseline22", variant(ratios=(2, 2), n=128, fine=16, inner=8, baseline=True))
    equal = valid_case("ratios22", ratios=(2, 2), n=128, fine=16, inner=8)
    a = np.fromfile(baseline / "raw-state.bin", dtype=np.float64)
    b = np.fromfile(equal / "raw-state.bin", dtype=np.float64)
    np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-13)
    np.testing.assert_array_equal(np.loadtxt(baseline / "NuRe_2DOpenaccMultiblock.dat"),
                                  np.loadtxt(equal / "NuRe_2DOpenaccMultiblock.dat"))
    REPORT["tests"]["ratios22"]["baseline_full_state_max_abs_difference"] = float(np.max(np.abs(a - b)))
    REPORT["tests"]["ratios22"]["baseline_bitwise_equal"] = bool(np.array_equal(a, b))
    save()

    mixed42 = valid_case("ratios42")
    mixed24 = valid_case("ratios24", ratios=(2, 4))
    valid_case("one_level", ratios=(), n=32, fine=8, inner=8)
    valid_case("four_levels", ratios=(2, 2, 2), n=256, fine=16, inner=8)
    for label, opts in [("conduction42", {}), ("conduction24", {"ratios": (2, 4)})]:
        folder = valid_case(label, conduction=True, **opts)
        values = np.loadtxt(folder / "FinalNuRe_2DOpenaccMultiblock.dat")
        assert abs(values[1] - 1) < 1e-10
        assert max(abs(values[3:6] - 1)) < 1e-8
        assert abs(values[6] - 256 * 256) < 1e-8 * 256 * 256
        REPORT["tests"][label]["NuVol_hot_cold_middle"] = values[[1, 3, 4, 5]].tolist()
        save()
    restart_case("restart42", reference=mixed42)
    restart_case("restart24", ratios=(2, 4), reference=mixed24)

    # [4,2] 和 [2,4] 的总比均为8，仍必须检查每个接口而非只检查累计总比。
    folder = ROOT / "restart_ratio_mismatch"
    folder.mkdir(exist_ok=True)
    for path in mixed42.iterdir():
        if path.suffix in (".meta", ".bin", ".dat"):
            shutil.copy2(path, folder / path.name)
    h.compile_case(folder, variant(ratios=(2, 4), steps=40, restart=True))
    result = h.execute(folder, expect=1)
    assert "configuration mismatch" in result["last_output"].lower(), result["last_output"]
    REPORT["tests"]["restart_ratio_mismatch"] = result
    save()
    print("restart_ratio_mismatch PASS", flush=True)

    for label, opts, message in [
        ("invalid_domain_division", {"n": 255}, "divisible"),
        ("invalid_inner_division", {"fine": 31}, "divisible"),
        ("invalid_outer_thickness", {"fine": 8}, "thickness"),
        ("invalid_inner_thickness", {"ratios": (2, 4), "inner": 8}, "thickness"),
        ("invalid_ratio", {"ratios": (3, 2), "n": 252, "fine": 36}, "ratio"),
    ]:
        folder = h.run_case(label, variant(**opts), expect=1)
        assert message in (folder / "run.log").read_text(errors="replace").lower()

    assert hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_HASH, "Protected baseline changed"
    REPORT["status"] = "passed"
    save()
    print("ALL MIXED-RATIO CHECKS PASSED", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        REPORT["status"] = "failed"
        REPORT["failure"] = str(exc)
        save()
        raise
