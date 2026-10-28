"""多级Luo副本的本地回归；仅修改ASCII验证目录内的算例副本。"""
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import tempfile
import time

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SOURCE_PATH = REPO / "多块网格/2DRBOpenaccMultilevel.F90"
PARENT_PATH = REPO / "多块网格/2DRBOpenaccMultiblock.F90"
ROOT = Path(os.environ.get("MULTILEVEL_VALIDATION_ROOT", str(Path(tempfile.gettempdir()) / "2drb-multilevel-checks")))
REPORT_PATH = REPO / "运行脚本/run_cases/Multilevel_implementation_20261027/validation.json"
COMPILER = Path("C:/msys64/ucrt64/bin/gfortran.exe")
ENV = os.environ.copy()
ENV["PATH"] = str(COMPILER.parent) + os.pathsep + ENV["PATH"]
ENV["ACC_DEVICE_TYPE"] = "host"
ENV["OMP_NUM_THREADS"] = "1"
SOURCE = SOURCE_PATH.read_text(encoding="utf-8")
PARENT = PARENT_PATH.read_text(encoding="utf-8")
REPORT = {
    "compiler": str(COMPILER),
    "backend": "GCC OpenACC host; no NVIDIA GPU validation",
    "build_root": str(ROOT),
    "parent_sha256": hashlib.sha256(PARENT_PATH.read_bytes()).hexdigest(),
    "source_sha256": hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest(),
    "tests": {},
}


def save():
    REPORT_PATH.write_text(json.dumps(REPORT, ensure_ascii=False, indent=2), encoding="utf-8")


def assign(text, name, value):
    text, n = re.subn(r"\b" + name + r"=\d+", name + "=" + str(value), text, count=1)
    assert n, name
    return text


def variant(levels=3, n=128, fine=16, inner=8, steps=32, periodic=False, steady=False, side=False,
            conduction=False, restart=False, old=False):
    text = PARENT if old else SOURCE
    for name, value in [("nx", n), ("ny", n)]:
        text = assign(text, name, value)
    if old:
        for name, value in [("refineRatio", 1 if levels == 1 else 2),
                            ("fineLayerCellsLeft", fine), ("fineLayerCellsRight", fine + 1),
                            ("fineLayerCellsBottom", fine), ("fineLayerCellsTop", fine + 1)]:
            text = assign(text, name, value)
        root_step = "refineRatio"
    else:
        for name, value in [("nGridLevels", levels),
                            ("finestLayerCellsLeftRight", fine), ("finestLayerCellsBottomTop", fine),
                            ("innerLayerCellsLeftRight", inner), ("innerLayerCellsBottomTop", inner)]:
            text = assign(text, name, value)
        root_step = "coarsestStep"
    # 用细步数指定短算例的时间间隔，正式源文件的t_ff参数不变。
    for name, value in [("outputSnapshotInterval", 8), ("reloadFileInterval", 16), ("outputPltFileInterval", 16)]:
        text, count = re.subn(r"\b" + name + r"=[0-9.]+d0", name + f"={value}.0d0/timeUnit", text)
        assert count, name
    text = re.sub(r"\bunsteadyRunDuration=[0-9.]+d0", f"unsteadyRunDuration={steps}.0d0/timeUnit", text)
    text = text.replace(f"ceiling(2000.0d0/dble({root_step}))", f"ceiling(8.0d0/dble({root_step}))")
    text = text.replace(f"ceiling(20000000.0d0/dble({root_step}))", f"ceiling({steps}.0d0/dble({root_step}))")
    if steady:
        text = re.sub(r"^!#define steadyFlow$", "#define steadyFlow", text, flags=re.M)
        text = re.sub(r"^#define unsteadyFlow$", "!#define unsteadyFlow", text, flags=re.M)
    if side:
        for name in ["RayleighBenardCell", "HorizontalWallsConstT", "VerticalWallsAdiabatic"]:
            text = re.sub(r"^#define " + name + r"$", "!#define " + name, text, flags=re.M)
        for name in ["SideHeatedCell", "HorizontalWallsAdiabatic", "VerticalWallsConstT"]:
            text = re.sub(r"^!#define " + name + r"$", "#define " + name, text, flags=re.M)
    if periodic:
        for name in ["VerticalWallsNoslip", "VerticalWallsAdiabatic"]:
            text = re.sub(r"^#define " + name + r"$", "!#define " + name, text, flags=re.M)
        for name in ["VerticalWallsPeriodicalU", "VerticalWallsPeriodicalT"]:
            text = re.sub(r"^!#define " + name + r"$", "#define " + name, text, flags=re.M)
    if conduction:
        # 只在验证副本的流场碰撞中关闭浮力，保留原始物理换算和输出标度。
        if old:
            text = text.replace("        gBeta &", "        0.0d0 &")
            text = text.replace("Sq,gBeta &", "Sq,0.0d0 &")
        else:
            text = text.replace("SqLevel(l),gBeta &", "SqLevel(l),0.0d0 &")
            # 精确线性导热解含恒定非平衡通量；只用geq起步会有正常的热流建立过程。
            # g_a=omega_a*T-omega_a*dx/Qk*(e_a dot grad T)，检验所有接口能否保持该解。
            direction, length = ("ex(a)", "nx") if side else ("ey(a)", "ny")
            text = text.replace(
                "block(b)%g(i,j,a)=omegaT(a)*block(b)%T(i,j)",
                "block(b)%g(i,j,a)=omegaT(a)*block(b)%T(i,j)-omegaT(a)*block(b)%dx/"
                f"QkLevel(block(b)%level)*{direction}*(Tcold-Thot)/dble({length})"
            )
    if restart:
        text = assign(text, "loadInitField", 1)
    # 导出完整f/g与宏观场，与旧两级/均匀分支逐项比较。
    if old:
        text = text.replace("    deallocate(f_coarse", "    call dump_reference()\n    deallocate(f_coarse", 1)
        dump = ["\nsubroutine dump_reference()\nuse commondata\nimplicit none\ninteger :: k\n",
                "open(newunit=k,file='raw-state.bin',access='stream',form='unformatted',status='replace')\n"]
        dump.append("if(refineRatio>1) then\n")
        for region in ["left", "right", "bottom", "top"]:
            dump.append("write(k) " + ",".join(v + "_" + region for v in ["f", "g", "rho", "u", "v", "T", "Fx", "Fy"]) + "\n")
        dump.append("endif\nwrite(k) " + ",".join(v + "_coarse" for v in ["f", "g", "rho", "u", "v", "T", "Fx", "Fy"]))
        dump.append("\nclose(k)\nend subroutine dump_reference\n")
    else:
        text = text.replace("        call exit_data_2d_openacc()", "        call update_host_all(.true.)\n"
                            "        call dump_reference()\n        call exit_data_2d_openacc()", 1)
        dump = ["\nsubroutine dump_reference()\nuse commondata\nimplicit none\ninteger :: k,b\n",
                "open(newunit=k,file='raw-state.bin',access='stream',form='unformatted',status='replace')\n",
                "do b=1,nBlocks\nwrite(k) block(b)%f,block(b)%g,block(b)%rho,block(b)%u,block(b)%v,block(b)%T,block(b)%Fx,block(b)%Fy\n",
                "enddo\nclose(k)\nend subroutine dump_reference\n"]
    return text + "".join(dump)


def compile_case(folder, text):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "solver.F90").write_text(text, encoding="utf-8")
    args = [str(COMPILER), "-cpp", "-ffree-line-length-none", "-fopenacc", "-O0", "-g",
            "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "solver.F90"]
    syntax = subprocess.run(args + ["-fsyntax-only"], cwd=folder, env=ENV, capture_output=True)
    (folder / "syntax.log").write_bytes(syntax.stdout + syntax.stderr)
    assert syntax.returncode == 0, syntax.stderr.decode(errors="replace")
    result = subprocess.run(args + ["-o", "solver.exe"], cwd=folder, env=ENV, capture_output=True)
    (folder / "compile.log").write_bytes(result.stdout + result.stderr)
    assert result.returncode == 0, result.stderr.decode(errors="replace")


def execute(folder, expect=0):
    t = time.monotonic()
    result = subprocess.run([str(folder / "solver.exe")], cwd=folder, env=ENV,
                            capture_output=True, timeout=180)
    (folder / "run.log").write_bytes(result.stdout + result.stderr)
    output = (result.stdout + result.stderr).decode(errors="replace")
    if expect == 0:
        assert result.returncode == 0, str(folder) + "\n" + output[-5000:]
    else:
        assert result.returncode != 0, "Invalid input unexpectedly ran"
    return {"exit_code": result.returncode, "seconds": time.monotonic() - t, "last_output": output[-1200:]}


def run_case(label, text, expect=0):
    folder = ROOT / label
    compile_case(folder, text)
    result = execute(folder, expect)
    REPORT["tests"][label] = result
    save()
    print(label, "PASS", round(result["seconds"], 2), flush=True)
    return folder


def checkpoint(folder):
    name = (folder / "reloadFile2DOpenaccMultilevel-latest.meta").read_text().splitlines()[1]
    return folder / name


def inspect_checkpoint(path):
    raw = path.read_bytes()
    offset = 0

    def read(fmt):
        nonlocal offset
        result = struct.unpack_from("<" + fmt, raw, offset)
        offset += struct.calcsize("<" + fmt)
        return result

    assert read("16s")[0].rstrip() == b"ML2DRESTART0001"
    config = read("16i")
    physics = read("12d")
    total, samples = read("2i")
    read("2d")
    (nb,) = read("i")
    assert nb == 4 * (config[2] - 1) + 1
    shapes = []
    for _ in range(nb):
        level, ni, nj, gi, gj = read("5i")
        dx, xoff, yoff = read("3d")
        assert dx == 2 ** (level - 1)
        count = ni * nj * (23 if config[-1] == 1 else 20)
        data = np.frombuffer(raw, dtype="<f8", count=count, offset=offset)
        assert np.isfinite(data).all()
        offset += count * 8
        shapes.append([level, ni, nj, dx, xoff, yoff])
    history = []
    for fine_level in range(1, config[2]):
        nf, cycles = read("2i")
        assert cycles == total // 2 ** fine_level
        data = np.frombuffer(raw, dtype="<f8", count=20 * nf * 3, offset=offset)
        assert np.isfinite(data).all()
        offset += data.size * 8
        history.append({"fine_level": fine_level, "receivers": nf, "cycles": cycles})
    assert offset == len(raw)
    return {"steps": total, "samples": samples, "blocks": shapes, "interfaces": history,
            "finest_Luo_Qk": physics[8], "finest_Luo_Qnu": physics[9]}


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    # 源码里的块内数值算法和逐矩打包/还原必须与保护的原文件一致。
    for name in ["collision", "streaming", "bounceback", "macro", "collisionT", "streamingT",
                 "bouncebackT", "macroT", "pack_interface_values", "unpack_interface_values",
                 "coarse_source_stencil", "lagrange_weights"]:
        pattern = r"(?ims)^  subroutine " + name + r"\b.*?^  end subroutine " + name + r"\b[^\n]*"
        def numerical_lines(body):
            return "\n".join(line for line in body.splitlines()
                             if not line.lstrip().startswith("!") or line.lstrip().startswith("!$acc"))
        assert numerical_lines(re.search(pattern, SOURCE).group()) == numerical_lines(re.search(pattern, PARENT).group()), name
    REPORT["unchanged_numerical_routines"] = 12
    for levels, n, fine in [(1, 32, 8), (2, 64, 8), (3, 128, 16), (4, 256, 16)]:
        folder = run_case(f"levels{levels}", variant(levels=levels, n=n, fine=fine))
        REPORT["tests"][f"levels{levels}"]["checkpoint"] = inspect_checkpoint(checkpoint(folder))
        save()
    for label, periodic, steady, side in [
        ("periodic3", True, False, False), ("steady3", False, True, True),
        ("side3", False, False, True)
    ]:
        run_case(label, variant(periodic=periodic, steady=steady, side=side))
    for levels, n, fine in [(1, 32, 8), (2, 64, 8)]:
        old = run_case(f"old{levels}", variant(levels=levels, n=n, fine=fine, old=True))
        new = ROOT / f"levels{levels}"
        a = np.fromfile(old / "raw-state.bin", dtype=np.float64)
        b = np.fromfile(new / "raw-state.bin", dtype=np.float64)
        np.testing.assert_allclose(a, b, rtol=1e-11, atol=1e-13)
        REPORT["tests"][f"baseline_comparison{levels}"] = {
            "full_state_max_abs_difference": float(np.max(np.abs(a - b))),
            "bitwise_equal": bool(np.array_equal(a, b)),
        }
        # 原统计定义下的体积、壁面和中线Nu/Re也应一致。
        h_old = np.loadtxt(old / "NuRe_2DOpenaccMultiblock.dat")
        h_new = np.loadtxt(new / "NuRe_2DOpenaccMultiblock.dat")
        np.testing.assert_allclose(h_old, h_new, rtol=1e-10, atol=1e-11)
        save()
    for levels, n, fine in [(3, 128, 16), (4, 256, 16)]:
        folder = run_case(f"conduction{levels}", variant(levels=levels, n=n, fine=fine, conduction=True))
        values = np.loadtxt(folder / "FinalNuRe_2DOpenaccMultiblock.dat")
        assert abs(values[1] - 1) < 1e-10
        assert max(abs(values[3:6] - 1)) < 1e-8
        assert abs(values[6] - n * n) < 1e-8 * n * n
        REPORT["tests"][f"conduction{levels}"]["NuVol_hot_cold_middle"] = values[[1, 3, 4, 5]].tolist()
        save()
    for label, levels, n, fine, steady, periodic in [
        ("restart3", 3, 128, 16, False, False), ("restart4", 4, 256, 16, False, False),
        ("restart_steady", 3, 128, 16, True, False), ("restart_periodic", 3, 128, 16, False, True)
    ]:
        opts = dict(levels=levels, n=n, fine=fine, steady=steady, periodic=periodic, side=steady)
        reference = ROOT / ("steady3" if steady else "periodic3" if periodic else f"levels{levels}")
        folder = run_case(label, variant(steps=16, **opts))
        compile_case(folder, variant(steps=32, restart=True, **opts))
        result = execute(folder)
        assert checkpoint(folder).read_bytes() == checkpoint(reference).read_bytes(), label + " checkpoint differs"
        np.testing.assert_array_equal(np.loadtxt(folder / "NuRe_2DOpenaccMultiblock.dat"),
                                      np.loadtxt(reference / "NuRe_2DOpenaccMultiblock.dat"))
        REPORT["tests"][label]["restart"] = "bitwise identical checkpoint and equal history"
        REPORT["tests"][label]["continuation"] = result
        save()
        print(label, "RESTART PASS", flush=True)
    # 重启必须拒绝篡改的有效历史前缀，同时允许丢弃检查点之后的残缺尾行。
    reference = ROOT / "levels3"
    for label, corrupt, mismatch in [("restart_bad_prefix", True, False),
                                      ("restart_bad_tail", False, False),
                                      ("restart_bad_mesh", False, True)]:
        folder = ROOT / label
        folder.mkdir(exist_ok=True)
        for name in ["reloadFile2DOpenaccMultilevel-latest.meta", checkpoint(reference).name,
                     "NuRe_2DOpenaccMultiblock.dat", "OutputFiles_2DOpenaccMultiblock.dat"]:
            shutil.copy2(reference / name, folder / name)
        history = folder / "NuRe_2DOpenaccMultiblock.dat"
        content = history.read_text()
        if corrupt:
            content = content.replace("1 8 ", "99 8 ", 1)
        else:
            content += "discarded incomplete tail\\n"
        history.write_text(content)
        opts = dict(levels=3, n=128, fine=16, steps=40, restart=True)
        if mismatch:
            opts["levels"] = 2
        compile_case(folder, variant(**opts))
        result = execute(folder, expect=1 if corrupt or mismatch else 0)
        REPORT["tests"][label] = result
        if corrupt:
            assert "Invalid retained Nu/Re history row" in result["last_output"]
        if mismatch:
            assert "configuration mismatch" in result["last_output"]
        save()
        print(label, "PASS", flush=True)
    for label, opts, error in [
        ("invalid_division", dict(levels=3, n=127), "divisible"),
        ("invalid_thickness", dict(levels=3, fine=2), "thickness"),
        ("invalid_center", dict(levels=4, n=128, fine=16), "central"),
        ("invalid_levels", dict(levels=0), "positive")
    ]:
        folder = run_case(label, variant(**opts), expect=1)
        assert error.lower() in (folder / "run.log").read_text(errors="replace").lower()
    default = variant(levels=3, n=1024, fine=64, inner=32, steps=8)
    folder = run_case("default1024", default)
    REPORT["tests"]["default1024"]["checkpoint"] = inspect_checkpoint(checkpoint(folder))
    # 只改变级数即可使用默认厚度切到2级、4级，短算例验证其配置与执行路径。
    for levels in [2, 4]:
        folder = run_case(f"default1024_{levels}", variant(levels=levels, n=1024, fine=64, inner=32, steps=8))
        REPORT["tests"][f"default1024_{levels}"]["checkpoint"] = inspect_checkpoint(checkpoint(folder))
    assert hashlib.sha256(PARENT_PATH.read_bytes()).hexdigest() == REPORT["parent_sha256"]
    REPORT["status"] = "passed"
    save()
    print("ALL CHECKS PASSED", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        REPORT["status"] = "failed"
        REPORT["failure"] = str(exc)
        save()
        raise
