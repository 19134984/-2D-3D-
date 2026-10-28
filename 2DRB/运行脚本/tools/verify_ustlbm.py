"""USTLBM 本地验证：只修改隔离副本，不修改正式 Fortran 源文件。

GCC OpenACC host 会真正解析并执行 OpenACC 指令；这些检查不能证明
NVHPC/NVIDIA GPU 的设备捕获、异步同步或性能正确。未连接或提交服务器。

示例：
  python 运行脚本/tools/verify_ustlbm.py --stages syntax,smoke,conduction,analytic,restart,postprocess
  python 运行脚本/tools/verify_ustlbm.py --stages benchmark --benchmark-grids 32,64

所有构建、运行、输入和日志保存在独立 ASCII 临时目录。JSON 报告记录正式源
及每个验证副本的 SHA-256；分阶段调用保留同一正式源的其他阶段结果。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile
import time
import traceback

import numpy as np


REPO = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = REPO / "STLBM/2DRBOpenaccUSTLBM.F90"
DEFAULT_REPORT = REPO / "运行脚本/run_cases/USTLBM_implementation_20261028/validation.json"
STAGES = ("syntax", "smoke", "conduction", "analytic", "benchmark", "restart", "postprocess")
FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
BENCHMARK_REFERENCE = {
    "author": "G. de Vahl Davis", "year": 1983,
    "title": "Natural convection of air in a square cavity: A bench mark numerical solution",
    "journal": "International Journal for Numerical Methods in Fluids 3:249-264",
    "doi": "10.1002/fld.1650030305",
    "url": "https://onlinelibrary.wiley.com/doi/10.1002/fld.1650030305",
    "role": "independent classical benchmark for the same differentially heated square cavity, Pr=0.71",
    "metric_note": "continuum mean Nu reference; finite-grid Nu_hot/Nu_cold/Nu_global/Nu_middle are reported separately",
    "provenance_note": "these low-Ra reference values are not claimed to be copied item by item from the 2017 USTLBM paper",
}
CONFIGS = {
    "steady_sideheated": (True, True, False),
    "unsteady_RB_adiabatic": (False, False, False),
    "unsteady_RB_periodic": (False, False, True),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fnum(value: float) -> str:
    return f"{value:.16e}".replace("e", "d")


def set_parameter(text: str, name: str, value: str, *, required: bool = True) -> str:
    """仅替换 commondata 内 parameter 声明；不碰运行时赋值和诊断字符串。"""
    match = re.search(r"(?ims)^[ \t]*module[ \t]+commondata\b.*?^[ \t]*end[ \t]+module\b[^\n]*", text)
    if not match:
        raise ValueError("commondata module was not found")
    lines = match.group(0).splitlines(keepends=True)
    output = []
    count = 0
    index = 0
    while index < len(lines):
        block = [lines[index]]
        index += 1
        if not re.match(r"[ \t]*(?:real|integer|logical)\s*\([^)]*\)\s*,\s*parameter\s*::", block[0], re.I):
            output.extend(block)
            continue
        while block[-1].split("!", 1)[0].rstrip().endswith("&") and index < len(lines):
            block.append(lines[index])
            index += 1
        statement = " ".join(line.split("!", 1)[0].strip().strip("&").strip() for line in block)
        prefix, declarations = statement.split("::", 1)
        parts, begin, depth = [], 0, 0
        for pos, char in enumerate(declarations):
            depth += (char in "([") - (char in ")]")
            if char == "," and depth == 0:
                parts.append(declarations[begin:pos].strip())
                begin = pos + 1
        parts.append(declarations[begin:].strip())
        found = False
        for pos, part in enumerate(parts):
            if re.match(re.escape(name) + r"\s*=", part, re.I):
                parts[pos] = name + "=" + value
                count += 1
                found = True
        if found:
            indent = re.match(r"[ \t]*", block[0]).group(0)
            output.append(indent + prefix.strip() + " :: " + ", ".join(parts) + "\n")
        else:
            output.extend(block)
    module = "".join(output)
    if required and count == 0:
        raise ValueError(f"parameter declaration not found: {name}")
    return text[:match.start()] + module + text[match.end():]


def set_macro(text: str, name: str, enabled: bool) -> str:
    pattern = r"(?m)^[ \t]*!?[ \t]*#define[ \t]+" + re.escape(name) + r"\b[^\n]*"
    text, count = re.subn(pattern, ("" if enabled else "!") + "#define " + name, text)
    if count != 1:
        raise ValueError(f"expected one macro declaration for {name}, found {count}")
    return text


def configure(source: str, config: str, *, n: int = 24, rayleigh: float = 1e3,
              steps: int = 64, restart: bool = False, conduction: bool = False,
              benchmark: bool = False) -> str:
    steady, side, periodic = CONFIGS[config]
    selections = {
        "steadyFlow": steady, "unsteadyFlow": not steady,
        "SideHeatedCell": side, "RayleighBenardCell": not side,
        "HorizontalWallsNoslip": True, "VerticalWallsNoslip": not periodic,
        "VerticalWallsPeriodicalU": periodic,
        "HorizontalWallsConstT": not side, "VerticalWallsConstT": side,
        "HorizontalWallsAdiabatic": side,
        "VerticalWallsAdiabatic": not side and not periodic,
        "VerticalWallsPeriodicalT": periodic,
    }
    text = source
    for name, enabled in selections.items():
        text = set_macro(text, name, enabled)
    replacements = {
        "nx": str(n), "ny": str(n), "Rayleigh": fnum(rayleigh),
        "Prandtl": "0.71d0", "Mach": "0.05d0",
        "loadInitField": "1" if restart else "0",
        "epsU": "1.0d-8", "epsT": "1.0d-8",
        "outputSnapshotInterval": "8.0d0/timeUnit",
        "reloadFileInterval": "32.0d0/timeUnit",
        "outputPltFileInterval": "32.0d0/timeUnit",
        "unsteadyRunDuration": f"{steps}.0d0/timeUnit",
        "outputSnapshotFile": "1", "outputReloadFile": "1", "outputPltFile": "1",
    }
    # 非稳态 itc_max 继续从目标时长推导，避免破坏重启的绝对时间含义。
    if steady:
        replacements["itc_max"] = str(steps)
    for name, value in replacements.items():
        text = set_parameter(text, name, value)
    if conduction:
        text = set_parameter(text, "gBeta", "0.0d0")
    if benchmark:
        # 不为长稳态算例存储每8步一幅快照；最终输出路径仍由正式主程序执行。
        for name in ("outputSnapshotInterval", "reloadFileInterval", "outputPltFileInterval"):
            text = set_parameter(text, name, "2000.0d0/timeUnit")
        text = set_parameter(text, "outputSnapshotFile", "0")
        text = set_parameter(text, "outputPltFile", "0")
        # 最多100个history点，避免只改采样频率而溢出原稳态历史缓存。
        text = set_parameter(text, "dimensionlessTimeMax", str(max(256, steps // 2000 + 8)))
    return text


def replace_main(text: str, program: str) -> str:
    pattern = r"(?ims)^[ \t]*program[ \t]+main\b.*?^[ \t]*end[ \t]+program(?:[ \t]+main)?[^\n]*"
    text, count = re.subn(pattern, lambda _: program.strip() + "\n", text)
    if count != 1:
        raise ValueError(f"expected one main program, found {count}")
    return text


def conduction_driver(steps: int) -> str:
    # xp/yp由正式initial产生；该程序只提供精确初值并调用实际ustlbm_step。
    return f"""
program main
    use openacc
    use commondata
    implicit none
    integer :: i,j,k
    call acc_init(acc_device_host)
    call initial()
    rho=1.0d0
    u=0.0d0
    v=0.0d0
    do j=1,ny
        do i=1,nx
#ifdef SideHeatedCell
            T(i,j)=Thot+(Tcold-Thot)*(xp(i)-xp(0))/(xp(nx+1)-xp(0))
#else
            T(i,j)=Thot+(Tcold-Thot)*(yp(j)-yp(0))/(yp(ny+1)-yp(0))
#endif
        enddo
    enddo
    call enter_data_2d_openacc()
    do k=1,{steps}
        call ustlbm_step()
    enddo
    !$acc update self(u,v,T,rho) async(1)
    !$acc wait(1)
    open(99,file='validation_fields.bin',form='unformatted',access='stream',status='replace')
    write(99) nx,ny,{steps}
    write(99) timeUnit,diffusivity,Thot,Tcold
    write(99) xp,yp
    write(99) u,v,T,rho
    close(99)
    call exit_data_2d_openacc()
end program main
"""


def periodic_ghosts(text: str) -> str:
    """只在解析测试副本替换两个物理边界routine；推进公式保持原样。"""
    for suffix in ("old", "star"):
        name = "build_ustlbm_" + suffix + "_ghosts_acc"
        dest = [field + "_" + suffix for field in ("rho", "u", "v", "T")]
        src = ["rho", "u", "v", "T"] if suffix == "old" else dest
        arrays = ",".join(dict.fromkeys(src + dest))
        assignments = "\n".join(f"                {d}(i,j)={s}(ic,jc)" for d, s in zip(dest, src))
        condition = ".true." if suffix == "old" else "i.EQ.0.OR.i.EQ.nx+1.OR.j.EQ.0.OR.j.EQ.ny+1"
        replacement = f"""subroutine {name}()
    use commondata
    implicit none
    integer :: i,j,ic,jc
    !$acc parallel loop collapse(2) default(none) present({arrays}) async(1) private(ic,jc)
    do j=0,ny+1
        do i=0,nx+1
            if({condition}) then
                ic=modulo(i-1,nx)+1
                jc=modulo(j-1,ny)+1
{assignments}
            endif
        enddo
    enddo
end subroutine {name}
"""
        pattern = r"(?ims)^[ \t]*subroutine[ \t]+" + name + r"\([^)]*\).*?^[ \t]*end[ \t]+subroutine[ \t]+" + name + r"[^\n]*"
        text, count = re.subn(pattern, lambda _: replacement, text)
        if count != 1:
            raise ValueError(f"expected one {name}, found {count}")
    return text


def periodic_driver(steps: int, mx: int = 1, my: int = 1) -> str:
    program = conduction_driver(steps)
    pattern = r"(?s)#ifdef SideHeatedCell\n.*?#endif"
    initial = f"            T(i,j)=0.1d0*sin(2.0d0*pi*{mx}*(dble(i)-0.5d0)/dble(nx)) &\n" \
              f"                *cos(2.0d0*pi*{my}*(dble(j)-0.5d0)/dble(ny))"
    return re.sub(pattern, lambda _: initial, program, count=1)


def field_dump(directory: Path) -> dict:
    raw = (directory / "validation_fields.bin").read_bytes()
    nx, ny, steps = struct.unpack_from("<3i", raw)
    timeunit, diffusivity, hot, cold = struct.unpack_from("<4d", raw, 12)
    values = np.frombuffer(raw, dtype="<f8", offset=44)
    if values.size != nx + ny + 4 + 4 * nx * ny or not np.isfinite(values).all():
        raise AssertionError("invalid field dump")
    fields = values[nx+ny+4:].reshape(4, -1)
    return {"nx": nx, "ny": ny, "steps": steps, "timeunit": timeunit,
            "diffusivity": diffusivity, "hot": hot, "cold": cold,
            "xp": values[:nx+2], "yp": values[nx+2:nx+ny+4],
            "fields": dict(zip(("u", "v", "T", "rho"),
                               (a.reshape((nx, ny), order="F") for a in fields)))}


def read_records(path: Path) -> list[bytes]:
    """读取GCC默认4字节record marker；损坏和多余数据均报错。"""
    data = path.read_bytes()
    offset = 0
    records = []
    while offset < len(data):
        if offset + 8 > len(data):
            raise ValueError(f"truncated record marker: {path}")
        size, = struct.unpack_from("<i", data, offset)
        if size < 0 or offset + size + 8 > len(data):
            raise ValueError(f"invalid record size {size}: {path}")
        tail, = struct.unpack_from("<i", data, offset + 4 + size)
        if size != tail:
            raise ValueError(f"record markers disagree: {path}")
        records.append(data[offset + 4:offset + 4 + size])
        offset += size + 8
    return records


def checkpoint(path: Path, n: int) -> dict:
    records = read_records(path)
    if len(records) != 5 or len(records[0]) != 24:
        raise ValueError(f"expected one signature/dimensions record and four state records: {path}")
    signature = records[0][:16]
    if b"USTLBM" not in signature:
        raise ValueError(f"invalid USTLBM signature: {signature!r}")
    if struct.unpack("<2i", records[0][16:]) != (n, n):
        raise ValueError(f"checkpoint dimensions differ: {path}")
    fields = {}
    for name, raw in zip(("rho", "u", "v", "T"), records[1:]):
        values = np.frombuffer(raw, dtype="<f8")
        if values.size != n * n or not np.isfinite(values).all():
            raise ValueError(f"nonfinite/wrong-size {name} record: {path}")
        fields[name] = values.reshape((n, n), order="F")
    return {"signature": signature.decode("ascii").strip(), "fields": fields}


def metadata(directory: Path) -> tuple[Path, dict[str, str]]:
    candidates = list(directory.glob("*latest.meta"))
    if len(candidates) != 1:
        raise ValueError(f"expected one latest metadata file in {directory}")
    meta = {}
    for line in candidates[0].read_text(encoding="utf-8").splitlines():
        row = line.split(maxsplit=1)
        if len(row) == 2:
            meta[row[0]] = row[1].strip()
    return candidates[0], meta


def latest_checkpoint(directory: Path, n: int) -> tuple[Path, dict, dict]:
    _, meta = metadata(directory)
    if "reloadFileName" not in meta or "itc_total" not in meta:
        raise ValueError("reload metadata has no filename or cumulative step")
    name = meta["reloadFileName"].strip("'\"")
    path = directory / name
    if not path.is_file():
        # 兼容旧元数据仅写编号的布局，但不会猜测多个文件。
        candidates = list(directory.glob(f"reload*-{name}.bin"))
        if len(candidates) != 1:
            raise ValueError(f"metadata target does not exist: {name}")
        path = candidates[0]
    return path, checkpoint(path, n), meta


def histories(directory: Path) -> dict[str, np.ndarray]:
    result = {}
    for prefix in ("Nu_VolAvg_", "Re_VolAvg_"):
        paths = list(directory.glob(prefix + "*.dat"))
        if len(paths) != 1:
            raise ValueError(f"expected one {prefix} history in {directory}")
        data = np.atleast_2d(np.loadtxt(paths[0]))
        if data.shape[1] != 2 or not np.isfinite(data).all():
            raise ValueError(f"invalid history: {paths[0]}")
        if data.shape[0] > 1 and not np.all(np.diff(data[:, 0]) > 0):
            raise ValueError(f"nonincreasing/duplicate history time: {paths[0]}")
        result[prefix] = data
    np.testing.assert_array_equal(result["Nu_VolAvg_"][:, 0], result["Re_VolAvg_"][:, 0])
    return result


def diagnostics(directory: Path) -> dict[str, float]:
    log = (directory / "run.log").read_text(encoding="utf-8", errors="replace")
    settings = list(directory.glob("SimulationSettings*.txt"))
    if settings:
        log += "\n" + settings[0].read_text(encoding="utf-8", errors="replace")
    values = {}
    for name in ("Nu_global", "Nu_hot", "Nu_cold", "Nu_middle", "u_mid_max", "v_mid_max", "NuVolAvg"):
        matches = re.findall(re.escape(name) + r"\s*=\s*(" + FLOAT_RE + ")", log)
        if matches:
            values[name] = float(matches[-1].replace("D", "E").replace("d", "e"))
    return values


class Validation:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.source_path = args.source.resolve()
        self.source_bytes = self.source_path.read_bytes()
        self.source = self.source_bytes.decode("utf-8").replace("\r\n", "\n")
        self.source_hash = hashlib.sha256(self.source_bytes).hexdigest()
        self.root = Path(tempfile.mkdtemp(prefix="ustlbm-verify-", dir=args.work_root))
        if not str(self.root).isascii():
            raise ValueError(f"build directory is not ASCII: {self.root}")
        (self.root / "formal_source.F90").write_bytes(self.source_bytes)
        self.env = os.environ.copy()
        self.env["PATH"] = str(Path(args.compiler).parent) + os.pathsep + self.env.get("PATH", "")
        self.env["ACC_DEVICE_TYPE"] = "host"
        self.env["OMP_NUM_THREADS"] = "1"
        self.flags = [args.compiler, "-cpp", "-fopenacc", "-ffree-line-length-none"]
        self.report_path = args.report.resolve()
        self.report = {}
        if self.report_path.is_file():
            old = json.loads(self.report_path.read_text(encoding="utf-8"))
            if old.get("source_sha256") == self.source_hash:
                self.report = old
            elif old.get("source_sha256"):
                archive = self.report_path.with_name(self.report_path.stem + "." + old["source_sha256"][:12] + ".json")
                if not archive.exists():
                    archive.write_bytes(self.report_path.read_bytes())
        version = subprocess.run([args.compiler, "--version"], env=self.env, capture_output=True, text=True)
        self.report.update({
            "source": str(self.source_path), "source_sha256": self.source_hash,
            "compiler": args.compiler, "compiler_version": version.stdout.splitlines()[0],
            "backend": "GCC OpenACC host; NVIDIA GPU and NVHPC behavior remain unverified",
            "build_root": str(self.root), "updated_at_local": time.strftime("%Y-%m-%d %H:%M:%S"),
            "benchmark_reference": BENCHMARK_REFERENCE,
        })
        self.report.setdefault("tests", {})

    def save(self):
        current_hash = sha256(self.source_path)
        # 一个阶段始终从self.source的不可变快照构建。协作方可能同时修改正式源；
        # 保留旧快照证据并显式标注失配，不能把旧结果冒称为当前正式源验证。
        self.report["current_source_sha256_at_save"] = current_hash
        self.report["matches_current_source"] = current_hash == self.source_hash
        self.report_path.parent.mkdir(parents=True, exist_ok=True)
        self.report_path.write_text(json.dumps(self.report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    def command(self, command: list[str], directory: Path, logname: str, timeout: int | None = None):
        with (directory / logname).open("wb") as log:
            try:
                proc = subprocess.run(command, cwd=directory, env=self.env, stdout=log,
                                      stderr=subprocess.STDOUT, timeout=timeout or self.args.timeout)
            except subprocess.TimeoutExpired:
                raise RuntimeError(f"timeout; see {directory / logname}") from None
        if proc.returncode:
            tail = (directory / logname).read_text(encoding="utf-8", errors="replace")[-3500:]
            raise RuntimeError(f"exit={proc.returncode}; {directory / logname}\n{tail}")

    def build(self, label: str, text: str, *, syntax: bool = False, directory: Path | None = None) -> Path:
        directory = directory or self.root / label
        directory.mkdir(parents=True, exist_ok=True)
        src = directory / "solver.F90"
        src.write_text(text, encoding="utf-8", newline="\n")
        flags = self.flags + (["-fsyntax-only"] if syntax else ["-O2", "-fcheck=bounds", "-ffpe-trap=invalid,zero,overflow"])
        cmd = flags + [str(src)] + ([] if syntax else ["-o", str(directory / "solver.exe")])
        self.command(cmd, directory, "compile.log")
        (directory / "source.sha256").write_text(sha256(src) + "\n", encoding="ascii")
        return directory

    def run(self, directory: Path):
        self.command([str(directory / "solver.exe")], directory, "run.log")
        log = (directory / "run.log").read_text(encoding="utf-8", errors="replace")
        if re.search(r"(?im)\b(?:nan|infinity)\b|^[ \t]*error:", log):
            raise RuntimeError(f"invalid runtime diagnostic; see {directory / 'run.log'}")

    def syntax(self):
        result = {}
        for label in CONFIGS:
            directory = self.build("syntax_" + label, configure(self.source, label), syntax=True)
            result[label] = {"status": "PASS", "directory": str(directory), "source_sha256": sha256(directory / "solver.F90")}
        return result

    def smoke(self):
        result = {}
        for label in CONFIGS:
            directory = self.build("smoke_" + label, configure(self.source, label))
            self.run(directory)
            path, state, meta = latest_checkpoint(directory, 24)
            total = int(meta["itc_total"])
            if not 64 <= total <= 65:
                raise AssertionError(f"unexpected completed horizon: {total}")
            hist = histories(directory)
            snapshots = sorted(directory.glob("*Snapshot-*.bin"))
            if not snapshots:
                raise AssertionError("no postprocessing snapshot")
            records = read_records(snapshots[-1])
            if len(records) != 4 or any(len(raw) != 24 * 24 * 8 for raw in records):
                raise AssertionError("snapshot must preserve four u/v/T/rho records")
            for raw in records:
                if not np.isfinite(np.frombuffer(raw, dtype="<f8")).all():
                    raise AssertionError("nonfinite snapshot")
            result[label] = {"status": "PASS", "directory": str(directory), "completed_steps": total,
                             "checkpoint_sha256": sha256(path), "signature": state["signature"],
                             "history_rows": int(hist["Nu_VolAvg_"].shape[0]),
                             "snapshot_records": "u_scaled,v_scaled,T,rho; four n*n float64 records"}
        return result

    def conduction(self):
        result = {}
        for label in CONFIGS:
            text = configure(self.source, label, conduction=True)
            text = replace_main(text, conduction_driver(64))
            directory = self.build("conduction_" + label, text)
            self.run(directory)
            dump = field_dump(directory)
            nx, ny = dump["nx"], dump["ny"]
            xp, yp = dump["xp"], dump["yp"]
            hot, cold = dump["hot"], dump["cold"]
            u, v, temp, rho = [dump["fields"][name] for name in ("u", "v", "T", "rho")]
            side = CONFIGS[label][1]
            coord = xp if side else yp
            line = hot + (cold-hot)*(coord[1:-1]-coord[0])/(coord[-1]-coord[0])
            exact = np.broadcast_to(line[:, None] if side else line[None, :], (nx, ny))
            terr = float(np.max(np.abs(temp-exact)))
            speed = float(max(np.max(np.abs(u)), np.max(np.abs(v))))
            density = float(np.max(np.abs(rho-1)))
            tline = temp.mean(axis=1 if side else 0)
            slope = float(np.polyfit(coord[1:-1], tline, 1)[0])
            nu = -slope*(coord[-1]-coord[0])/(hot-cold)
            if terr > 1e-10 or speed > 1e-12 or density > 1e-12 or abs(nu-1) > 1e-9:
                raise AssertionError(f"conduction: T={terr}, speed={speed}, rho={density}, Nu={nu}")
            result[label] = {"status": "PASS", "directory": str(directory), "steps": dump["steps"],
                             "T_max_abs_error": terr, "velocity_max_abs": speed,
                             "rho_max_abs_error": density, "independent_Nu_from_temperature_slope": nu,
                             "diffusivity": dump["diffusivity"], "timeUnit": dump["timeunit"]}
        return result

    def analytic(self):
        # 连续解和单步傅里叶乘子由论文权重独立给出，不另写推进算法作参考。
        # a0=0,a_cardinal=1/6,a_diagonal=1/12 ->
        # E=(cos(theta_x)+cos(theta_y)+cos(theta_x)*cos(theta_y))/3,
        # G=E+(tauT-1)*(E**2-1), tauT=0.5+1.5*kappa.
        result = {}
        errors = []
        cases = [("single_step", 24, 1, 2, 3)] + [
            (f"continuous_N{n}", n, 64*(n//32)**2, 1, 1) for n in (32, 64, 128)]
        for label, n, steps, mx, my in cases:
            text = configure(self.source, "unsteady_RB_periodic", n=n, steps=max(8, steps), conduction=True)
            text = set_parameter(text, "diffusivity", "0.1d0")
            text = periodic_ghosts(text)
            text = replace_main(text, periodic_driver(steps, mx, my))
            directory = self.build("analytic_" + label, text)
            self.run(directory)
            dump = field_dump(directory)
            field = dump["fields"]["T"]
            x = (np.arange(n)+0.5)/n
            initial = 0.1*np.sin(2*np.pi*mx*x)[:, None]*np.cos(2*np.pi*my*x)[None, :]
            tx, ty = 2*np.pi*mx/n, 2*np.pi*my/n
            e = (np.cos(tx)+np.cos(ty)+np.cos(tx)*np.cos(ty))/3
            tau = 0.5+1.5*dump["diffusivity"]
            g = e+(tau-1)*(e*e-1)
            discrete = initial*g**steps
            exact = initial*np.exp(-dump["diffusivity"]*(tx*tx+ty*ty)*steps)
            discrete_error = float(np.max(np.abs(field-discrete)))
            continuous_l2 = float(np.sqrt(np.mean((field-exact)**2))/np.sqrt(np.mean(exact**2)))
            if discrete_error > 3e-12:
                raise AssertionError(f"Fourier multiplier/index/sign mismatch: {label}: {discrete_error}")
            for name in ("u", "v"):
                if np.max(np.abs(dump["fields"][name])) > 1e-12:
                    raise AssertionError("zero-force analytic diffusion induced velocity")
            if np.max(np.abs(dump["fields"]["rho"]-1)) > 1e-12:
                raise AssertionError("analytic diffusion changed constant density")
            result[label] = {"status": "PASS", "directory": str(directory), "grid": n,
                             "steps": steps, "Fourier_E": float(e), "Fourier_G": float(g),
                             "max_error_against_independent_discrete_symbol": discrete_error,
                             "relative_L2_error_against_continuous_diffusion": continuous_l2,
                             "diffusion_time_kappa_t_over_L2": dump["diffusivity"]*steps/n**2}
            if label != "single_step":
                errors.append(continuous_l2)
        orders = [float(np.log2(errors[i]/errors[i+1])) for i in range(2)]
        result["continuous_refinement"] = {
            "status": "PASS" if min(orders) > 1.8 else "FAIL", "observed_orders": orders,
            "scope": "continuous transient diffusion, joint spatial/time refinement under diffusive scaling; not isolated temporal order"}
        return result

    def restart(self):
        result = {}
        for label in ("unsteady_RB_adiabatic", "unsteady_RB_periodic"):
            full = self.build("restart_full_" + label, configure(self.source, label, steps=64))
            self.run(full)
            fullpath, fullstate, fullmeta = latest_checkpoint(full, 24)
            for split_step in (32, 35):
                key = label + f"_checkpoint{split_step}"
                split = self.build("restart_split_" + key, configure(self.source, label, steps=split_step))
                self.run(split)
                _, _, firstmeta = latest_checkpoint(split, 24)
                if int(firstmeta["itc_total"]) != split_step:
                    raise AssertionError(f"split run did not stop at checkpoint step {split_step}")
                original_histories = {}
                if split_step == 35:
                    for prefix in ("Nu_VolAvg_", "Re_VolAvg_"):
                        paths = list(split.glob(prefix + "*.dat"))
                        if len(paths) != 1:
                            raise AssertionError("missing history before damaged-tail restart test")
                        with paths[0].open("ab") as handle:
                            handle.write(b"incomplete discarded tail\n")
                        original_histories[paths[0].name] = paths[0].read_bytes()
                self.build("unused", configure(self.source, label, steps=64, restart=True), directory=split)
                self.run(split)
                splitpath, splitstate, splitmeta = latest_checkpoint(split, 24)
                if int(fullmeta["itc_total"]) != 64 or int(splitmeta["itc_total"]) != 64:
                    raise AssertionError("restart did not reach the same absolute horizon")
                diffs = {}
                for name in ("rho", "u", "v", "T"):
                    diffs[name] = float(np.max(np.abs(fullstate["fields"][name] - splitstate["fields"][name])))
                    np.testing.assert_array_equal(fullstate["fields"][name], splitstate["fields"][name])
                for name, array in histories(full).items():
                    np.testing.assert_array_equal(array, histories(split)[name])
                backups = {}
                for filename, original in original_histories.items():
                    matches = [path for path in split.iterdir()
                               if path.is_file() and path.name != filename and path.read_bytes() == original]
                    if not matches:
                        raise AssertionError(f"no byte-exact backup of pre-trim damaged history: {filename}")
                    backups[filename] = [path.name for path in matches]
                result[key] = {"status": "PASS", "continuous_directory": str(full), "restart_directory": str(split),
                               "split_step": split_step, "completed_steps": 64, "field_max_abs_differences": diffs,
                               "checkpoint_bitwise_equal": fullpath.read_bytes() == splitpath.read_bytes(),
                               "history": "bitwise equal numerical values and absolute time columns",
                               "damaged_tail_byte_exact_backups": backups}
        label = "steady_sideheated"
        full = self.build("steady_restart_full", configure(self.source, label, steps=64))
        self.run(full)
        _, fullstate, _ = latest_checkpoint(full, 24)
        split = self.build("steady_restart_split35", configure(self.source, label, steps=35))
        self.run(split)
        for prefix in ("Nu_VolAvg_", "Re_VolAvg_"):
            for path in split.glob(prefix + "*.dat"):
                path.unlink()  # 仅删除本次隔离测试产生的最终诊断文件，模拟稳态中断。
        self.build("unused", configure(self.source, label, steps=64, restart=True), directory=split)
        self.run(split)
        _, splitstate, splitmeta = latest_checkpoint(split, 24)
        if int(splitmeta["itc_total"]) != 64:
            raise AssertionError("steady restart failed to reach step 64")
        for name in ("rho", "u", "v", "T"):
            np.testing.assert_array_equal(fullstate["fields"][name], splitstate["fields"][name])
        final_hist = histories(split)
        result["steady_without_final_NuRe_files"] = {
            "status": "PASS", "directory": str(split), "split_step": 35,
            "completed_steps": 64, "field_comparison": "bitwise identical to uninterrupted 64 steps",
            "regenerated_final_history_rows": int(final_hist["Nu_VolAvg_"].shape[0])}
        return result

    def postprocess(self):
        text = configure(self.source, "steady_sideheated", steps=64)
        ref = self.build("postprocess_default_scale", text)
        self.run(ref)
        changed = self.build("postprocess_freefall_scale", set_parameter(text, "velocityScaleCompare", "1.0d0/velocityUnit"))
        self.run(changed)
        _, a, _ = latest_checkpoint(ref, 24)
        _, b, _ = latest_checkpoint(changed, 24)
        for name in ("rho", "u", "v", "T"):
            np.testing.assert_array_equal(a["fields"][name], b["fields"][name])
        da, db = diagnostics(ref), diagnostics(changed)
        diffs = {}
        for name in ("Nu_global", "Nu_hot", "Nu_cold", "Nu_middle", "NuVolAvg"):
            if name not in da or name not in db:
                raise AssertionError(f"missing {name} for output-scale regression")
            diffs[name] = abs(da[name]-db[name])
            np.testing.assert_allclose(da[name], db[name], atol=1e-13, rtol=1e-13)
        np.testing.assert_array_equal(histories(ref)["Nu_VolAvg_"], histories(changed)["Nu_VolAvg_"])
        return {"Nu_invariant_under_velocity_output_scale": {
            "status": "PASS", "reference_directory": str(ref), "changed_directory": str(changed),
            "Nu_absolute_differences": diffs, "field_comparison": "bitwise identical",
            "velocity_output": {"default_u_mid_max": da["u_mid_max"], "freefall_u_mid_max": db["u_mid_max"]}}}

    def benchmark(self):
        # de Vahl Davis方腔值；测试前应与本次复现论文对应表核对。
        # 速度采用热扩散标度UL/kappa，不能与自由落体速度混比。
        references = {1000.0: {"Nu": 1.118, "u": 3.649, "v": 3.697},
                      10000.0: {"Nu": 2.243, "u": 16.178, "v": 19.617}}
        result = {}
        for rayleigh, reference in references.items():
            if rayleigh not in self.args.benchmark_rayleigh:
                continue
            for n in self.args.benchmark_grids:
                label = f"Ra{int(rayleigh)}_N{n}"
                text = configure(self.source, "steady_sideheated", n=n, rayleigh=rayleigh,
                                 steps=self.args.benchmark_steps, benchmark=True)
                directory = self.build("benchmark_" + label, text)
                self.run(directory)
                values = diagnostics(directory)
                if any(name not in values for name in ("Nu_hot", "Nu_cold", "u_mid_max", "v_mid_max")):
                    raise AssertionError("missing final benchmark diagnostics")
                convfiles = [p for p in directory.glob("convergence*.plt") if "all" not in p.name]
                converged = False
                last = None
                if convfiles:
                    rows = []
                    for line in convfiles[0].read_text(encoding="utf-8").splitlines():
                        row = line.split()
                        if len(row) == 3:
                            try:
                                rows.append([float(x.replace("D", "E")) for x in row])
                            except ValueError:
                                pass
                    if rows:
                        last = rows[-1]
                        converged = last[1] <= 1e-8 and last[2] <= 1e-8
                relative = {"Nu_hot": abs(values["Nu_hot"]/reference["Nu"]-1),
                            "u_mid_max": abs(values["u_mid_max"]/reference["u"]-1),
                            "v_mid_max": abs(values["v_mid_max"]/reference["v"]-1)}
                wall_imbalance = abs(values["Nu_hot"]-values["Nu_cold"])/max(abs(values["Nu_hot"]), 1e-30)
                _, final_state, _ = latest_checkpoint(directory, n)
                density = final_state["fields"]["rho"]
                state_text = "\n".join(path.read_text(encoding="utf-8", errors="replace")
                                       for path in directory.glob("SimulationSettings*.txt"))
                state_pattern = r"state_check itc/rho_mean/rho_min/rho_max/maxMa:\s*(\d+)" + (r"\s*(" + FLOAT_RE + ")")*4
                state_rows = [[float(x.replace("D", "E").replace("d", "e")) for x in row]
                              for row in re.findall(state_pattern, state_text)]
                state_ranges = {}
                if state_rows:
                    state_values = np.asarray(state_rows)
                    state_ranges = {"state_check_records": len(state_rows),
                                    "distinct_checked_steps": len(set(state_values[:, 0])),
                                    "rho_mean_min": float(state_values[:, 1].min()),
                                    "rho_mean_max": float(state_values[:, 1].max()),
                                    "rho_local_min_over_checks": float(state_values[:, 2].min()),
                                    "rho_local_max_over_checks": float(state_values[:, 3].max()),
                                    "max_local_Mach_over_checks": float(state_values[:, 4].max())}
                result[label] = {"status": "PASS" if converged and max(relative.values()) <= 0.05 and wall_imbalance < 0.01 else "FAIL",
                                 "directory": str(directory), "converged": converged, "last_convergence_row": last,
                                 "diagnostics": values, "reference": reference, "relative_errors": relative,
                                 "reference_source": BENCHMARK_REFERENCE,
                                 "Nu_global_relative_error": abs(values["Nu_global"]/reference["Nu"]-1),
                                 "Nu_middle_relative_error": abs(values["Nu_middle"]/reference["Nu"]-1),
                                 "hot_cold_relative_imbalance": wall_imbalance,
                                 "rho_final": {"mean": float(density.mean()), "min": float(density.min()), "max": float(density.max())},
                                 "state_check_ranges": state_ranges,
                                 "scope": "low-Ra benchmark agreement; does not establish high-Ra stability or grid independence"}
                print(f"{label}: {result[label]['status']}; converged={converged}; "
                      f"Nu_hot={values['Nu_hot']:.8f}; relative_error={relative['Nu_hot']:.4%}", flush=True)
        return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--compiler", default=shutil.which("gfortran") or "C:/msys64/ucrt64/bin/gfortran.exe")
    parser.add_argument("--stages", default="syntax,smoke,conduction,analytic,restart,postprocess")
    parser.add_argument("--work-root", type=Path, help="ASCII父目录；省略时使用系统临时目录")
    parser.add_argument("--timeout", type=int, default=240, help="每次编译/运行的超时秒数")
    parser.add_argument("--benchmark-grids", default="32,64")
    parser.add_argument("--benchmark-rayleigh", default="1000,10000", help="仅选择已有参考值的Ra，逗号分隔")
    parser.add_argument("--benchmark-steps", type=int, default=200000)
    args = parser.parse_args()
    args.benchmark_grids = [int(x) for x in args.benchmark_grids.split(",")]
    args.benchmark_rayleigh = [float(x) for x in args.benchmark_rayleigh.split(",")]
    if not args.benchmark_rayleigh or any(x not in (1000.0, 10000.0) for x in args.benchmark_rayleigh):
        parser.error("--benchmark-rayleigh must select 1000 and/or 10000")
    selected = [x.strip() for x in args.stages.split(",") if x.strip()]
    if not selected or any(x not in STAGES for x in selected):
        parser.error("--stages must use: " + ",".join(STAGES))
    if not args.source.is_file():
        parser.error(f"source not found: {args.source}")
    run = Validation(args)
    print("Build root:", run.root, flush=True)
    failed = False
    for stage in selected:
        print("Stage:", stage, flush=True)
        begin = time.monotonic()
        try:
            details = getattr(run, stage)()
            if stage == "benchmark":
                previous = run.report.get("tests", {}).get(stage, {}).get("cases", {})
                details = {**previous, **details}  # 追加单例时必须保留已有粗网格FAIL证据。
                for case in details.values():
                    case["reference_source"] = BENCHMARK_REFERENCE
            nested_failed = any(value.get("status") == "FAIL" for value in details.values())
            outcome = {"status": "FAIL" if nested_failed else "PASS", "cases": details}
            failed |= nested_failed
        except Exception as exc:
            outcome = {"status": "FAIL", "error": str(exc), "traceback": traceback.format_exc()}
            failed = True
        outcome["elapsed_seconds"] = round(time.monotonic() - begin, 3)
        run.report["tests"][stage] = outcome
        run.save()
        print(stage + ": " + outcome["status"], flush=True)
        if outcome["status"] == "FAIL":
            print(outcome.get("error", "one or more cases failed; see report"), flush=True)
    print("Report:", run.report_path, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
