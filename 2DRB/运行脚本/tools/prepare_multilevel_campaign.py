"""Freeze the seven user-requested original-Luo, three-level P100 cases."""
from pathlib import Path
import hashlib
import json
import re
import tarfile

ROOT = Path(__file__).resolve().parents[2]
NAME = 'Multilevel_original_L3_20261027_v3'
OUT = ROOT / '运行脚本' / 'run_cases' / NAME
REMOTE = '/data2/XLLi/Multiblock/' + NAME
SOURCE = ROOT / '多块网格' / '2DRBOpenaccMultilevel.F90'
POST = ROOT / '多块网格' / '2DRBOpenaccMultilevelPostprocess.F90'
raw = SOURCE.read_bytes()
post_raw = POST.read_bytes()
sha = lambda data: hashlib.sha256(data).hexdigest()
text = raw.decode('utf-8-sig').replace('\r\n', '\n')
post_text = post_raw.decode('utf-8-sig').replace('\r\n', '\n')


def macro(s, name, enabled):
    s, count = re.subn(r'^!?#define ' + name + r'[ \t]*$',
                       ('#define ' if enabled else '!#define ') + name, s, flags=re.M)
    assert count == 1, (name, count)
    return s


def parameter(s, name, value):
    s, count = re.subn(r'(\b' + name + r'\s*=\s*)[\d.]+(?:[dDeE][+-]?\d+)?',
                       lambda m: m[1] + str(value), s, count=1)
    assert count == 1, name
    return s


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data.encode('utf-8') if isinstance(data, str) else data)


def create_case(name, s, post, meta, diagnostic=False):
    folder = OUT / name
    folder.mkdir(parents=True, exist_ok=False)
    for child in ['source', 'build', 'results']:
        (folder / child).mkdir()
    data = s.encode('utf-8')
    write(folder / 'source/solver.F90', data)
    write(folder / 'source/local_master.F90', raw)
    write(folder / 'source/postprocess.F90', post)
    meta.update(case=name, remote=REMOTE + '/' + name,
                local_master_sha256=sha(raw), case_source_sha256=sha(data),
                local_postprocess_sha256=sha(post_raw),
                postprocess_case_sha256=sha(post.encode('utf-8')),
                numerical_kernels_changed=False, loadInitField=0,
                nGridLevels=3, adjacent_ratio=2, node='node05',
                start='from_zero', diagnostic=diagnostic)
    write(folder / 'manifest.json', json.dumps(meta, indent=2, ensure_ascii=False) + '\n')
    run = r'''#!/bin/bash
set -euo pipefail
CASE='__CASE__'
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
cd "$CASE/results"
test ! -e run.started
date -Is > run.started
trap 'rc=$?; echo "$rc" > "$CASE/results/wrapper.status"; date -Is > "$CASE/results/wrapper.finished"' EXIT
test "$(hostname -s)" = node05
printf 'job=%s\nhost=%s\ncuda_visible_devices=%s\n' "${PBS_JOBID:-unknown}" "$(hostname)" "$CUDA_VISIBLE_DEVICES" > job_identity.txt
sha256sum "$CASE/source/solver.F90" > runtime_source.sha256
test "$(cut -d ' ' -f1 runtime_source.sha256)" = '__HASH__'
nvfortran --version > compiler_version.txt 2>&1
printf '%s\n' '-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend' > build_flags.txt
cd "$CASE/build"
set +e
nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend "$CASE/source/solver.F90" -o solver.exe > "$CASE/results/compile.log" 2>&1
rc=$?
set -e
echo "$rc" > "$CASE/results/compile.status"
test "$rc" -eq 0 || exit "$rc"
cd "$CASE/results"
nvidia-smi > nvidia_smi_before.log 2>&1 || true
date -Is > solver.started
set +e
"$CASE/build/solver.exe" > solver.stdout 2>&1
rc=$?
set -e
echo "$rc" > run.status
date -Is > run.finished
test "$rc" -eq 0 || { echo runtime_failed > validation.status; exit "$rc"; }
if grep -Eqi '^[[:space:]]*(Error:|ERROR STOP|FATAL)' solver.stdout; then
    echo fatal_error > validation.status; exit 87
fi
if grep -Eqi '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)' solver.stdout ./*.dat; then
    echo nonfinite > validation.status; exit 86
fi
grep -q 'Completed total fine steps:' solver.stdout
test -s FinalNuRe_2DOpenaccMultiblock.dat
echo finite_run_completed > validation.status
nvidia-smi > nvidia_smi_after.log 2>&1 || true
__POSTPROCESS__
'''.replace('__CASE__', meta['remote']).replace('__HASH__', sha(data))
    if diagnostic:
        run = run.replace('__POSTPROCESS__', '')
    else:
        run = run.replace('__POSTPROCESS__', r'''cd "$CASE/build"
set +e
nvfortran -cpp -O2 -Mextend "$CASE/source/postprocess.F90" -o postprocess.exe > "$CASE/results/postprocess_compile.log" 2>&1
rc=$?
set -e
echo "$rc" > "$CASE/results/postprocess_compile.status"
test "$rc" -eq 0 || exit "$rc"
cd "$CASE/results"
set +e
"$CASE/build/postprocess.exe" > postprocess.stdout 2>&1
rc=$?
set -e
echo "$rc" > postprocess.status
exit "$rc"''')
    write(folder / 'run.sh', run)
    return meta


assert not OUT.exists(), 'Do not overwrite a prepared or submitted campaign'
OUT.mkdir(parents=True)
cases = []
specs = [('steady', 6, 1024), ('steady', 7, 1024), ('steady', 8, 1024),
         ('unsteady', 6, 1024), ('unsteady', 7, 1024),
         ('unsteady', 8, 4096), ('unsteady', 9, 4096)]
for mode, exponent, n in specs:
    steady = mode == 'steady'
    name = f'{mode}_Ra1e{exponent}_N{n}_L3'
    s = text
    for key, enabled in [('steadyFlow', steady), ('unsteadyFlow', not steady),
                         ('SideHeatedCell', steady), ('RayleighBenardCell', not steady),
                         ('HorizontalWallsConstT', not steady), ('VerticalWallsAdiabatic', not steady),
                         ('HorizontalWallsAdiabatic', steady), ('VerticalWallsConstT', steady)]:
        s = macro(s, key, enabled)
    changes = dict(nx=n, ny=n, nGridLevels=3, Rayleigh=f'1.0d{exponent}',
                   Prandtl='0.71d0' if steady else '0.7d0', Mach='0.1d0',
                   finestLayerCellsLeftRight=n // 16, finestLayerCellsBottomTop=n // 16,
                   innerLayerCellsLeftRight=n // 32, innerLayerCellsBottomTop=n // 32)
    for key, value in changes.items():
        s = parameter(s, key, value)
    post = macro(macro(post_text, 'steadyFlow', steady), 'unsteadyFlow', not steady)
    cases.append(create_case(name, s, post, dict(mode=mode,
        geometry='SideHeatedCell' if steady else 'RayleighBenardCell',
        ra=10**exponent, pr=0.71 if steady else 0.7, mach=0.1,
        finest_n=n, equivalent_grids=[n, n//2, n//4],
        fine_band_cells=n//16, middle_band_local_cells=n//32,
        source_parameter_changes=changes,
        epsU=1e-7, epsT=1e-7, steady_max_steps=20000000 if steady else None,
        steady_check_steps=2000, target_tff=None if steady else 1000,
        sample_and_snapshot_tff=10 if steady else 0.5,
        checkpoint_tff=100, tecplot_tff=100,
        statistics_tff=None if steady else [500, 750, 1000],
        walltime='120:00:00', relaxation='original raw scaling from finest Luo',
        regularization=False)))

# A separate short, synchronized GPU check precedes production submissions.
diagnostics = []
for base, n, steps in [(cases[0], 128, 40), (cases[3], 128, 40),
                        (cases[3], 1024, 8), (cases[5], 4096, 8)]:
    s = (OUT / base['case'] / 'source/solver.F90').read_text(encoding='utf-8')
    for key, value in dict(nx=n, ny=n, finestLayerCellsLeftRight=max(8,n//16),
                           finestLayerCellsBottomTop=max(8,n//16),
                           innerLayerCellsLeftRight=max(8,n//32),
                           innerLayerCellsBottomTop=max(8,n//32),
                           outputSnapshotFile=0, outputPltFile=0, outputReloadFile=0).items():
        s = parameter(s, key, value)
    old = 'do while(reloadItcOffset+itc<itc_max)'
    assert s.count(old) == 1
    s = s.replace(old, f'do while(reloadItcOffset+itc<min(itc_max,{steps}))')
    name = f'preflight/{base["mode"]}_N{n}_{steps}steps'
    diagnostics.append(create_case(name, s, post_text,
        dict(mode=base['mode'], finest_n=n, expected_fine_steps=steps,
             source_case=base['case'], synchronized=True), diagnostic=True))


def pbs(name, jobname, command, walltime):
    folder = OUT / name
    (folder / 'results').mkdir(parents=True, exist_ok=True)
    remote = REMOTE + '/' + name
    data = f'''#!/bin/bash
#PBS -N {jobname}
#PBS -q batch
#PBS -l nodes=node05:ppn=1
#PBS -l walltime={walltime}
#PBS -o {remote}/results/pbs.stdout
#PBS -e {remote}/results/pbs.stderr
set -euo pipefail
test "$(hostname -s)" = node05
{command}
'''
    write(folder / 'run.pbs', data)


pbs('preflight', 'ML3_GPUcheck', 'export NVCOMPILER_ACC_SYNCHRONOUS=1\n' + '\n'.join(
    f'bash {d["remote"]}/run.sh' for d in diagnostics), '00:30:00')
for case in cases:
    jobname = f'ML3_{"S" if case["mode"] == "steady" else "U"}_R{len(str(case["ra"]))-1}'
    pbs(case['case'], jobname, 'bash ' + case['remote'] + '/run.sh', case['walltime'])

manifest = dict(name=NAME, remote=REMOTE, master_source=str(SOURCE),
    master_sha256=sha(raw), postprocess_sha256=sha(post_raw),
    cases=cases, diagnostics=diagnostics,
    order='three side-heated steady, then four RB unsteady; independent afterany chain',
    chosen_node='node05', selection_initial_observation={
        'time_server':'2026-10-27T22:48:27+08:00',
        'node05':{'pbs_job':'7249.master', 'gpu_util_percent':0, 'gpu_used_MiB':448,
                  'gpu_free_MiB':15821, 'other_gpu_process_MiB':446},
        'node07':{'pbs_job':'7400.master', 'gpu_util_percent':100, 'gpu_used_MiB':968,
                  'gpu_free_MiB':15301},
        'reason':'node05 lower GPU occupancy; preserve existing processes; serialize campaign'},
    submission_requires_fresh_probe=True)
write(OUT / 'manifest.json', json.dumps(manifest, indent=2, ensure_ascii=False) + '\n')
write(OUT / 'README.md', '''# 原始 Luo 三级网格 P100 算例

本目录是本地母版的冻结快照。最细级使用原始 Luo 参数，其他级别采用原始粗细网格缩放，无 RR/HRR。

稳态：侧壁加热，Pr=0.71，Ra=1e6、1e7、1e8，1024/512/256。
非稳态：RB，Pr=0.7，Ra=1e6、1e7 用 1024/512/256，Ra=1e8、1e9 用 4096/2048/1024；目标 1000 TFF。
非稳态后处理窗口 500–1000 TFF，同时比较 500–750 与 750–1000。

输出频率、文件名、列格式保持母版设置。完整快照预计产生约 0.8 TB 数据；提交时另查可用空间。
GPU 短测通过后才提交正式任务；正式任务以 afterany 串行依赖，独立算例失败不阻止后续算例。
每个任务申请 120 小时。达到时间限制不等于完成物理时长或收敛，需依现有有限检查点另行续算。
manifest.json 保存参数及三层哈希的前两层；results/runtime_source.sha256 保存实际编译源哈希。
submission.json 和 node_checks/ 保存提交编号、依赖及每次提交前的实时节点检查。
''')
with tarfile.open(OUT / 'cases.tar', 'w') as archive:
    for p in sorted(OUT.iterdir()):
        if p.name != 'cases.tar':
            archive.add(p, arcname=NAME + '/' + p.name)
print(json.dumps(dict(local=str(OUT), remote=REMOTE,
                      cases=len(cases), diagnostics=len(diagnostics), master_sha256=sha(raw))))
