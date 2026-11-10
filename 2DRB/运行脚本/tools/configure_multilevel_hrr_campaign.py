"""仅为本次全新算例生成冻结的P100脚本与GPU短测输入。"""
from pathlib import Path
import hashlib
import json
import re
import shutil

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'运行脚本/run_cases/Multilevel_HRR0995_Ra9_Ra10_N2048_20261103'
WORK=Path('C:/Users/User/.codex/visualizations/2026/10/26/01a176f3-942f-7de2-84f7-0618b0d3d607/multilevel_hrr_20261103')


def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text,encoding='utf-8',newline='\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assign(text,name,value):
    text,n=re.subn(r'(\b'+name+r'=)[\d.]+(?:d[\d+-]+)?',lambda m:m[1]+str(value),text,count=1)
    assert n==1,name
    return text


def main():
    m=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
    proof=json.loads((OUT/'local_validation.json').read_text(encoding='utf-8'))
    assert proof['status']=='passed'
    assert not (OUT/'submission.json').exists(), 'Do not alter a submitted campaign'
    for c in m['cases']: assert not (OUT/c['case']/'job.id').exists()
    small=OUT/'verification/small/solver.F90'
    small.parent.mkdir(exist_ok=True)
    shutil.copyfile(WORK/'levels3/solver.F90',small)
    shutil.copyfile(WORK/'levels3/raw-state.bin',OUT/'verification/small/cpu-state.bin')
    count=(OUT/'verification/small/cpu-state.bin').stat().st_size//8
    write(OUT/'verification/validate_smoke.F90',f'''program validate_smoke
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
integer,parameter :: count={count}
real(8),allocatable :: actual(:),reference(:)
real(8) :: error
allocate(actual(count),reference(count))
open(10,file='raw-state.bin',access='stream',form='unformatted',status='old')
read(10)actual
close(10)
open(10,file='cpu-state.bin',access='stream',form='unformatted',status='old')
read(10)reference
close(10)
if(.not.all(ieee_is_finite(actual))) error stop 'Nonfinite GPU smoke field/distribution'
error=maxval(abs(actual-reference)/max(1.0d0,abs(reference)))
if(error>1d-10)error stop 'CPU/GPU full-state comparison failed'
write(*,*) 'PASS CPU/GPU f,g,rho,u,v,T,Fx,Fy comparison; max scaled difference:',error
end program
''')
    diagnostics=[dict(case='verification/small',remote=m['remote']+'/verification/small',
        expected_fine_steps=32,source_sha256=sha(small),full_state_values=count)]
    for c in m['cases']:
        c.update(node='node05',walltime='120:00:00',permitted_existing_gpu_pids=[4112],
            node_reason='Live 2026-11-03 check: node07 GPU100% with FPC job7524; node05 GPU0% in two samples, existing PID4112 uses446MiB, 15821MiB free. Choose lower GPU occupancy and serialize the two new cases; preserve all existing jobs.')
        if c['Ra']==10**9:
            c['predecessor_case']='Ra1e10_N2048_L3_HRR0995'
        folder=OUT/c['case']; text=(folder/'source/solver.F90').read_text(encoding='utf-8')
        assert sha(folder/'source/solver.F90')==c['source_sha256']
        library,n=re.subn(r'^    program main\b.*?^    end program main\b','',text,flags=re.M|re.S)
        assert n==1
        write(folder/'source/solver_library.F90',library)
        actual=text.replace('do while(reloadItcOffset+itc<itc_max)',
                            'do while(reloadItcOffset+itc<min(itc_max,8))')
        for key in ['outputSnapshotFile','outputPltFile','outputReloadFile']:
            actual=assign(actual,key,0)
        write(folder/'source/actual_grid_smoke.F90',actual)
        c['library_sha256']=sha(folder/'source/solver_library.F90')
        c['actual_smoke_sha256']=sha(folder/'source/actual_grid_smoke.F90')
        diagnostics.append(dict(case=c['case']+'/source/actual_grid_smoke',remote=c['remote']+'/source',
            expected_fine_steps=8,source_filename='actual_grid_smoke.F90',source_sha256=c['actual_smoke_sha256']))
        script=r'''#!/bin/bash
set -euo pipefail
CASE='__CASE__'
CAMPAIGN='__CAMPAIGN__'
EXPECTED_NODE=node05
NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 ACC_DEVICE_TYPE=nvidia
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
mkdir -p "$CASE/build/main" "$CASE/build/kernel" "$CASE/build/small" "$CASE/build/actual"
cd "$CASE/results"
test ! -e run.started
date -Is > run.started
phase=preflight
finish() {
    rc=$?
    echo "$rc" > "$CASE/results/wrapper.status"
    if test "$phase" = preflight && test "$rc" -ne 0; then echo failed > "$CASE/results/preflight.status"; fi
    date -Is > "$CASE/results/wrapper.finished"
}
trap finish EXIT
test "$(hostname -s)" = "$EXPECTED_NODE"
printf 'job=%s\nhost=%s\ncuda_visible_devices=%s\n' "${PBS_JOBID:-unknown}" "$(hostname -s)" "$CUDA_VISIBLE_DEVICES" > job_identity.txt
sha256sum "$CASE/source/solver.F90" > runtime_source.sha256
test "$(cut -d ' ' -f1 runtime_source.sha256)" = '__SOURCE_HASH__'
echo '__LIBRARY_HASH__  __CASE__/source/solver_library.F90' | sha256sum -c - > preflight_hashes.log
echo '__ACTUAL_HASH__  __CASE__/source/actual_grid_smoke.F90' | sha256sum -c - >> preflight_hashes.log
echo '__SMALL_HASH__  __CAMPAIGN__/verification/small/solver.F90' | sha256sum -c - >> preflight_hashes.log
nvfortran --version > compiler_version.txt 2>&1
nvidia-smi > nvidia_smi_before.log 2>&1
printf '%s\n' '-cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend' > build_flags.txt
cd "$CASE/build/main"
set +e
nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend "$CASE/source/solver.F90" -o solver.exe > "$CASE/results/compile.log" 2>&1
rc=$?
set -e
echo "$rc" > "$CASE/results/compile.status"
test "$rc" -eq 0 || exit "$rc"
# 非零力、速度梯度及接口回环；对照已有独立Hermite参考值。
cd "$CASE/build/kernel"
nvfortran -cpp -DAUDIT_GPU -O3 -acc -gpu=cc60 -Mextend "$CASE/source/solver_library.F90" "$CAMPAIGN/verification/verify.F90" -o verify.exe > "$CASE/results/kernel_compile.log" 2>&1
mkdir "$CASE/results/kernel"
cp "$CAMPAIGN/verification/input.dat" "$CAMPAIGN/verification/expected.dat" "$CASE/results/kernel/"
cd "$CASE/results/kernel"
NVCOMPILER_ACC_SYNCHRONOUS=1 "$CASE/build/kernel/verify.exe" > "$CASE/results/kernel_verify.stdout" 2>&1
cp verification_error.dat "$CASE/results/verification_error.dat"
# 小网格完整三级推进，全部f/g和宏观数组与GCC主机参考逐值比较。
cd "$CASE/build/small"
nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend "$CAMPAIGN/verification/small/solver.F90" -o small.exe > "$CASE/results/small_compile.log" 2>&1
nvfortran -O2 -Mextend "$CAMPAIGN/verification/validate_smoke.F90" -o compare.exe > "$CASE/results/small_compare_compile.log" 2>&1
mkdir "$CASE/results/small"
cp "$CAMPAIGN/verification/small/cpu-state.bin" "$CASE/results/small/"
cd "$CASE/results/small"
NVCOMPILER_ACC_SYNCHRONOUS=1 "$CASE/build/small/small.exe" > solver.stdout 2>&1
grep -q 'Completed total fine steps: 32; levels: 3' solver.stdout
"$CASE/build/small/compare.exe" > "$CASE/results/small_comparison.stdout" 2>&1
# 正式Ra及2048/1024/512实际几何短测：包含所有接口和真实显存分配。
cd "$CASE/build/actual"
nvfortran -cpp -O3 -acc -gpu=cc60 -Minfo=accel -Mextend "$CASE/source/actual_grid_smoke.F90" -o actual.exe > "$CASE/results/actual_compile.log" 2>&1
mkdir "$CASE/results/actual"
cd "$CASE/results/actual"
NVCOMPILER_ACC_SYNCHRONOUS=1 "$CASE/build/actual/actual.exe" > solver.stdout 2>&1
grep -q 'Completed total fine steps: 8; levels: 3' solver.stdout
test -s FinalNuRe_2DOpenaccMultiblock.dat
echo passed > "$CASE/results/preflight.status"
phase=production
cd "$CASE/results"
date -Is > solver.started
set +e
"$CASE/build/main/solver.exe" > solver.stdout 2>&1
rc=$?
set -e
echo "$rc" > run.status
date -Is > run.finished
test "$rc" -eq 0 || exit "$rc"
grep -q 'Completed total fine steps:' solver.stdout
test -s FinalNuRe_2DOpenaccMultiblock.dat
if grep -Eqi '(^|[^[:alpha:]])(nan|infinity|inf)([^[:alpha:]]|$)' solver.stdout NuRe_2DOpenaccMultiblock.dat FinalNuRe_2DOpenaccMultiblock.dat; then
    echo nonfinite > validation.status
    exit 86
fi
echo finite_run_completed > validation.status
nvidia-smi > nvidia_smi_after.log 2>&1
cd "$CASE/build/main"
nvfortran -cpp -O2 -Mextend "$CASE/source/postprocess.F90" -o postprocess.exe > "$CASE/results/postprocess_compile.log" 2>&1
cd "$CASE/results"
set +e
"$CASE/build/main/postprocess.exe" > postprocess.stdout 2>&1
rc=$?
set -e
echo "$rc" > postprocess.status
exit "$rc"
'''
        for token,value in dict(CASE=c['remote'],CAMPAIGN=m['remote'],SOURCE_HASH=c['source_sha256'],
                LIBRARY_HASH=c['library_sha256'],ACTUAL_HASH=c['actual_smoke_sha256'],SMALL_HASH=sha(small)).items():
            script=script.replace('__'+token+'__',value)
        assert '__' not in script
        write(folder/'run.sh',script)
        write(folder/'run.pbs',f'''#!/bin/bash
#PBS -N MLH_R{9 if c['Ra']==10**9 else 10}_2048
#PBS -q batch
#PBS -l nodes=node05:ppn=1
#PBS -l walltime=120:00:00
#PBS -o {c['remote']}/results/pbs.stdout
#PBS -e {c['remote']}/results/pbs.stderr
set -euo pipefail
bash {c['remote']}/run.sh
''')
        (folder/'results').mkdir(exist_ok=True)
    m['diagnostics']=diagnostics
    m['submission_order']=['Ra1e10_N2048_L3_HRR0995','Ra1e9_N2048_L3_HRR0995']
    m['submission_order_reason']='Check the more demanding Ra10 first; Ra9 follows afterany in its own from-zero directory.'
    write(OUT/'manifest.json',json.dumps(m,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status='configured',node='node05',order=m['submission_order'],
        snapshot_GiB_per_case=m['cases'][0]['estimated_snapshot_payload_bytes']/2**30),ensure_ascii=False))


if __name__=='__main__': main()
