"""Set concrete walltimes from the successful preflight; write serial PBS jobs."""
from pathlib import Path
import hashlib
import json
import math
import tarfile

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"运行脚本/run_cases/DUGKS_P100_node05_20261030"
cfg=json.loads((LOCAL/"campaign.json").read_text(encoding="utf-8"))
gate=json.loads((LOCAL/"preflight/assessment.json").read_text(encoding="utf-8"))
assert gate["passed"] and len(gate["cases"])==7
remote=cfg["remote_dir"]
timings={x["case"]:x for x in gate["cases"]}

template="""#!/usr/bin/env bash
#PBS -N {job_name}
#PBS -q batch
#PBS -l nodes=node05:ppn=1
#PBS -l walltime={walltime}
#PBS -o {directory}/results/pbs.stdout
#PBS -e {directory}/results/pbs.stderr
set -euo pipefail
CAMPAIGN={remote}
CASE_ROOT={directory}
RESULT="$CASE_ROOT/results"
BUILD="$CASE_ROOT/formal_build"
mkdir -p "$RESULT" "$BUILD"
cd "$RESULT"
test "$(hostname -s)" = node05
test ! -e run.started
date -Is > run.started
finish() {{ rc=$?; printf '%s\\n' "$rc" > "$RESULT/job.exit.status"; date -Is > "$RESULT/run.finished"; }}
trap finish EXIT
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3
export PYTHONHOME="$CAMPAIGN/runtime/python36"
export PATH="$PYTHONHOME/bin:$NVHPC_ROOT/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$PYTHONHOME/lib64:$NVHPC_ROOT/compilers/lib:$(printenv LD_LIBRARY_PATH || true)"
export OMP_NUM_THREADS=1
export ACC_DEVICE_TYPE=nvidia
export CUDA_VISIBLE_DEVICES="$(printenv CUDA_VISIBLE_DEVICES || printf 0)"
{{ date -Is; hostname; printf 'PBS_JOBID=%s\\nPBS_NODEFILE=%s\\n' "$PBS_JOBID" "$PBS_NODEFILE"; cat "$PBS_NODEFILE"; }} > job_identity.txt
nvidia-smi > gpu.start.txt 2>&1 || true
python3 "$CAMPAIGN/gpu_gate.py" > gpu_gate.formal.json
python3 -c 'import json,sys; a=json.load(open(sys.argv[1])); assert a["passed"]' "$RESULT/smoke_validation.json"
nvfortran --version > compiler_version.txt 2>&1
printf '%s\\n' '-cpp -O3 -acc -gpu=cc60 -Mextend -Minfo=accel' > build_flags.txt
sha256sum "$CASE_ROOT/source/solver.F90" > runtime_source.sha256
printf '%s  %s\\n' '{sha}' "$CASE_ROOT/source/solver.F90" | sha256sum -c -
cd "$BUILD"
set +e
nvfortran -cpp -O3 -acc -gpu=cc60 -Mextend -Minfo=accel "$CASE_ROOT/source/solver.F90" -o solver.exe > "$RESULT/compile.log" 2>&1
rc=$?
set -e
printf '%s\\n' "$rc" > "$RESULT/compile.status"
test "$rc" -eq 0
cd "$RESULT"
date -Is > solver.started
set +e
/usr/bin/time -v -o timing.txt "$BUILD/solver.exe" > solver.stdout 2> solver.stderr
rc=$?
set -e
printf '%s\\n' "$rc" > run.status
date -Is > solver.finished
nvidia-smi > gpu.end.txt 2>&1 || true
set +e
python3 "$CAMPAIGN/assess.py" --case "$CASE_ROOT" > assessment.json
assessment_rc=$?
set -e
printf '%s\\n' "$assessment_rc" > validation.status
test "$rc" -eq 0
exit "$assessment_rc"
"""
for case in cfg["cases"]:
    expected=timings[case["name"]]["estimated_1000_tff_hours"]
    # Long-run output/CPU diagnostics and shared-GPU contention need margin.
    hours=max(240,int(math.ceil((2*expected+24)/24))*24)
    case["walltime_hours"]=hours
    case["preflight_estimate_1000_tff_hours"]=expected
    case["walltime"]="%d:00:00" % hours
    text=template.format(job_name="DG_"+("S" if case["case_type"]=="side" else "R")+str(int(math.log10(case["Ra"]))),
        directory=case["remote_dir"],remote=remote,walltime=case["walltime"],sha=case["source_sha256"])
    directory=LOCAL/case["name"]
    (directory/"run.pbs").write_text(text,encoding="utf-8",newline="\n")
    (directory/"case.json").write_text(json.dumps(case,indent=2)+"\n",encoding="utf-8",newline="\n")
cfg["preflight_passed"]=True
cfg["formal_walltime_policy"]="max(240h, round_up_24h(2*warm_start_estimate+24h)); not a guaranteed completion prediction"
(LOCAL/"campaign.json").write_text(json.dumps(cfg,indent=2,ensure_ascii=False)+"\n",encoding="utf-8",newline="\n")

submit='''"""Submit once, checking current PBS and both GPUs before each qsub."""
from pathlib import Path
import datetime,json,re,subprocess
root=Path(__file__).resolve().parent
cfg=json.loads((root/'campaign.json').read_text())
assert cfg['preflight_passed'] and cfg['node']=='node05'
marker=root/'submission.started'
with marker.open('x') as file:
    file.write(datetime.datetime.now().isoformat()+'\\n')
submitted=[]
previous=None
for case in cfg['cases']:
    work=root/case['name']
    commands=[['date','-Is'],['qstat','-a'],['pbsnodes','node05','node07']]
    for node in ['node05','node07']:
        commands.append(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10',node,
            'date -Is; nvidia-smi; nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv'])
    observations=[]
    for command in commands:
        p=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        observations.append('COMMAND: '+' '.join(command)+'\\n'+p.stdout)
        if p.returncode:
            raise RuntimeError('Node pre-submit check failed: '+repr(command))
    (work/'pre_submit_node_check.txt').write_text('\\n'.join(observations))
    (work/'node_selection_reason.txt').write_text('User explicitly selected node05 and authorized sharing; GPU free-memory gate is checked at startup. This case starts only after the preceding independent DUGKS job ends. No other task is modified.\\n')
    subprocess.check_call(['bash','-n',str(work/'run.pbs')])
    args=['qsub']
    if previous:
        args+=['-W','depend=afterany:'+previous]
    args.append(str(work/'run.pbs'))
    output=subprocess.check_output(args,universal_newlines=True).strip()
    if not re.fullmatch(r'\\d+\\.master',output):
        raise RuntimeError('Unexpected qsub response: '+output)
    item=dict(case=case['name'],job_id=output,dependency=previous,node='node05',
        walltime=case['walltime'],directory=str(work),source_sha256=case['source_sha256'])
    submitted.append(item)
    with (root/'submitted.tsv').open('a') as file:
        file.write('\\t'.join([case['name'],output,previous or '-',case['walltime'],str(work)])+'\\n')
    (root/'submission.json').write_text(json.dumps(submitted,indent=2)+'\\n')
    detail=subprocess.check_output(['qstat','-f',output],universal_newlines=True)
    (work/'qstat_after_submit.txt').write_text(detail)
    assert 'Resource_List.nodes = node05:ppn=1' in detail
    print(json.dumps(item),flush=True)
    previous=output
(root/'submission.finished').write_text(datetime.datetime.now().isoformat()+'\\n')
'''
(LOCAL/"submit_chain.py").write_text(submit,encoding="utf-8",newline="\n")
archive=ROOT/"tmp/dugks_p100_20261030/formal_scripts.tar.gz"
with tarfile.open(archive,"w:gz") as tar:
    for name in ["campaign.json","submit_chain.py","assess.py","gpu_gate.py"]:
        tar.add(LOCAL/name,arcname=name)
    for case in cfg["cases"]:
        for name in ["run.pbs","case.json"]:
            tar.add(LOCAL/case["name"]/name,arcname=case["name"]+"/"+name)
print(json.dumps([dict(case=x["name"],estimated_hours=x["preflight_estimate_1000_tff_hours"],walltime=x["walltime"])
                  for x in cfg["cases"]],indent=2))
