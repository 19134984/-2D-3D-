"""P100 campaign transport; credentials are read only from temporary SSHPASS."""
from pathlib import Path
import datetime
import io
import json
import os
import re
import shlex
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '运行脚本/run_cases/Multilevel_original_L3_20261027_v3'
MANIFEST = json.loads((OUT / 'manifest.json').read_text(encoding='utf-8'))
REMOTE = MANIFEST['remote']
SSH = ['C:/msys64/usr/bin/sshpass.exe', '-e', 'C:/msys64/usr/bin/ssh.exe',
       '-o', 'ConnectTimeout=15', 'xlli@10.169.174.118']


def remote(command, data=None, check=True):
    assert os.environ.get('SSHPASS'), 'Provide a temporary session credential'
    result = subprocess.run(SSH + [command], input=data, capture_output=True)
    if check and result.returncode:
        raise RuntimeError(result.stdout.decode(errors='replace') + result.stderr.decode(errors='replace'))
    return result


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def collect_results(case, names):
    folder = REMOTE + '/' + case + '/results'
    quoted = ' '.join(shlex.quote(n) for n in names)
    command = ('cd ' + shlex.quote(folder) +
               '; for f in ' + quoted + '; do test ! -f "$f" || printf "%s\\0" "$f"; done | tar --null -T - -cf -')
    result = remote(command)
    target = OUT / case / 'results'
    target.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(result.stdout), mode='r:') as archive:
        for member in archive:
            assert member.isfile() and member.name in names, member.name
            (target / member.name).write_bytes(archive.extractfile(member).read())


def probe(label):
    # Query all users' jobs; a PBS "free" state alone does not mean an idle GPU.
    command = '''date -Is
hostname
qstat -n
pbsnodes node05 node07
df -h /data2/XLLi
for n in node05 node07; do
    echo GPU_NODE=$n
    ssh -o BatchMode=yes -o ConnectTimeout=10 "$n" 'nvidia-smi --query-gpu=index,name,uuid,utilization.gpu,memory.used,memory.free,memory.total --format=csv; nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv; sleep 1; nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.free --format=csv'
done
'''
    result = remote(command)
    text = result.stdout.decode(errors='replace')
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    record = dict(timestamp_utc=stamp, label=label, output=text)
    save_json(OUT / 'node_checks' / (label + '.json'), record)
    # Fail if either node's complete GPU probe failed or available memory is low.
    for node in ['node05', 'node07']:
        section = text.split('GPU_NODE=' + node + '\n', 1)[1].split('GPU_NODE=', 1)[0]
        match = re.search(r'^0, Tesla P100.*?, (\d+) %, (\d+) MiB, (\d+) MiB,', section, re.M)
        assert match, (node, section)
        record[node] = dict(gpu_percent=int(match[1]), used_MiB=int(match[2]), free_MiB=int(match[3]))
    assert record['node05']['free_MiB'] >= 8192, 'node05 no longer has sufficient free GPU memory'
    # Allow our running predecessor on node05, but require a low-occupancy initial choice.
    if label == 'preflight':
        assert record['node05']['gpu_percent'] < 20, 'Initial node05 GPU is now busy'
    record['chosen_node'] = 'node05'
    record['reason'] = ('node05 had the lower external GPU occupancy and sufficient memory; '
                        'serialize all campaign jobs on this node, leaving existing jobs intact')
    save_json(OUT / 'node_checks' / (label + '.json'), record)
    remote('cat > ' + shlex.quote(REMOTE + '/node_check_' + label + '.json'),
           json.dumps(record, indent=2).encode('utf-8'))
    print(json.dumps({k:v for k,v in record.items() if k != 'output'}, ensure_ascii=False), flush=True)
    return record


def submit(name, label, dependency=None):
    jobfile = OUT / name / 'job.id'
    dep = (' -W depend=afterany:' + dependency) if dependency else ''
    case = REMOTE + '/' + name
    # Write the identifier remotely in the same connection as qsub. An unknown
    # connection outcome is recovered by reading job.id, never by blind resubmit.
    command = ('cd ' + shlex.quote(case) + '; test ! -s job.id || exit 91; '
               'job=$(qsub' + dep + ' run.pbs) || exit $?; '
               'printf "%s\\n" "$job" > job.id; printf "%s\\n" "$job"')
    if jobfile.exists():
        job = jobfile.read_text().strip()
        recovered = remote('cat ' + shlex.quote(case + '/job.id'), check=False)
        assert job == recovered.stdout.decode().strip(), 'Local and remote job IDs differ'
    else:
        assert remote('test ! -s ' + shlex.quote(case + '/job.id'), check=False).returncode == 0
        probe(label)
        result = remote(command, check=False)
        job = result.stdout.decode().strip()
        if not re.fullmatch(r'\d+\.master', job):
            recovered = remote('cat ' + shlex.quote(case + '/job.id'), check=False)
            job = recovered.stdout.decode().strip()
            assert re.fullmatch(r'\d+\.master', job), result.stderr.decode(errors='replace')
    jobfile.write_text(job + '\n', encoding='ascii')
    status = remote('qstat -f ' + job).stdout.decode()
    assert re.search(r'Resource_List\.nodes\s*=\s*node05:ppn=1', status), status
    if dependency:
        assert 'afterany:' + dependency in status, status
    assigned = re.search(r'exec_host\s*=\s*(\S+)', status)
    if assigned:
        assert assigned[1].split('/')[0] == 'node05', status
    (OUT / name / 'qstat_submission.txt').write_text(status, encoding='utf-8')
    item = dict(case=name, job_id=job, dependency_afterany=dependency,
                requested_node='node05', actual_node=assigned[1] if assigned else None,
                remote=case, state=re.search(r'job_state\s*=\s*(\S+)', status)[1])
    existing = OUT / 'submission.json'
    submitted = json.loads(existing.read_text(encoding='utf-8')) if existing.exists() else []
    submitted = [s for s in submitted if s['job_id'] != job]
    submitted.append(item)
    save_json(existing, submitted)
    remote('cat > ' + shlex.quote(REMOTE + '/submission.json'), existing.read_bytes())
    print(json.dumps(item), flush=True)
    return job


action = sys.argv[1]
if action == 'upload':
    remote('mkdir -p /data2/XLLi/Multiblock; test ! -e ' + shlex.quote(REMOTE) +
           ' && tar -xf - -C /data2/XLLi/Multiblock', (OUT / 'cases.tar').read_bytes())
    paths = [d['remote'] for d in MANIFEST['cases'] + MANIFEST['diagnostics']]
    script = 'set -eu\n'
    for folder in paths:
        script += f'dos2unix {folder}/run.sh >/dev/null 2>&1\nbash -n {folder}/run.sh\n'
    for folder in [REMOTE + '/preflight'] + [d['remote'] for d in MANIFEST['cases']]:
        script += f'dos2unix {folder}/run.pbs >/dev/null 2>&1\nbash -n {folder}/run.pbs\n'
    for case in MANIFEST['cases'] + MANIFEST['diagnostics']:
        script += f'echo "{case["case_source_sha256"]}  {case["remote"]}/source/solver.F90" | sha256sum -c -\n'
    result = remote('bash -s', script.encode())
    (OUT / 'upload_verification.txt').write_bytes(result.stdout)
    print(result.stdout.decode())
elif action == 'preflight':
    submit('preflight', 'preflight')
elif action == 'preflight-status':
    job = (OUT / 'preflight/job.id').read_text().strip()
    result = remote('qstat -f ' + job, check=False)
    print(result.stdout.decode())
    for diagnostic in MANIFEST['diagnostics']:
        folder = diagnostic['remote'] + '/results'
        command = (f'echo CASE={diagnostic["case"]}; ls {folder}; '
                   f'cat {folder}/compile.status {folder}/run.status {folder}/wrapper.status 2>/dev/null; '
                   f'tail -8 {folder}/compile.log 2>/dev/null; tail -5 {folder}/solver.stdout 2>/dev/null')
        print(remote(command, check=False).stdout.decode())
elif action == 'submit-cases':
    for diagnostic in MANIFEST['diagnostics']:
        folder = diagnostic['remote'] + '/results'
        command = (f'test "$(cat {folder}/wrapper.status)" = 0 && '
                   f'test "$(cat {folder}/run.status)" = 0 && '
                   f'grep -q "Completed total fine steps: {diagnostic["expected_fine_steps"]};" {folder}/solver.stdout')
        remote(command)
    previous = None
    for index, case in enumerate(MANIFEST['cases'], 1):
        previous = submit(case['case'], f'{index:02d}_{case["case"]}', previous)
elif action == 'status':
    submissions = json.loads((OUT / 'submission.json').read_text(encoding='utf-8'))
    for item in submissions:
        job = item['job_id']
        result = remote('qstat -f ' + job, check=False)
        write_to = OUT / item['case'] / 'qstat_latest.txt'
        write_to.write_bytes(result.stdout + result.stderr)
        status = result.stdout.decode()
        assigned = re.search(r'exec_host\s*=\s*(\S+)', status)
        state = re.search(r'job_state\s*=\s*(\S+)', status)
        if assigned:
            assert assigned[1].split('/')[0] == 'node05', status
            item['actual_node'] = assigned[1]
        item['state'] = state[1] if state else 'not_in_current_queue'
        print(json.dumps(item), flush=True)
        folder = item['remote'] + '/results'
        collect_results(item['case'], ['job_identity.txt', 'runtime_source.sha256', 'compiler_version.txt', 'compile.status',
                     'compile.log', 'run.status', 'solver.stdout', 'wrapper.status', 'validation.status',
                     'SimulationSettings2DOpenaccMultilevel.txt', 'FinalNuRe_2DOpenaccMultiblock.dat',
                     'Convergence_2DOpenaccMultiblock.dat', 'NuRe_2DOpenaccMultiblock.dat',
                     'postprocess.status', 'NuReStatistics_2DOpenaccMultiblock.dat'])
    save_json(OUT / 'submission.json', submissions)
    for diagnostic in MANIFEST['diagnostics']:
        collect_results(diagnostic['case'], ['job_identity.txt', 'runtime_source.sha256', 'compiler_version.txt',
                     'compile.status', 'compile.log', 'run.status', 'solver.stdout', 'wrapper.status',
                     'validation.status', 'SimulationSettings2DOpenaccMultilevel.txt', 'FinalNuRe_2DOpenaccMultiblock.dat'])
elif action == 'probe':
    probe('manual_' + datetime.datetime.now().strftime('%H%M%S'))
else:
    raise SystemExit('Unknown action: ' + action)
