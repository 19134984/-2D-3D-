"""P100 transport; authentication is supplied only through the process environment."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
exponent = int(sys.argv[2]) if len(sys.argv) > 2 else 6
assert exponent in (6, 7, 8)
n = 256 if exponent == 6 else 512
LOCAL = ROOT / f'p100/Ra1e{exponent}/sideheated_N{n}_20260918'
CASE = f'/data2/XLLi/PTRT/Ra1e{exponent}/sideheated_N{n}_20260918'
SSH = ['C:/msys64/usr/bin/sshpass.exe', '-e', 'C:/msys64/usr/bin/ssh.exe',
       '-o', 'ConnectTimeout=15', 'xlli@10.169.174.118']


def remote(command, data=None):
    result = subprocess.run(SSH+[command], input=data, capture_output=True, timeout=90)
    if result.returncode:
        print(result.stderr.decode(errors='replace'))
        print(result.stdout.decode(errors='replace'))
        raise SystemExit(result.returncode)
    return result.stdout


action = sys.argv[1]
if action == 'upload':
    manifest = json.loads((LOCAL / 'manifest.json').read_text())
    assert hashlib.sha256((ROOT/'2DRBOpenaccPTRT.F90').read_bytes()).hexdigest() == manifest['local_master_sha256']
    command = f"ssh node07 'test ! -e {CASE} && tar -xf - -C /data2/XLLi/PTRT'"
    print(remote(command, (ROOT/f'p100/ra1e{exponent}_submission.tar').read_bytes()).decode())
    command = f"dos2unix {CASE}/run.pbs && bash -n {CASE}/run.pbs && sha256sum {CASE}/source/solver.F90 {CASE}/source/local_master.F90"
    result = remote(command)
    (LOCAL/'upload_verification.txt').write_bytes(result)
    text = result.decode()
    assert manifest['case_source_sha256'] in text
    assert manifest['local_master_sha256'] in text
    print(text)
elif action == 'submit':
    # An uncertainty or failure must be inspected; never blindly retry qsub.
    dependency = sys.argv[3] if len(sys.argv) > 3 else None
    if dependency:
        assert dependency.endswith('.master') and dependency.split('.')[0].isdigit()
    args = f'-W depend=afterany:{dependency} ' if dependency else ''
    result = remote(f'cd {CASE} && test ! -e job.id && mkdir submission.lock && qsub {args}run.pbs > job.id && cat job.id')
    job = result.decode().strip()
    assert job.endswith('.master') and job.split('.')[0].isdigit(), job
    (LOCAL/'job.id').write_text(job+'\n')
    (LOCAL/'submission.json').write_text(json.dumps(dict(job=job,node='node07',
        remote_case=CASE,dependency=('afterany:'+dependency) if dependency else None),indent=2))
    print(job)
elif action == 'status':
    job = (LOCAL/'job.id').read_text().strip()
    result = remote(f'qstat -f {job}; for f in compile.status run.status validation.status job_identity.txt runtime_source.sha256; do if test -f {CASE}/results/$f; then echo FILE=$f; cat {CASE}/results/$f; fi; done; if test -f {CASE}/results/solver.stdout; then tail -8 {CASE}/results/solver.stdout; fi')
    (LOCAL/'submission_status.txt').write_bytes(result)
    print(result.decode(errors='replace'))
else:
    raise SystemExit('usage: upload | submit | status')
