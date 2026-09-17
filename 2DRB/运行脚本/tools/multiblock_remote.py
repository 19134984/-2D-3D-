"""Transport prepared test files through the P100 login node; password via SSHPASS only."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'运行脚本/Multiblock_sideheated_20260916_v3'
SSH = ['C:/msys64/usr/bin/sshpass.exe','-e','C:/msys64/usr/bin/ssh.exe',
       '-o','ConnectTimeout=15','xlli@10.169.174.118']
action=sys.argv[1]
if action=='upload':
    command="ssh node05 'test ! -e /data2/XLLi/Multiblock/campaign_20260916_v3/run.pbs && tar -xf - -C /data2/XLLi/Multiblock'"
    subprocess.run(SSH+[command],input=(OUT/'cases.tar').read_bytes(),check=True)
elif action=='compile':
    command="ssh node05 'bash /data2/XLLi/Multiblock/Ra1e6/ratio2/sideheated_N256_20260916_v3/run.sh compile-only && bash /data2/XLLi/Multiblock/Ra1e7/ratio2/sideheated_N512_20260916_v3/run.sh compile-only'"
    subprocess.run(SSH+[command],check=True)
elif action=='submit':
    r=subprocess.run(SSH+["cd /data2/XLLi/Multiblock/campaign_20260916_v3 && test ! -s job.id && qsub run.pbs"],capture_output=True,check=True)
    job=r.stdout.decode().strip()
    assert job.endswith('.master') and job.split('.')[0].isdigit(),job
    (OUT/'campaign_20260916_v3/job.id').write_text(job)
    subprocess.run(SSH+[f"printf '%s\\n' '{job}' > /data2/XLLi/Multiblock/campaign_20260916_v3/job.id"],check=True)
    print(job)
elif action=='inspect':
    subprocess.run(SSH+["ssh node05 'for c in /data2/XLLi/Multiblock/Ra*/ratio2/*; do echo CASE=$c; cat $c/results/compile.status; tail -15 $c/results/compile.log; ls $c/results; tail -8 $c/results/solver.stdout 2>/dev/null; done'; qstat -u xlli"],check=True)
elif action=='collect-first':
    remote='/data2/XLLi/Multiblock/Ra1e6/ratio2/sideheated_N256_20260916_v3/results'
    local=OUT/'Ra1e6/ratio2/sideheated_N256_20260916_v3/results'
    for name in ['compile.status','compile.log','runtime_source.sha256','run.status',
                 'solver.stdout','job_identity.txt','timing.txt','validation.status',
                 'SimulationSettings2DOpenaccMultiblock.txt','pbs.stdout','pbs.stderr',
                 'Convergence_2DOpenaccMultiblock.dat','NuRe_2DOpenaccMultiblock.dat',
                 'reloadFile2DOpenaccMultiblock-latest.meta','run.started','run.finished']:
        result=subprocess.run(SSH+[f'cat {remote}/{name}'],capture_output=True)
        if result.returncode==0:
            (local/name).write_bytes(result.stdout)
    job=(OUT/'campaign_20260916_v3/job.id').read_text().strip()
    result=subprocess.run(SSH+['qstat -f '+job],capture_output=True)
    (local/'qstat.txt').write_bytes(result.stdout)
    print(result.stdout.decode(errors='replace'))
    print((local/'solver.stdout').read_text()[-4000:])
else:
    raise SystemExit('unknown action')
