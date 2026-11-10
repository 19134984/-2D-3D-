"""Accept the frozen mixed-ratio GPU smoke tests and their CPU comparison."""
from pathlib import Path
import datetime, hashlib, json, math, re, sys

root=Path(__file__).resolve().parents[2]
campaign=root/'运行脚本/run_cases'/sys.argv[1]
manifest=json.loads((campaign/'manifest.json').read_text(encoding='utf-8'))
cpu=json.loads((campaign/'cpu_preflight_validation.json').read_text(encoding='utf-8'))
assert cpu['status']=='passed' and cpu['compile_returncode']==cpu['run_returncode']==0
assert cpu['source_sha256']==cpu['copied_source_sha256']
checks=[]
for case in manifest['diagnostics']:
    folder=campaign/case['case']; results=folder/'results'
    source_hash=hashlib.sha256((folder/'source/solver.F90').read_bytes()).hexdigest()
    assert source_hash==case['case_source_sha256']
    runtime_hash=(results/'runtime_source.sha256').read_text().split()[0]
    assert runtime_hash==source_hash
    for name in ['compile.status','run.status','wrapper.status']:
        assert (results/name).read_text().strip()=='0', (case['case'],name)
    assert (results/'validation.status').read_text().strip()=='finite_run_completed'
    identity=(results/'job_identity.txt').read_text()
    assert re.search(r'^host=node05$',identity,re.M),identity
    log=(results/'solver.stdout').read_text()
    steps=case['expected_fine_steps']
    assert f'Completed total fine steps: {steps}; levels: 3' in log
    final_text=(results/'FinalNuRe_2DOpenaccMultiblock.dat').read_text()
    assert f'# itc_total = {steps}\n' in final_text
    values=[float(x) for line in final_text.splitlines() if line and not line.startswith('#') for x in line.split()]
    assert len(values)==12 and all(math.isfinite(v) for v in values)
    settings=(results/'SimulationSettings2DOpenaccMultilevelMixedRatio.txt').read_text()
    assert 'Interface refinement ratios (fine to coarse): 4 2' in settings
    check=dict(case=case['case'],source_sha256=source_hash,runtime_source_sha256=runtime_hash,steps=steps,all_final_metrics_finite=True,FinalNuRe_values=values)
    if case['case']=='preflight/small':
        assert source_hash==cpu['source_sha256']
        ref=cpu['FinalNuRe_values']
        assert len(ref)==len(values)
        differences=[abs(a-b) for a,b in zip(values,ref)]
        scaled=[d/max(1.0,abs(b)) for d,b in zip(differences,ref)]
        assert max(scaled)<1e-11,scaled
        check['cpu_comparison']=dict(max_absolute_difference=max(differences),max_difference_scaled_by_max_1_abs_reference=max(scaled),acceptance_tolerance=1e-11,per_metric_absolute_difference=differences)
    else:
        assert 'Configuration: 4096 4096 3 ' in settings
        assert 'Level 1 dx,Snu,Sq,Qk,Qnu,area,tau:  1.' in settings
    checks.append(check)
proof=dict(status='passed',timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),master_sha256=manifest['master_sha256'],backend='NVHPC 24.3 OpenACC Tesla P100, node05',checks=checks,scope='64 fine steps: implementation and initial numerical sanity only; not long-time Ra1e9 stability, benchmark agreement, or grid independence')
(campaign/'gpu_preflight_validation.json').write_text(json.dumps(proof,indent=2)+'\n',encoding='utf-8')
print(json.dumps(proof))
