"""Read-only result assessment and the Ra9 -> Ra10 release gate; Python stdlib."""
from pathlib import Path
import json
import math
import statistics
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
cfg=json.loads((ROOT/'campaign.json').read_text())
name=sys.argv[1]
if len(sys.argv)>2 and sys.argv[2]=='--gate':
    evidence=ROOT/name/'results/assessment.json'
    report=json.loads(evidence.read_text()) if evidence.exists() else {}
    print(json.dumps({'Ra9_case':name,'release_Ra10':report.get('release_Ra10',False)}))
    sys.exit(0 if report.get('release_Ra10',False) else 91)

case=cfg['cases'][name]; folder=ROOT/name/'results'
result={'case':name,'target_tff':case['target_tff'],'structural_pass':False,
        'release_Ra10':False,'grid_independence_verified':False}
errors=[]


def rows(pattern):
    paths=list(folder.glob(pattern))
    assert len(paths)==1,pattern
    data=[]
    for line in paths[0].read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'): continue
        row=[float(x.replace('D','E').replace('d','e')) for x in line.split()]
        assert all(math.isfinite(x) for x in row),str(paths[0])+' contains nonfinite data'
        data.append(row)
    assert data,pattern+' empty'
    return data


try:
    nu=rows('Nu_VolAvg_*.dat'); re=rows('Re_VolAvg_*.dat')
    walls=rows('PTRTWallHistory.dat'); local=rows('PTRTCollisionHistory.dat')
    assert len(nu)==len(re)==len(walls)==1001,'Expected complete t=0..1000 histories'
    assert all(abs(a[0]-i)<1e-6 and abs(b[0]-i)<1e-6 and abs(w[0]-i)<1e-4
               for i,(a,b,w) in enumerate(zip(nu,re,walls))),'Mismatched history clocks'
    assert int((folder/'solver.status').read_text())==0,'Nonzero solver exit'
    settings=(folder/'SimulationSettings2DOpenaccLBMCDE_D2Q5.txt').read_text()
    assert 'Flow collision = Luo P-TRT' in settings,'P-TRT runtime identity missing'
    meta=list(folder.glob('*-latest.meta')); assert len(meta)==1
    ledger=dict(line.split(None,1) for line in meta[0].read_text().splitlines() if line.strip())
    time_final=float(ledger['time_tf'])
    tol=2/(case['nx']*math.sqrt(3)/case['Ma'])
    assert abs(time_final-1000)<tol,'Stopping checkpoint did not reach target'
    maxmass=max(x[2] for x in local); maxforce=max(x[3] for x in local); maxghost=max(x[4] for x in local)
    assert max(maxmass,maxforce,maxghost)<1e-10,'Collision invariants exceed tolerance'
    assert all(x[-1]==0 for x in local),'Nonfinite populations detected'
    result.update(finite_histories=True,target_reached=True,last_checkpoint_tff=time_final,
        last_checkpoint_step=int(ledger['itc_total']),samples=len(nu),collision_residual_max={
            'mass':maxmass,'force':maxforce,'ghost':maxghost})
    a=[(x[1],y[1]) for x,y in zip(nu,re) if 500<=x[0]<750]
    b=[(x[1],y[1]) for x,y in zip(nu,re) if 750<=x[0]<=1000]
    def avg(data): return statistics.mean(x[0] for x in data),math.sqrt(statistics.mean(x[1]**2 for x in data))
    n1,r1=avg(a); n2,r2=avg(b)
    dn=100*abs(n1-n2)/abs(n2); dr=100*abs(r1-r2)/abs(r2)
    wall=tuple(statistics.mean(row[i] for row in walls if 750<=row[0]<=1000) for i in [1,2])
    wallrel=max(abs(wall[0]-n2),abs(wall[1]-n2),abs(wall[0]-wall[1]))/abs(n2)
    policy=cfg['release_policy']
    result.update(Nu=n2,Re_rms=r2,Nu_half_window_percent=dn,Re_half_window_percent=dr,
        statistics_converged=dn<policy['half_window_percent'] and dr<policy['half_window_percent'],
        Nu_bottom_top=wall,wall_budget_max_relative=wallrel,
        Nu_reference_relative=(n2-case['Nu_reference'])/case['Nu_reference'],
        Re_reference_relative=(r2-case['Re_reference'])/case['Re_reference'])
    assert n2>1 and r2>0,'No physically reasonable convection in final window'
except (AssertionError,OSError,ValueError,KeyError,IndexError,ZeroDivisionError) as exc:
    errors.append(str(exc))

audit=subprocess.run([sys.executable,str(ROOT/'audit_flow_history.py'),'--case',str(folder),'--strict'],
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
(folder/'history_audit.txt').write_text(audit.stdout)
if audit.returncode: errors.append('Complete-history/restart-ledger audit failed')
result['structural_errors']=errors
result['structural_pass']=not errors
if result['structural_pass']:
    p=cfg['release_policy']
    result['release_Ra10']=bool(result['statistics_converged'] and result['wall_budget_max_relative']<p['wall_budget_relative']
        and abs(result['Nu_reference_relative'])<p['Nu_reference_relative']
        and abs(result['Re_reference_relative'])<p['Re_reference_relative'])
result['note']='Numerical horizon, finite values, statistical convergence, reference agreement and grid independence are separate claims.'
(folder/'assessment.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(result,indent=2,allow_nan=False))
sys.exit(0 if result['structural_pass'] else 90)
