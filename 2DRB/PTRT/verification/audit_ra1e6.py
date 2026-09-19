"""Independently reconstruct the final uniform-grid checkpoint and compare Table 2."""
from pathlib import Path
import hashlib
import json
import re
import struct
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
exponent = int(sys.argv[1]) if len(sys.argv) > 1 else 6
assert exponent in (6,7)
expected_n = 256 if exponent == 6 else 512
CASE = ROOT/f'p100/Ra1e{exponent}/sideheated_N{expected_n}_20260918'
R = CASE/'results'
meta = dict(line.split(maxsplit=1) for line in (R/'reloadFile2DOpenacc-latest.meta').read_text().splitlines())
n = int(meta['nx'])
assert n == int(meta['ny']) == expected_n
assert meta['flowMode'] == 'steadyFlow'
name = 'reloadFile2DOpenacc-'+meta['reloadFileName']+'.bin'
blob = (R/name).read_bytes()
offset = 0
fields = []
for q in (9, 5):
    size = struct.unpack_from('<i', blob, offset)[0]
    assert size == n*n*q*8
    offset += 4
    a = np.frombuffer(blob, dtype='<f8', count=n*n*q, offset=offset).reshape((n,n,q), order='F')
    assert np.isfinite(a).all()
    offset += size
    assert struct.unpack_from('<i', blob, offset)[0] == size
    offset += 4
    fields.append(a)
assert offset == len(blob)
f,g = fields
ex = np.array([0,1,0,-1,0,1,-1,-1,1])
ey = np.array([0,0,1,0,-1,1,1,-1,-1])
rho = f.sum(axis=2)
assert np.all(rho > 0)
T = g.sum(axis=2)
nu = .1*n*np.sqrt(3*.71/10**exponent)/3
alpha = nu/.71
gbeta = 10**exponent*nu*alpha/n**3
u = (f@ex)/rho
v = (f@ey)/rho + .5*gbeta*T
final = dict(NuVolAvg=float(1+np.mean(u*T)*n/alpha),
             ReVolRMS=float(np.sqrt(np.mean(u*u+v*v))*n/nu),
             Nu_hot=float(np.mean((4-9*T[0,:]+T[1,:])*n/3)),
             Nu_cold=float(np.mean((4+9*T[-1,:]-T[-2,:])*n/3)))
l,r = n//2-1,n//2
final['Nu_middle'] = float(np.mean(n/alpha*.5*(u[l,:]*T[l,:]+u[r,:]*T[r,:])+n*(T[l,:]-T[r,:])))
checks = np.loadtxt(R/'Convergence_PTRT.dat')
assert np.isfinite(checks).all()
assert np.all(checks[-1,1:3] <= 1e-7)
assert checks[-1,0] == int(meta['itc_total'])
# The final diagnostic call intentionally repeats the last converged check.
assert np.array_equal(checks[-1],checks[-2])
assert np.all(np.diff(checks[:-1,0]) == 2000)
assert (R/'run.status').read_text().strip() == '0'
assert (R/'compile.status').read_text().strip() == '0'
source_hash = hashlib.sha256((CASE/'source/solver.F90').read_bytes()).hexdigest()
assert source_hash == (R/'runtime_source.sha256').read_text().split()[0]
stdout = (R/'solver.stdout').read_text()
assert not re.search(r'(?i)(?<![a-z])(nan|inf|infinity)(?![a-z])|CUDA_ERROR|ERROR STOP|FATAL|Error:', stdout)
reported = {}
for key, label in [('NuVolAvg','NuVolAvg'),('ReVolRMS','ReVolAvg'),('Nu_hot','Nu_hot'),('Nu_cold','Nu_cold'),('Nu_middle','Nu_middle')]:
    reported[key] = float(re.findall(r'^'+label+r'\s*=\s*(\S+)',stdout,re.M)[-1])
    assert abs(final[key]-reported[key]) < 1e-8, (key,final[key],reported[key])
refs = {'paper_257':dict(NuVolAvg=8.8282,Nu_hot=8.8509,ReVolRMS=99.6110),
        'paper_513':dict(NuVolAvg=8.8202,Nu_hot=8.8317,ReVolRMS=99.3839),
        'paper_1025':dict(NuVolAvg=8.8188,Nu_hot=8.8246,ReVolRMS=99.3576),
        'paper_extrapolated':dict(NuVolAvg=8.8186,Nu_hot=8.8206,ReVolRMS=99.3542)}
for key,value in zip(refs, (8.8277,8.8200,8.8188,8.8186)):
    refs[key]['Nu_middle'] = value
if exponent == 7:
    refs = {'paper_513':dict(NuVolAvg=16.5245,Nu_hot=16.5461,Nu_middle=16.5243,ReVolRMS=230.6601),
            'paper_1025':dict(NuVolAvg=16.5133,Nu_hot=16.5242,Nu_middle=16.5133,ReVolRMS=230.2448),
            'paper_2049':dict(NuVolAvg=16.5114,Nu_hot=16.5169,Nu_middle=16.5114,ReVolRMS=230.1948),
            'paper_extrapolated':dict(NuVolAvg=16.5110,Nu_hot=16.5132,Nu_middle=16.5110,ReVolRMS=230.1879)}
report = dict(final_checkpoint_metrics=final, reported=reported, references=refs,
    deviations_percent={name:{k:100*(final[k]-v)/v for k,v in row.items()} for name,row in refs.items()},
    step=int(meta['itc_total']),t_ff=float(meta['time_tf']),errors=checks[-1,1:3].tolist(),
    checkpoint_f_g_finite=True, checkpoint_sha256=hashlib.sha256(blob).hexdigest(),
    source_sha256=source_hash, mass_drift_percent=float(100*(rho.mean()-1)),
    hot_cold_difference_percent=100*abs(final['Nu_hot']-final['Nu_cold'])/final['Nu_hot'])
(CASE/'文献对比.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2))
