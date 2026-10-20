from pathlib import Path
import math,json,re
p=Path('.')
settings=(p/'SimulationSettings2DOpenaccMultiblock.txt').read_text()
assert 'Unsteady target duration completed at a synchronized step.' in settings, 'Missing completion marker'
rows=[s.split() for s in (p/'NuRe_2DOpenaccMultiblock.dat').read_text().splitlines() if s.strip() and not s.startswith('#')]
assert rows, 'Empty history'
previous_step=0
for i,row in enumerate(rows,1):
 assert len(row)==14 and int(row[0])==i, 'Invalid history structure'
 assert int(row[1])>previous_step, 'Nonmonotone history'
 previous_step=int(row[1])
 assert all(math.isfinite(float(v.replace('D','E'))) for v in row[2:]), 'Nonfinite history'
last=rows[-1]; final_step=int(last[1]); final_time=float(last[2])
assert final_time>=1000.0, 'Target duration not reached'
meta=(p/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().splitlines()
assert meta[0].strip()=='MB2D_RESTART 1'
checkpoint=p/meta[1].strip()
assert checkpoint.is_file() and checkpoint.stat().st_size>0, 'Missing final checkpoint'
def after(prefix):
 return meta[next(i for i,s in enumerate(meta) if s.startswith(prefix))+1].split()
state=after('# itc_total NuReSampleCount'); target=after('# requested_target_itc')
assert int(state[0])==int(target[0])==final_step, 'Final checkpoint/history/target step mismatch'
assert int(state[1])==int(last[0]), 'Checkpoint sample count mismatch'
assert abs(float(state[2])-final_time)<1e-8, 'Checkpoint/history time mismatch'
assert all(math.isfinite(float(x)) for x in state[2:]), 'Nonfinite checkpoint metadata'
assert final_time-1000.0<0.51, 'Unexpected target rounding'
report=dict(completed=True,final_step=final_step,final_tff=final_time,samples=len(rows),checkpoint=checkpoint.name,statistical_convergence_assessed=False)
(p/'completion.json').write_text(json.dumps(report,indent=2)+'\n')
(p/'validation.status').write_text('duration_completed_history_finite\n')
print(json.dumps(report))
