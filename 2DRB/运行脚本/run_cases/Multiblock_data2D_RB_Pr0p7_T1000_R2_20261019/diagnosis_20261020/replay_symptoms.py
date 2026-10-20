from pathlib import Path
import json,sys
root=Path(__file__).resolve().parents[1]/'postprocess_20261020_0905'
a=json.loads((root/'verified_audit.json').read_text())
red=[]
for c in a['cases']:
 if c['case']=='Ra1e6' and max(c['relative_error_Nu_percent'],c['relative_error_Re_percent'])>5:
  red.append(c['case']+': benchmark mismatch Nu='+str(c['relative_error_Nu_percent'])+'%, Re='+str(c['relative_error_Re_percent'])+'%')
 if c['case'] in ['Ra1e8','Ra1e9','Ra1e10']:
  log=(root/c['case']/'solver.stdout').read_text()
  if c['solver_exit']!=0 and 'ERROR STOP Nonfinite or nonpositive density in owned cells' in log:
   red.append(c['case']+': reproduced saved failure at step '+str(c['failure_step']))
print('\n'.join('RED '+s for s in red))
sys.exit(1 if red else 0)
