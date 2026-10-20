from pathlib import Path
import subprocess,shutil,hashlib,datetime,json,math,tarfile
root=Path('/data2/XLLi/Multiblock/data2D/RB')
out=root/'postprocess_20261020_0905'
src=out/'postprocess.F90'
audit={'checked':datetime.datetime.now().isoformat(),'job_id':'7266.master','postprocess_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'window_tff':[500,1000],'compiler':'nvfortran 24.3 CPU executable','cases':[]}
for e in range(6,11):
 name='Ra1e'+str(e); r=root/name/'results'; d=out/name; d.mkdir()
 for filename in ['NuRe_2DOpenaccMultiblock.dat','solver.stdout','SimulationSettings2DOpenaccMultiblock.txt','reloadFile2DOpenaccMultiblock-latest.meta','run.started','run.finished','compile.status','run.status','validation.status','completion.json','completion_validation.log','runtime_source.sha256']:
  if (r/filename).exists(): shutil.copy2(str(r/filename),str(d/filename))
 rows=[x.split() for x in (d/'NuRe_2DOpenaccMultiblock.dat').read_text().splitlines() if x.strip() and not x.startswith('#')]
 info={'case':name,'samples':len(rows),'last_history_tff':float(rows[-1][2]),'solver_exit':int((d/'run.status').read_text()),'nonfinite_history_rows':sum(not all(math.isfinite(float(v)) for v in x[2:]) for x in rows),'history_sha256':hashlib.sha256((d/'NuRe_2DOpenaccMultiblock.dat').read_bytes()).hexdigest()}
 with (d/'postprocess.stdout').open('w') as log:
  run=subprocess.run([str(out/'postprocess.exe')],cwd=str(d),stdout=log,stderr=subprocess.STDOUT)
 info['postprocess_exit']=run.returncode
 f=d/'NuReStatistics_2DOpenaccMultiblock.dat'
 if f.exists(): info['statistics']=f.read_text()
 audit['cases'].append(info)
 print(json.dumps(info))
for name in ['case_exit_status.tsv','submission.json','submission_verification.json','run.status','campaign.started','campaign.finished']:
 f=root/'serial_node07'/name
 if f.exists():
  (out/'serial_node07').mkdir(exist_ok=True); shutil.copy2(str(f),str(out/'serial_node07'/name))
(out/'audit.json').write_text(json.dumps(audit,indent=2))
with tarfile.open(str(out)+'.tar.gz','w:gz') as tf: tf.add(str(out),arcname=out.name)
print('POSTPROCESS_BUNDLE',str(out)+'.tar.gz')
