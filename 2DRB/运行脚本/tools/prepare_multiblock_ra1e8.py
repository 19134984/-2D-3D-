"""Parameterize the latest local multiblock master for the authorized Ra=1e8 run."""
from pathlib import Path
import hashlib,json,re

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'运行脚本/Multiblock_sideheated_20260916_v3'
OLD=BASE/'Ra1e7/ratio2/sideheated_N512_20260916_v3'
REL='Ra1e8/ratio2/sideheated_N512_20260916_v3'
CASE=BASE/REL
raw=(ROOT/'多块网格/2DRBOpenaccMultiblock.F90').read_bytes()
s=raw.decode('utf-8-sig').replace('\r\n','\n')
assert 'EnableUseG' not in s
for macro,on in [('steadyFlow',True),('unsteadyFlow',False),('RayleighBenardCell',False),
                 ('HorizontalWallsConstT',False),('VerticalWallsAdiabatic',False),
                 ('SideHeatedCell',True),('HorizontalWallsAdiabatic',True),('VerticalWallsConstT',True)]:
    s,count=re.subn(r'^!?#define '+macro+r'[ \t]*$',('#define ' if on else '!#define ')+macro,s,flags=re.M)
    assert count==1,(macro,count)
changes={'nx':'512','ny':'512','refineRatio':'2','Rayleigh':'1.0d8','Prandtl':'0.71d0',
         'fineLayerCellsLeft':'64','fineLayerCellsRight':'65','fineLayerCellsBottom':'64','fineLayerCellsTop':'65'}
for key,val in changes.items():
    s,count=re.subn(r'(\b'+key+r'\s*=\s*)[\d.]+(?:[dDeE][+-]?\d+)?',lambda m:m[1]+val,s,count=1)
    assert count==1,key
assert re.search(r'loadInitField\s*=\s*0',s)
CASE.mkdir(parents=True,exist_ok=False)
for child in ['source','results']:
    (CASE/child).mkdir()
data=s.replace('\r\n','\n').encode('utf-8')
(CASE/'source/solver.F90').write_bytes(data)
(CASE/'source/local_master.F90').write_bytes(raw)
meta=json.loads((OLD/'manifest.json').read_text())
old_hash=meta['case_source_sha256']
meta.update(case=REL,ra=100000000,local_master_sha256=hashlib.sha256(raw).hexdigest(),
            case_source_sha256=hashlib.sha256(data).hexdigest(),source_changes=changes,
            reference_grid='Paper Table 2 minimum 513^2; user convention fine-equivalent 512^2')
(CASE/'manifest.json').write_text(json.dumps(meta,indent=2),encoding='utf-8',newline='\n')
run=(OLD/'run.sh').read_text().replace('Ra1e7/','Ra1e8/').replace(old_hash,meta['case_source_sha256'])
(CASE/'run.sh').write_text(run,encoding='utf-8',newline='\n')
helper=(ROOT/'运行脚本/tools/multiblock_ra1e7.py').read_text(encoding='utf-8')
helper=helper.replace('Ra=1e7','Ra=1e8').replace('Ra1e7/','Ra1e8/').replace('MB7N512R2','MB8N512R2')
helper=helper.replace("command=f'test ! -e {REMOTE}/job.id", "command=f'mkdir -p {REMOTE}/results && test ! -e {REMOTE}/job.id")
(ROOT/'运行脚本/tools/multiblock_ra1e8.py').write_text(helper,encoding='utf-8',newline='\n')
print(json.dumps(meta,indent=2))
