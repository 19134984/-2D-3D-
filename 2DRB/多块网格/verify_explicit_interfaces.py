"""Verify explicit coarse/fine interfaces, three histories, diagnostics and restarts.
Usage: python verify_explicit_interfaces.py PATH_TO_PRE_CHANGE_SOURCE
Builds and numerical outputs stay in a temporary directory; no solver settings are edited.
Uses gfortran OpenACC host with bounds checks, not a GPU validation.
"""
from pathlib import Path
import sys,re,json,hashlib,argparse
import numpy as np
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here))
import verify_multiblock as v
import verify_compact as c
parser=argparse.ArgumentParser(description='Compare explicit interfaces against the previous shared-buffer source.')
parser.add_argument('baseline',type=Path)
args=parser.parse_args()
old=args.baseline.read_text(encoding='utf-8-sig');new=v.SOURCE.read_text(encoding='utf-8-sig')

def macro(s,key,enabled):
    return re.sub(r'^!?#define '+key+r'\s*$',('#define ' if enabled else '!#define ')+key,s,flags=re.M)

def config(source,ratio,side,steady,ha,end,restart=False):
    source=v.variant(source,ratio=ratio,restart=restart)
    for key,val in [('steadyFlow',steady),('unsteadyFlow',not steady),('SideHeatedCell',side),
                    ('RayleighBenardCell',not side),('HorizontalWallsConstT',not side),
                    ('VerticalWallsAdiabatic',not side),('HorizontalWallsAdiabatic',side),
                    ('VerticalWallsConstT',side),('SideHeatedHa',ha)]: source=macro(source,key,val)
    source=re.sub(r'(::\s*itc_max\s*=\s*ceiling\()\d+\.0d0',rf'\g<1>{end}.0d0',source)
    source=re.sub(r'(::\s*itc_max\s*=)\s*\d+',rf'\g<1> {end}',source)
    source=re.sub(r'(::\s*unsteadyRunDuration\s*=)\s*[^\n!]+',rf'\g<1> {end}.0d0/timeUnit ',source)
    source=re.sub(r'(\beps[UT]\s*=)\s*[^,\n!]+',r'\g<1> 0.0d0',source)
    for name,steps in [('outputSnapshotInterval',5.5),('reloadFileInterval',9.5),('outputPltFileInterval',7.5)]:
        source=re.sub(r'(::\s*'+name+r'\s*=)\s*[^\n!]+',rf'\g<1> {max(steps,float(ratio))}d0/timeUnit ',source)
    source=re.sub(r'(^[ \t]*steadyCheckIntervalItc[ \t]*=[ \t]*)[^\n]+',rf'\g<1> {((10+ratio-1)//ratio)*ratio}',source,flags=re.M)
    dump="    call update_host_all(.true.)\n    open(unit=79,file='all_fields.bin',form='unformatted',access='stream',status='replace')\n"
    for r in ['coarse','left','right','bottom','top']:
        if r=='left':dump+='    if(refineRatio>1) then\n'
        dump+='    write(79) '+','.join(f'{a}_{r}' for a in ['f','g','u','v','T','rho','Fx','Fy'])+'\n'
    dump+="    endif\n    close(79)\n    open(unit=79,file='all_history.bin',form='unformatted',access='stream',status='replace')\n"
    if 'coarseHistoryPrevious' in source:
        dump+='    write(79) coarseHistoryPrevious,coarseHistoryCurrent,coarseHistoryNext\n'
    else:
        dump+='    write(79) coarseHistory(:,0,:)\n    if(refineRatio>1) write(79) coarseHistory(:,1,:),coarseHistory(:,2,:)\n'
    dump+='    close(79)\n'
    source=source.replace('    call exit_data_2d_openacc()',dump+'    call exit_data_2d_openacc()',1)
    return source

checks=[]
cases=[('r1_rb',1,False,True,False,(16,17,16,17)),('r2_rb',2,False,False,False,(16,17,16,17)),
('r2_side',2,True,True,False,(16,17,16,17)),('r4_side',4,True,False,False,(16,17,16,17)),
('r8_rb',8,False,True,False,(32,33,32,33)),('r2_ha',2,True,True,True,(16,17,16,17)),
('r2_asymmetric',2,True,False,False,(44,13,12,13))]
print('BUILD',v.BUILD,flush=True)
for label,ratio,side,steady,ha,walls in cases:
    dims=dict(n=128,ny=112,ra=10000,walls=walls) if ratio==8 else dict(n=96,ny=80,ra=10000,walls=walls)
    def build(src,tag,end,restart=False):
        return v.compile_source(label+tag,config(src,ratio,side,steady,ha,end,restart),**dims)
    before,exe=build(old,'_before',40);v.run([exe],before)
    after,exe=build(new,'_after',40);v.run([exe],after)
    for file in ['all_fields.bin']+(['all_history.bin'] if ratio>1 else []):
        x=np.fromfile(before/file,dtype='<f8');y=np.fromfile(after/file,dtype='<f8')
        assert np.isfinite(y).all()
        assert np.array_equal(x,y),(label,file,float(np.max(np.abs(x-y))))
    assert (before/'NuRe_2DOpenaccMultiblock.dat').read_bytes()==(after/'NuRe_2DOpenaccMultiblock.dat').read_bytes()
    for path in after.glob('*Snapshot-*.bin'):c.validate_snapshot(path)
    half=(22//ratio)*ratio
    split,exe=build(new,'_split',half);v.run([exe],split)
    _,exe=build(new,'_resume',40-half if steady else 40,True);v.run([exe],split)
    assert c.checkpoint(after).read_bytes()==c.checkpoint(split).read_bytes(),label+' restart checkpoint'
    for file in ['all_fields.bin','all_history.bin','NuRe_2DOpenaccMultiblock.dat']:
        assert (after/file).read_bytes()==(split/file).read_bytes(),(label,file,'restart')
    checks.append(dict(case=label,fields_exact=True,history_exact=ratio>1,NuRe_exact=True,restart_exact=True))
    print('PASS',label,'fields, history, diagnostics, current-format restart',flush=True)
report=dict(source_sha256=hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),device='gfortran OpenACC host',build_directory=str(v.BUILD),checks=checks)
(here/'explicit_interfaces_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('REPORT',here/'explicit_interfaces_verification.json',flush=True)
