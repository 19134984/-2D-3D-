"""Verify direct interface-moment collision, fields, diagnostics and exact new-format restarts.
Usage: python verify_interface_moments.py PATH_TO_PRE_CHANGE_SOURCE [--candidate SOURCE] [--report REPORT.json]
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
parser=argparse.ArgumentParser(description='Compare direct interface-moment collision against the preceding source.')
parser.add_argument('baseline',type=Path)
parser.add_argument('--candidate',type=Path,default=v.SOURCE)
parser.add_argument('--report',type=Path,default=here/'interface_moments_verification.json')
args=parser.parse_args()
old=args.baseline.read_text(encoding='utf-8-sig');new=args.candidate.read_text(encoding='utf-8-sig')

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


def close_values(x,y,label,rtol=2.0e-11,atol=2.0e-12):
    assert x.shape==y.shape,(label,x.shape,y.shape)
    assert np.isfinite(x).all() and np.isfinite(y).all(),(label,'nonfinite')
    maximum=float(np.max(np.abs(x-y))) if x.size else 0.0
    assert np.allclose(x,y,rtol=rtol,atol=atol),(label,maximum)
    return maximum

def numeric_text(path):
    rows=[]
    for line in path.read_text().splitlines():
        try: row=[float(x.replace('D','E')) for x in line.split()]
        except ValueError: continue
        if row: rows.extend(row)
    return np.asarray(rows)

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
    before,exe=build(old,'_before',39);v.run([exe],before)
    after,exe=build(new,'_after',39);v.run([exe],after)
    errors={}
    for file in ['all_fields.bin']+['all_history.bin']:
        x=np.fromfile(before/file,dtype='<f8');y=np.fromfile(after/file,dtype='<f8')
        assert np.isfinite(y).all()
        errors[file]=close_values(x,y,(label,file))
    errors['NuRe']=close_values(numeric_text(before/'NuRe_2DOpenaccMultiblock.dat'),numeric_text(after/'NuRe_2DOpenaccMultiblock.dat'),(label,'NuRe'),rtol=2e-10,atol=2e-10)
    old_checkpoint=c.checkpoint(before).read_bytes()
    new_checkpoint=c.checkpoint(after).read_bytes()
    assert len(old_checkpoint)==len(new_checkpoint),label+' checkpoint size'
    # Exact checkpoint validity and population storage are verified below by continuous/restart equivalence.

    for path in after.glob('*Snapshot-*.bin'):c.validate_snapshot(path)
    compared=[]
    for pattern in ['*Snapshot-*.bin','*.dat']:
        names=sorted(path.name for path in before.glob(pattern))
        assert names==sorted(path.name for path in after.glob(pattern)),(label,pattern,'file names')
        for name in names:
            if name.endswith('.dat'):
                close_values(numeric_text(before/name),numeric_text(after/name),(label,name),rtol=2e-9,atol=2e-9)
            compared.append(name)

    half=(22//ratio)*ratio
    split,exe=build(new,'_split',half);v.run([exe],split)
    _,exe=build(new,'_resume',39-half if steady else 39,True);v.run([exe],split)
    assert c.checkpoint(after).read_bytes()==c.checkpoint(split).read_bytes(),label+' restart checkpoint'
    for file in ['all_fields.bin','all_history.bin','NuRe_2DOpenaccMultiblock.dat']:
        assert (after/file).read_bytes()==(split/file).read_bytes(),(label,file,'restart')
    checks.append(dict(case=label,max_absolute_errors=errors,restart_exact=True,output_files_checked=compared))
    print('PASS',label,errors,'continuous/restart exact',flush=True)
report=dict(source_sha256=hashlib.sha256(args.candidate.read_bytes()).hexdigest(),baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),device='gfortran OpenACC host',build_directory=str(v.BUILD),checks=checks)
args.report.write_text(json.dumps(report,indent=2),encoding='utf-8')
print('REPORT',args.report,flush=True)
