from pathlib import Path
import tarfile, json, re, math, hashlib
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'server_evidence'
OUT.mkdir(exist_ok=True)
with tarfile.open(ROOT/'all_tests_review_20261023.tar.gz') as tf:
    for m in tf.getmembers():
        p=(OUT/m.name).resolve()
        assert OUT.resolve() in p.parents and m.isfile()
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes(tf.extractfile(m).read())
with tarfile.open(ROOT/'completion_reviews_20261023.tar.gz') as tf:
    for m in tf.getmembers():
        p=(OUT/'RB'/m.name).resolve()
        assert OUT.resolve() in p.parents and m.isfile()
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes(tf.extractfile(m).read())

def interp(rows,t,k,square=False):
    for a,b in zip(rows,rows[1:]):
        if a[2]<=t<=b[2]:
            va=a[k]**2 if square else a[k];vb=b[k]**2 if square else b[k]
            return va+(vb-va)*(t-a[2])/(b[2]-a[2])
    raise ValueError(t)

def avg(rows,lo,hi,k,square=False):
    pts=[(lo,interp(rows,lo,k,square))]+[(r[2],r[k]**2 if square else r[k]) for r in rows if lo<r[2]<hi]+[(hi,interp(rows,hi,k,square))]
    return math.fsum((b[0]-a[0])*(a[1]+b[1])/2 for a,b in zip(pts,pts[1:]))/(hi-lo)

res={}
for f in OUT.rglob('NuRe_2DOpenaccMultiblock.dat'):
    d=f.parent
    rows=[list(map(float,l.split())) for l in f.read_text().splitlines() if l.strip() and not l.lstrip().startswith('#')]
    assert all(len(r)==14 and all(math.isfinite(x) for x in r) for r in rows)
    assert all(int(b[0])==int(a[0])+1 and b[1]>a[1] and b[2]>a[2] for a,b in zip(rows,rows[1:]))
    src=d.parent/'source/solver.F90'
    if not src.exists():
        choices=list(d.parent.parent.glob('source/*.F90'))
        src=choices[0] if choices else None
    s=src.read_text(encoding='utf-8') if src else ''
    settings=(d/'SimulationSettings2DOpenaccMultiblock.txt').read_text() if (d/'SimulationSettings2DOpenaccMultiblock.txt').exists() else ''
    n=re.search(r'\bnx\s*=\s*(\d+)',s,re.I)
    N=int(n[1]) if n else None
    r={'N':N,'rows':len(rows),'last_tff':rows[-1][2],'status':(d/'run.status').read_text().strip() if (d/'run.status').exists() else 'running',
       'all_saved_values_finite':True,'history_contiguous':True}
    if src and (d/'runtime_source.sha256').exists():
        h=hashlib.sha256(src.read_bytes()).hexdigest()
        assert h in (d/'runtime_source.sha256').read_text()
        r['runtime_source_verified']=True
    if (d/'completion.json').exists():
        completion=json.loads((d/'completion.json').read_text())
        if completion.get('completed_target'):
            assert completion['final_checkpoint']['all_values_finite']
            assert completion['final_fields']['all_nodes_finite']
            r['server_final_checkpoint_and_fields_finite']=True
    if N:
        r['mass_change_percent']=100*(rows[-1][8]/N**2-1)
        r['max_mass_drift_percent']=100*max(abs(x[8]/N**2-1) for x in rows)
    text=(d/'solver.stdout').read_text(errors='replace')
    fail=re.search(r'Invalid state: time, xOffset,yOffset, i,j:\s*(\d+)',text)
    if fail:r['invalid_detected_tff']=int(fail[1])*rows[-1][2]/rows[-1][1]
    lo,hi=(500,1000) if rows[-1][2]>=1000 else (50,100)
    if rows[-1][2]>=hi:
        nu=avg(rows,lo,hi,3);re_st=math.sqrt(avg(rows,lo,hi,4,True))
        hot=avg(rows,lo,hi,5);cold=avg(rows,lo,hi,6)
        refnu,refre=(53.51,15101) if 'Ra1e9' in str(f) else (25.25,3662)
        r['window']=[lo,hi]
        r.update(Nu=nu,Re_space_time_RMS=re_st,Nu_hot=hot,Nu_cold=cold,
                 Nu_error_percent=100*abs(nu/refnu-1),Re_error_percent=100*abs(re_st/refre-1),
                 wall_volume_difference_percent=100*max(abs(hot/nu-1),abs(cold/nu-1)),
                 Nu_half_difference_percent=100*abs(avg(rows,lo,(lo+hi)/2,3)-avg(rows,(lo+hi)/2,hi,3))/abs(nu),
                 Re_half_difference_percent=100*abs(math.sqrt(avg(rows,lo,(lo+hi)/2,4,True))-math.sqrt(avg(rows,(lo+hi)/2,hi,4,True)))/re_st)
        if (d/'completion.json').exists() and completion.get('completed_target'):
            assert math.isclose(nu,completion['statistics']['whole'][0],rel_tol=1e-12)
            assert math.isclose(re_st,completion['statistics']['whole'][1],rel_tol=1e-12)
    res[str(d.relative_to(OUT))]=r
(ROOT/'independent_statistics.json').write_text(json.dumps(res,indent=2)+'\n')
for name,r in res.items():print(name,json.dumps(r))
