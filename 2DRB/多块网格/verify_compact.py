"""Compare five named rectangles with an explicitly supplied connected-ring source.

Runs small OpenACC-host cases with bounds checks. Build files stay in a temp folder.
Usage: python verify_compact.py PATH_TO_PRE_CHANGE_SOURCE
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import struct
import numpy as np
import verify_multiblock as v
from convert_restart_v10 import convert

HISTORY = ['rhoHistory', 'uHistory', 'vHistory', 'THistory', 'FxHistory', 'FyHistory',
           'flowNeqHistory', 'thermalNeqHistory']


def node_record(region, i='i', j='j', nh='0'):
    fields = [f'f_{region}({i},{j},:)', f'g_{region}({i},{j},:)']
    fields += [f'{n}_{region}({i},{j})' for n in ['rho', 'u', 'v', 'T', 'Fx', 'Fy']]
    history = [f'{n}_{region}({i},{j},'+(':'+',' if n in HISTORY[-2:] else '')+'k)' for n in HISTORY]
    return 'write(77) '+', '.join(fields)+', ('+', '.join(history)+f', k=0,{nh})\n'


def driver(compact, steps=40):
    ni, nj, nh = ('nxCoarse', 'nyCoarse', 'historyLastCoarse') if compact else (
        'xLocalCount(1)', 'yLocalCount(1)', 'historyLast(1)')
    body = f'''program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,k,ii,jj,step
    call acc_init(acc_device_host)
    call initial()
    call enter_data_2d_openacc()
    do step=1,{steps}/refineRatio
        call advance_multiblock()
    enddo
    call update_host_all(.true.)
    open(unit=77,file='canonical.bin',form='unformatted',access='stream',status='replace')
    do j=1,{nj}
        do i=1,{ni}
            {node_record('coarse', nh=nh)}
        enddo
    enddo
    if (refineRatio>1) then
        do j=1,ny
            do i=1,nx
'''
    if compact:
        for index, (r, condition) in enumerate([
            ('left', 'i<=nxLeft'), ('right', 'i>nx-nxRight'),
            ('bottom', 'j<=nyBottom'), ('top', 'j>ny-nyTop')]):
            body += ('                if' if index == 0 else '                else if')+f' ({condition}) then\n'
            body += f'                    ii=i-nint(xOffset{r.capitalize()}); jj=j-nint(yOffset{r.capitalize()})\n'
            body += '                    '+node_record(r, 'ii', 'jj')
        body += '                endif\n'
    else:
        body += '                if (.not.fine_active(i,j)) cycle\n                '+node_record('fine')
    body += '''            enddo
        enddo
    endif
    close(77)
    call calNuRe()
    call check()
    call output_SnapshotFile()
    call output_Tecplot()
    call output_ReloadFile()
    call exit_data_2d_openacc()
end program main'''
    return body


def compare_baseline(old, new, report):
    cases = [(1, False, False, False, None), (2, False, False, False, None),
             (2, True, True, False, None), (4, True, False, False, None),
             (8, False, True, False, None), (2, True, False, True, None),
             (2, True, False, False, (44,13,12,13))]
    for ratio, side, steady, ha, walls in cases:
        label=f'r{ratio}_side{int(side)}_steady{int(steady)}_ha{int(ha)}'+('_seam' if walls else '')
        dims=dict(n=128, ny=112, walls=(32,33,32,33), ra=10000) if ratio==8 else dict(
            n=96, ny=80, walls=(16,17,16,17), ra=10000)
        if walls: dims['walls']=walls
        kwargs=dict(ratio=ratio, side=side, steady=steady, ha=ha, legacy=True)
        results=[]
        for compact, source in [(False,old), (True,new)]:
            folder,exe=v.compile_source(label+('_compact' if compact else '_ring'),
                                        v.variant(source,**kwargs),driver(compact),**dims)
            v.run([exe],folder)
            results.append(folder)
        a,b=[np.fromfile(f/'canonical.bin',dtype='<f8') for f in results]
        if a.shape!=b.shape or not np.all(np.isfinite(b)):
            raise AssertionError((label,a.shape,b.shape,'nonfinite or incompatible state'))
        error=float(np.max(np.abs(a-b)))
        if not np.array_equal(a,b):
            np.savez(v.BUILD/(label+'_mismatch.npz'),before=a,after=b)
            raise AssertionError((label,'field/history mismatch',error,np.flatnonzero(a!=b)[:10]))
        for filename in ['NuRe_2DOpenaccMultiblock.dat']+(['Convergence_2DOpenaccMultiblock.dat'] if steady else []):
            a,b=[np.loadtxt(f/filename) for f in results]
            np.testing.assert_allclose(a,b,rtol=2e-12,atol=2e-12)
        converted=convert(checkpoint(results[0]),results[0]/'converted-v11.bin')
        # Header includes convergence sums, whose summation order has changed.
        assert converted.read_bytes()[272:]==checkpoint(results[1]).read_bytes()[272:]
        validate_snapshot(next(results[1].glob('*Snapshot-*.bin')))
        report['checks'].append(dict(case=label,fine_steps=40,fields_and_history_bitwise_equal=True,
                                     diagnostics_tolerance=2e-12,converted_v10_payload_bitwise_equal=True))
        print('PASS exact fields/history and diagnostics:',label,flush=True)


def checkpoint(folder):
    return folder/(folder/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().strip()


def validate_snapshot(path):
    with path.open('rb') as f:
        assert f.read(16).rstrip()==b'MB2DSNAPSHOT0003'
        count,nx,ny,itc=np.fromfile(f,dtype='<i4',count=4)
        np.fromfile(f,dtype='<f8',count=2)
        area=0.
        for _ in range(count):
            ni,nj=np.fromfile(f,dtype='<i4',count=2)
            np.fromfile(f,dtype='<f8',count=7+ni+nj)
            weights=np.fromfile(f,dtype='<f8',count=ni*nj)
            assert np.all(weights>=0)
            area+=weights.sum()
            fields=np.fromfile(f,dtype='<f8',count=4*ni*nj)
            assert fields.size==4*ni*nj and np.all(np.isfinite(fields))
        assert not f.read(1) and area==nx*ny


def actual_main(new,report):
    dims=dict(n=96,ny=80,walls=(16,17,16,17),ra=10000)
    def settings(steady,end,restart=False):
        s=v.variant(new,ratio=2,side=True,steady=steady,restart=restart)
        s=re.sub(r'(::\s*(?:outputSnapshotInterval|reloadFileInterval|outputPltFileInterval)\s*=)\s*[^\n!]+',
                 r'\1 0.02d0 ',s)
        if steady:
            s=re.sub(r'(::\s*itc_max\s*=)\s*\d+',r'\g<1> '+str(end),s)
            s=re.sub(r'(::\s*eps[UT]\s*=)\s*[^\n]+',r'\1 0.0d0',s)
            # Keep the short convergence-check cadence independent of source formatting.
            s=re.sub(r'\bmod\s*\(\s*itc\s*,\s*2000\s*\)', 'mod(itc, 20)', s, flags=re.I)
        else:
            s=re.sub(r'(::\s*unsteadyRunDuration\s*=)\s*[^\n!]+',r'\1 '+str(end)+'d0 ',s)
        return s
    for steady in [False,True]:
        label='steady' if steady else 'unsteady'
        half,full=(80,160) if steady else (0.05,0.1)
        continuous,exe=v.compile_source('main_'+label+'_full',settings(steady,full),**dims)
        v.run([exe],continuous)
        split,exe=v.compile_source('main_'+label+'_split',settings(steady,half),**dims)
        v.run([exe],split)
        _,exe=v.compile_source('main_'+label+'_resume',settings(steady,full,True),**dims)
        v.run([exe],split)
        assert checkpoint(continuous).read_bytes()[272:]==checkpoint(split).read_bytes()[272:]
        for filename in ['NuRe_2DOpenaccMultiblock.dat']+(['Convergence_2DOpenaccMultiblock.dat'] if steady else []):
            assert (continuous/filename).read_bytes()==(split/filename).read_bytes()
        for path in split.glob('*Snapshot-*.bin'): validate_snapshot(path)
        report['checks'].append(dict(actual_main_v11_restart=label,fields_history_and_diagnostics_exact=True))
        print('PASS actual main continuous versus resumed:',label,flush=True)
    # The diagnostic history is independent of all three file-output switches.
    reference=(continuous/'NuRe_2DOpenaccMultiblock.dat').read_bytes()
    for flag,pattern in [('outputSnapshotFile','*Snapshot-*.bin'),('outputPltFile','*Tecplot-*.dat'),
                         ('outputReloadFile','reloadFile2DOpenaccMultiblock-*')]:
        s=settings(True,160).replace(flag+' = 1',flag+' = 0')
        folder,exe=v.compile_source('disabled_'+flag,s,**dims)
        v.run([exe],folder)
        assert not list(folder.glob(pattern))
        assert (folder/'NuRe_2DOpenaccMultiblock.dat').read_bytes()==reference
    report['checks'].append(dict(independent_output_switches=True,snapshot_area_exact=True))
    print('PASS independent output switches and snapshot areas',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline',type=Path,nargs='?',help='Optional ring source for exact before/after comparison')
    args=parser.parse_args()
    new=v.SOURCE.read_text(encoding='utf-8-sig')
    report=dict(source_sha256=hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
                device='gfortran OpenACC host',
                build_directory=str(v.BUILD),checks=[])
    if args.baseline:
        old=args.baseline.read_text(encoding='utf-8-sig')
        report.update(baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
                      baseline_path=str(args.baseline.resolve()))
        compare_baseline(old,new,report)
    actual_main(new,report)
    report['status']='passed'
    (v.HERE/'compact_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
