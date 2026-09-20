"""Check asymmetric overlap geometry, interpolation, runtime and exact restart.

Usage: python verify_asymmetric.py [pre-change-compact.F90]
All builds/runs use temporary directories and OpenACC host bounds checking.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
import verify_multiblock as v
import verify_compact as c


def interface_probe():
    # Exercise every fine receiver, including corners and compact-array seams.
    code = '''
    open(unit=76,file='geometry.dat',status='replace')
    write(76,*) nxCoarse,nyCoarse,nxLeft,nyLeft,nxRight,nyRight,nxBottom,nyBottom,nxTop,nyTop
    close(76)
'''
    for region in ('Left', 'Right', 'Bottom', 'Top'):
        p = 'coarseTo'+region
        code += f'''
    do ii=1,{p}Count
        x=xOffset{region}+{p}Ti(ii)-0.5d0
        y=yOffset{region}+{p}Tj(ii)-0.5d0
        expected=(x/nx)**3*(y/ny)**3
        value=0.0d0
        if ({p}Same(ii)) then
            value=((xOffsetCoarse+({p}Si(ii)-0.5d0)*dxCoarse)/nx)**3 &
                 *((yOffsetCoarse+({p}Sj(ii)-0.5d0)*dxCoarse)/ny)**3
        else
            do b=1,4
                do a=1,4
                    value=value+{p}Wx(a,ii)*{p}Wy(b,ii) &
                      *((xOffsetCoarse+({p}Si(ii)+a-1.5d0)*dxCoarse)/nx)**3 &
                      *((yOffsetCoarse+({p}Sj(ii)+b-1.5d0)*dxCoarse)/ny)**3
                enddo
            enddo
        endif
        if (abs(value-expected)>1d-12) error stop 'Cubic interface reproduction failed'
    enddo
'''
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path, nargs='?')
    args = parser.parse_args()
    source = v.SOURCE.read_text(encoding='utf-8-sig')
    report = dict(source_sha256=hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
                  device='gfortran OpenACC host', build_directory=str(v.BUILD), checks=[])
    baseline = args.baseline.read_text(encoding='utf-8-sig') if args.baseline else None
    if baseline:
        report['baseline_sha256'] = hashlib.sha256(args.baseline.read_bytes()).hexdigest()
    driver = c.driver(True, steps=160)
    driver = driver.replace('integer :: i,j,k,ii,jj,step',
                            'integer :: i,j,k,ii,jj,step,a,b\n    real(8) :: x,y,value,expected')
    driver = driver.replace('call initial()', 'call initial()\n'+interface_probe())
    for ratio, side, steady in [(1,False,False),(2,False,False),(2,True,True),
                                (4,True,False),(8,False,True)]:
        dims = dict(n=128, ny=112, walls=(32,33,32,33), ra=10000)
        label = f'r{ratio}_side{int(side)}_steady{int(steady)}'
        settings = dict(ratio=ratio, side=side, steady=steady)
        # Uniform mode does not allocate coarse/fine interfaces.
        probe = driver if ratio>1 else c.driver(True,steps=160)
        folder, exe = v.compile_source(label,v.variant(source,**settings),probe,**dims)
        v.run([exe],folder)
        assert np.isfinite(np.fromfile(folder/'canonical.bin',dtype='<f8')).all()
        c.validate_snapshot(next(folder.glob('*Snapshot-*.bin')))
        if ratio>1:
            geom=np.loadtxt(folder/'geometry.dat',dtype=int)
            np.testing.assert_array_equal(geom, [64//ratio+5,48//ratio+5,
                                                 34,112,35,112,59,34,59,35])
        record=dict(case=label,fine_steps=160,finite=True,snapshot_area_exact=True,
                    cubic_interpolation_checked=ratio>1)
        if baseline:
            oldfolder, oldexe=v.compile_source(label+'_before',v.variant(baseline,**settings),
                                               c.driver(True,steps=160),**dims)
            v.run([oldexe],oldfolder)
            old=np.loadtxt(oldfolder/'NuRe_2DOpenaccMultiblock.dat')
            new=np.loadtxt(folder/'NuRe_2DOpenaccMultiblock.dat')
            assert np.isfinite(new).all()
            record['NuRe_before']=old.tolist()
            record['NuRe_after']=new.tolist()
            if ratio==1:
                assert (oldfolder/'canonical.bin').read_bytes()==(folder/'canonical.bin').read_bytes()
            if ratio==2 and not side:
                symmetric=re.sub(r'(fineOverlapCells\s*=)\s*2',r'\1 4',source)
                restored, ex=v.compile_source('restored_symmetric',v.variant(symmetric,**settings),
                                               c.driver(True,steps=160),**dims)
                v.run([ex],restored)
                assert (oldfolder/'canonical.bin').read_bytes()==(restored/'canonical.bin').read_bytes()
                record['restored_symmetric_bitwise_equal']=True
                # An old symmetric checkpoint must not be loaded into narrower fine arrays.
                resume=v.variant(source,**settings,restart=True)
                _, ex=v.compile_source('reject_old_geometry',resume,c.driver(True,steps=0),**dims)
                try:
                    v.run([ex],oldfolder)
                except RuntimeError as error:
                    assert 'Restart block layout mismatch' in str(error)
                else:
                    raise AssertionError('Old symmetric checkpoint was accepted')
                record['old_geometry_rejected']=True
        report['checks'].append(record)
        print('PASS',label,flush=True)
    c.actual_main(source,report)
    report['status']='passed'
    (v.HERE/'asymmetric_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('Report:',v.HERE/'asymmetric_verification.json',flush=True)


if __name__=='__main__':
    main()
