"""Compare direct in-device exchange against an explicitly supplied pre-refactor solver.

Usage: python verify_direct_exchange.py PATH_TO_PRE_CHANGE_SOURCE
Builds and cases stay in a temporary directory; requires gfortran OpenACC host and numpy.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

import numpy as np

import verify_compact as c
import verify_multiblock as v
import verify_simplified as s


REGIONS = ('coarse', 'left', 'right', 'bottom', 'top')


def driver(packed, steps):
    body = f'''program main
    use commondata
    use openacc
    implicit none
    integer :: step,k
    call acc_init(acc_device_host)
    call initial()
    call enter_data_2d_openacc()
    do step=1,{steps}/refineRatio
        call advance_multiblock()
    enddo
    call update_host_all(.true.)
    open(unit=77,file='canonical.bin',form='unformatted',access='stream',status='replace')
'''
    for i, region in enumerate(REGIONS):
        if i == 1:
            body += '    if(refineRatio>1) then\n'
        body += '    write(77) '+','.join(name+'_'+region for name in
                                        ('f', 'g', 'rho', 'u', 'v', 'T', 'Fx', 'Fy'))+'\n'
    body += '    endif\n'
    if packed:
        body += '    write(77) coarseBuffer\n'
    else:
        body += '''    write(77) (rhoTimeInterp_coarse(:,:,k),uTimeInterp_coarse(:,:,k),&
        vTimeInterp_coarse(:,:,k),TTimeInterp_coarse(:,:,k),FxOverDxTimeInterp_coarse(:,:,k),&
        FyOverDxTimeInterp_coarse(:,:,k),flowNeqTimeInterp_coarse(:,:,:,k),&
        thermalNeqTimeInterp_coarse(:,:,:,k),k=0,timeInterpMaxIndexCoarse)
'''
    body += '''    close(77)
'''
    if steps:
        body += '''    if(mod(reloadItcOffset+itc,outputSnapshotIntervalItc)==0) then
        call calNuRe()
#ifdef steadyFlow
        call check()
#endif
    endif
    call output_SnapshotFile()
    call output_Tecplot()
    call output_ReloadFile()
'''
    body += '''    call exit_data_2d_openacc()
    call deallocate_grid_arrays()
end program main'''
    return body


def copy_checkpoint(folder, target):
    for path in folder.iterdir():
        if path.suffix in ('.dat', '.bin', '.meta', '.txt'):
            shutil.copy2(path, target/path.name)


def canonical(folder):
    return np.fromfile(folder/'canonical.bin', dtype='<f8')


def polynomial_driver():
    body = '''program main
    use commondata
    use openacc
    implicit none
    integer :: i,j,k,count
    real(8) :: x,y,tau,error,expectedRho,expectedT
    logical :: skin
    call acc_init(acc_device_host)
    call initial()
    call enter_data_2d_openacc()
    !$acc wait(1)
    coarseBuffer=0.0d0
    do k=0,2
        tau=dble(k-1)
        do j=1,size(f_coarse,2)
            y=(yOffsetCoarse+(j-0.5d0)*dxCoarse)/dble(ny)
            do i=1,size(f_coarse,1)
                x=(xOffsetCoarse+(i-0.5d0)*dxCoarse)/dble(nx)
                coarseBuffer(i,j,0,k)=1.0d0+(0.01d0*x**3*y**3)*(1.0d0+tau+tau**2)
                coarseBuffer(i,j,3,k)=0.02d0*x**3+0.03d0*y**3+0.04d0*tau**2
            enddo
        enddo
    enddo
    !$acc update device(coarseBuffer) async(1)
'''
    for region in REGIONS[1:]:
        body += f'''    rho_{region}=-123.0d0
    T_{region}=-456.0d0
    !$acc update device(rho_{region},T_{region}) async(1)
'''
    body += '''    tau=0.37d0
    call coarse_to_fine([0.5d0*tau*(tau-1.0d0),1.0d0-tau**2,0.5d0*tau*(tau+1.0d0)])
    call update_host_all(.false.)
    error=0.0d0
    count=0
'''
    for region in REGIONS[1:]:
        cap = region.capitalize()
        body += f'''    do j=1,size(f_{region},2)
        y=yOffset{cap}+j-0.5d0
        do i=1,size(f_{region},1)
            x=xOffset{cap}+i-0.5d0
            call fine_skin(x,y,skin)
            if(skin) then
                count=count+1
                expectedRho=1.0d0+0.01d0*(x/nx)**3*(y/ny)**3*(1.0d0+tau+tau**2)
                expectedT=0.02d0*(x/nx)**3+0.03d0*(y/ny)**3+0.04d0*tau**2
                error=max(error,abs(rho_{region}(i,j)-expectedRho),abs(T_{region}(i,j)-expectedT))
            else
                if(rho_{region}(i,j)/=-123.0d0 .OR. T_{region}(i,j)/=-456.0d0) &
                    error stop 'Interpolation overwrote a non-receiver node'
            endif
        enddo
    enddo
'''
    body += '''    if(count==0 .OR. error>1.0d-12) error stop 'Cubic space/quadratic time interpolation failed'
    open(unit=77,file='polynomial.dat',status='replace')
    write(77,*) count,error
    close(77)
    call exit_data_2d_openacc()
    call deallocate_grid_arrays()
end program main'''
    return body


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    args = parser.parse_args()
    old = args.baseline.read_text(encoding='utf-8-sig')
    new = v.SOURCE.read_text(encoding='utf-8-sig')
    report = dict(source_sha256=hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
                  baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
                  baseline_path=str(args.baseline.resolve()), device='gfortran OpenACC host',
                  build_directory=str(v.BUILD), checks=[])
    cases = [
        ('uniform_rb', 1, False, False, False, (16, 17, 16, 17)),
        ('r2_rb', 2, False, False, False, (16, 17, 16, 17)),
        ('r2_side_steady', 2, True, True, False, (16, 17, 16, 17)),
        ('r4_side', 4, True, False, False, (16, 17, 16, 17)),
        ('r8_rb_steady', 8, False, True, False, (32, 33, 32, 33)),
        ('r2_side_ha', 2, True, True, True, (16, 17, 16, 17)),
        ('r2_side_seam', 2, True, False, False, (44, 13, 12, 13)),
    ]
    try:
        for label, ratio, side, steady, ha, walls in cases:
            dims = dict(n=128, ny=112, ra=10000) if ratio == 8 else dict(n=96, ny=80, ra=10000)
            dims['walls'] = walls
            settings = dict(ratio=ratio, side=side, steady=steady, ha=ha)
            folders = []
            for tag, source, packed in [('before', old, False), ('after', new, True)]:
                folder, exe = v.compile_source(label+'_'+tag, s.short_settings(source, **settings),
                                               driver(packed, 40), **dims)
                v.run([exe], folder)
                folders.append(folder)
            before, after = folders
            field_error = s.compare_values(canonical(before), canonical(after), label+' fields/coarse history')
            diagnostics = ['NuRe_2DOpenaccMultiblock.dat']
            if steady:
                diagnostics.append('Convergence_2DOpenaccMultiblock.dat')
            errors = {name: s.compare_values(np.loadtxt(before/name), np.loadtxt(after/name), label+' '+name)
                      for name in diagnostics}
            snapshot_error = s.compare_snapshots(next(before.glob('*Snapshot-*.bin')),
                                                 next(after.glob('*Snapshot-*.bin')))
            restored, exe = v.compile_source(label+'_read_old', s.short_settings(new, restart=True, **settings),
                                             driver(True, 0), **dims)
            copy_checkpoint(before, restored)
            v.run([exe], restored)
            s.compare_values(canonical(before), canonical(restored), label+' restore old', exact=True)

            # Start from an old checkpoint inside a run; perform actual new exchanges after restoring.
            half_steps = (20//ratio)*ratio
            split, exe = v.compile_source(label+'_old_half', s.short_settings(old, **settings),
                                         driver(False, half_steps), **dims)
            v.run([exe], split)
            _, exe = v.compile_source(label+'_continue_old', s.short_settings(new, restart=True, **settings),
                                     driver(True, 40-half_steps), **dims)
            v.run([exe], split)
            continued_error = s.compare_values(canonical(after), canonical(split), label+' continue old')

            constant, exe = v.compile_source(label+'_constant', s.short_settings(new, **settings),
                                             s.constant_driver(), **dims)
            v.run([exe], constant)
            actual, expected = np.loadtxt(constant/'constant_diagnostics.dat')
            s.compare_values(expected, actual, label+' constant field')
            entry = dict(case=label, fields_coarse_history_max_abs_error=field_error,
                         diagnostics=errors, snapshot_max_abs_error=snapshot_error,
                         geometry_and_integration_areas_exact=True, read_old_checkpoint_exact=True,
                         continue_old_checkpoint_max_abs_error=continued_error)
            if ratio > 1:
                poly, exe = v.compile_source(label+'_polynomial', s.short_settings(new, **settings),
                                            polynomial_driver(), **dims)
                v.run([exe], poly)
                count, error = np.loadtxt(poly/'polynomial.dat')
                entry.update(polynomial_receivers=int(count), polynomial_max_abs_error=float(error),
                             non_receiver_nodes_unchanged=True)
            report['checks'].append(entry)
            print('PASS', label, 'fields/history, diagnostics, restart, area and interpolation', flush=True)
        report['analytical_statistics'] = s.verify_statistics(new)
        report['status'] = 'passed'
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = str(error)
        raise
    finally:
        path = v.HERE/'direct_exchange_verification.json'
        path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print('REPORT', path, flush=True)
        print('BUILD', v.BUILD, flush=True)


if __name__ == '__main__':
    main()
