"""Check parameter/integration simplification against a supplied pre-change solver.

Usage: python verify_simplified.py PATH_TO_PRE_CHANGE_SOURCE [--report PATH]
Runs temporary, bounds-checked OpenACC host builds; this is not GPU validation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil

import numpy as np

import verify_compact as c
import verify_multiblock as v


STEPS = 40
TOL = 1e-12


def short_settings(source, **settings):
    source = v.variant(source, **settings)
    # Make the diagnostic emitted by c.driver coincide with a scheduled sample,
    # so the resulting checkpoint satisfies the real restart clock contract.
    return re.sub(
        r'(::\s*(?:outputSnapshotInterval|reloadFileInterval|outputPltFileInterval)\s*=)\s*[^\n!]+',
        rf'\1 {STEPS - 0.25}d0 / timeUnit ', source,
    )


def compare_values(before, after, label, exact=False):
    assert before.shape == after.shape, (label, before.shape, after.shape)
    assert np.isfinite(before).all() and np.isfinite(after).all(), label
    if exact:
        np.testing.assert_array_equal(before, after, err_msg=label)
    else:
        np.testing.assert_allclose(after, before, rtol=TOL, atol=TOL, err_msg=label)
    return float(np.max(np.abs(after - before))) if before.size else 0.0


def snapshot_records(path):
    c.validate_snapshot(path)
    records = []
    with path.open('rb') as stream:
        head = np.fromfile(stream, dtype='<i4', count=4)
        time = np.fromfile(stream, dtype='<f8', count=2)
        for _ in range(head[0]):
            shape = np.fromfile(stream, dtype='<i4', count=2)
            ni, nj = shape
            coords = np.fromfile(stream, dtype='<f8', count=7)
            widths = np.fromfile(stream, dtype='<f8', count=ni + nj)
            areas = np.fromfile(stream, dtype='<f8', count=ni * nj)
            fields = np.fromfile(stream, dtype='<f8', count=4 * ni * nj)
            records.append((shape, coords, widths, areas, fields))
        assert not stream.read(1)
    return head, time, records


def compare_snapshots(before, after):
    a, b = snapshot_records(before), snapshot_records(after)
    compare_values(a[0], b[0], 'snapshot header', exact=True)
    compare_values(a[1], b[1], 'snapshot time', exact=True)
    assert len(a[2]) == len(b[2])
    error = 0.0
    for old, new in zip(a[2], b[2]):
        for name, x, y in zip(('shape', 'coordinates', 'widths', 'cell areas'), old[:4], new[:4]):
            compare_values(x, y, 'snapshot ' + name, exact=True)
        error = max(error, compare_values(old[4], new[4], 'snapshot fields'))
    return error


def constant_driver():
    body = '''program main
    use commondata
    use openacc
    implicit none
    real(8) :: actual(11), expected(11), scale, speed
    call acc_init(acc_device_host)
    call initial()
    call enter_data_2d_openacc()
'''
    for region in ('coarse', 'left', 'right', 'bottom', 'top'):
        if region == 'left':
            body += '    if (refineRatio > 1) then\n'
        body += f'''    rho_{region}=1.25d0
    u_{region}=0.003d0
    v_{region}=0.004d0
    T_{region}=0.125d0
    !$acc update device(rho_{region},u_{region},v_{region},T_{region})
'''
    body += '''    endif
    call calculate_NuRe(actual(1),actual(2),actual(3),actual(4),actual(5),actual(6), &
        actual(7),actual(8),actual(9),actual(10),actual(11))
    scale=lengthUnit/(Thot-Tcold)
#ifdef SideHeatedCell
    speed=0.003d0
#else
    speed=0.004d0
#endif
    expected(1)=1d0+speed*0.125d0*scale/diffusivity
    expected(2)=0.005d0*lengthUnit/viscosity
    expected(3)=8d0*(Thot-0.125d0)*scale/3d0
    expected(4)=8d0*(0.125d0-Tcold)*scale/3d0
    expected(5)=speed*0.125d0*scale/diffusivity
    expected(6)=1.25d0*dble(nx)*dble(ny)
    expected(7:9)=0.125d0
    expected(10:11)=1.25d0
    open(unit=78,file='constant_diagnostics.dat',status='replace')
    write(78,'(11(ES25.16E3,1X))') actual
    write(78,'(11(ES25.16E3,1X))') expected
    close(78)
    call output_SnapshotFile()
    call exit_data_2d_openacc()
    call deallocate_grid_arrays()
end program main'''
    return body


def verify_statistics(source):
    # Exercise interpolation at a half-window boundary between two samples.
    # Re^2 is linear from 1 to 9, so its full/half means are 5, 3, and 7.
    for name, value in [('unsteadyRunDuration', '2.0d0'),
                        ('unsteadyAverageStartTf', '0.0d0'),
                        ('unsteadyAverageEndTf', '2.0d0')]:
        source, count = re.subn(r'(::\s*' + name + r'\s*=)\s*[^\n!]+',
                               rf'\1 {value} ', source)
        assert count == 1, name
    driver = '''program main
    implicit none
    real(8) :: actual(5), coverage
    open(unit=78,file='analytical_statistics.dat',status='replace')
    call average_window(0d0,2d0,actual,coverage)
    write(78,'(6(ES25.16E3,1X))') coverage,actual
    call average_window(0d0,1d0,actual,coverage)
    write(78,'(6(ES25.16E3,1X))') coverage,actual
    call average_window(1d0,2d0,actual,coverage)
    write(78,'(6(ES25.16E3,1X))') coverage,actual
    close(78)
    call output_unsteady_NuRe_postprocess()
end program main'''
    folder, exe = v.compile_source('analytical_statistics', source, driver)
    history = np.array([[0, 2, 1, 4, 6, 8, 0, 0, 0, 0, 0, 0],
                        [2, 6, 3, 8, 10, 12, 0, 0, 0, 0, 0, 0]], dtype=float)
    np.savetxt(folder / 'NuRe_2DOpenaccMultiblock.dat', history,
               header='t_ff NuVolAvg ReVolRMS Nu_hot Nu_cold Nu_middle and six unused fields')
    v.run([exe], folder)
    expected = np.array([[2, 4, np.sqrt(5), 6, 8, 10],
                         [1, 3, np.sqrt(3), 5, 7, 9],
                         [1, 5, np.sqrt(7), 7, 9, 11]])
    actual = np.loadtxt(folder / 'analytical_statistics.dat')
    error = compare_values(expected, actual, 'analytical full/half-window statistics')
    output = (folder / 'NuReStatistics_2DOpenaccMultiblock.dat').read_text()
    expected_output = {
        'Covered durations, full/first/last:': expected[:, 0],
        'whole:': expected[0, 1:],
        'first:': expected[1, 1:],
        'last :': expected[2, 1:],
        'relative half-window difference:': np.array([
            0.5, (np.sqrt(7) - np.sqrt(3)) / np.sqrt(5), 1 / 3, 0.25, 0.2]),
    }
    for label, values in expected_output.items():
        line = next(line.strip() for line in output.splitlines() if line.strip().startswith(label))
        observed = np.fromstring(line[len(label):], sep=' ')
        error = max(error, compare_values(values, observed, 'analytical postprocess ' + label))
    return dict(max_abs_error=error, full_first_last_actual=actual.tolist(),
                full_first_last_expected=expected.tolist(),
                relative_half_window_expected=expected_output['relative half-window difference:'].tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path, help='Saved solver immediately before this simplification')
    parser.add_argument('--report', type=Path, default=v.HERE / 'simplified_verification.json',
                        help='JSON report path (default: simplified_verification.json beside this script)')
    args = parser.parse_args()
    old = args.baseline.read_text(encoding='utf-8-sig')
    new = v.SOURCE.read_text(encoding='utf-8-sig')
    report = dict(source_sha256=hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
                  baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
                  baseline_path=str(args.baseline.resolve()), device='gfortran OpenACC host',
                  build_directory=str(v.BUILD), tolerance=dict(rtol=TOL, atol=TOL), checks=[])
    cases = [
        ('uniform_rb', 1, False, False, False, None),
        ('r2_rb', 2, False, False, False, None),
        ('r2_side_steady', 2, True, True, False, None),
        ('r4_side', 4, True, False, False, None),
        ('r8_rb_steady', 8, False, True, False, (32, 33, 32, 33)),
        ('r2_side_ha', 2, True, True, True, None),
        ('r2_side_seam', 2, True, False, False, (44, 13, 12, 13)),
    ]
    try:
        for label, ratio, side, steady, ha, walls in cases:
            dims = dict(n=128, ny=112, ra=10000) if ratio == 8 else dict(n=96, ny=80, ra=10000)
            dims['walls'] = walls or (16, 17, 16, 17)
            settings = dict(ratio=ratio, side=side, steady=steady, ha=ha)
            folders = []
            for name, source in [('before', old), ('after', new)]:
                folder, exe = v.compile_source(label + '_' + name, short_settings(source, **settings),
                                               c.driver(True, steps=STEPS), **dims)
                v.run([exe], folder)
                folders.append(folder)
            before, after = folders
            field_error = compare_values(np.fromfile(before / 'canonical.bin', dtype='<f8'),
                                         np.fromfile(after / 'canonical.bin', dtype='<f8'), label + ' fields/history')
            files = ['NuRe_2DOpenaccMultiblock.dat']
            if steady:
                files.append('Convergence_2DOpenaccMultiblock.dat')
            diagnostics = {}
            for filename in files:
                diagnostics[filename] = compare_values(np.loadtxt(before / filename), np.loadtxt(after / filename),
                                                        label + ' ' + filename)
            snapshot_error = compare_snapshots(next(before.glob('*Snapshot-*.bin')),
                                               next(after.glob('*Snapshot-*.bin')))
            # Load the old checkpoint with the new source without advancing.  All
            # distribution, macro and interpolation-history values must survive.
            resumed, exe = v.compile_source(label + '_cross_resume',
                                            short_settings(new, restart=True, **settings),
                                            c.driver(True, steps=0), **dims)
            for path in before.iterdir():
                if path.suffix in ('.dat', '.bin', '.meta', '.txt'):
                    shutil.copy2(path, resumed / path.name)
            v.run([exe], resumed)
            compare_values(np.fromfile(before / 'canonical.bin', dtype='<f8'),
                           np.fromfile(resumed / 'canonical.bin', dtype='<f8'),
                           label + ' cross-version restart', exact=True)

            constant, exe = v.compile_source(label + '_constant', short_settings(new, **settings),
                                             constant_driver(), **dims)
            v.run([exe], constant)
            actual, expected = np.loadtxt(constant / 'constant_diagnostics.dat')
            constant_error = compare_values(expected, actual, label + ' constant-field diagnostics')
            c.validate_snapshot(next(constant.glob('*Snapshot-*.bin')))
            report['checks'].append(dict(case=label, fine_steps=STEPS,
                                         fields_history_max_abs_error=field_error,
                                         diagnostic_max_abs_errors=diagnostics,
                                         snapshot_fields_max_abs_error=snapshot_error,
                                         snapshot_geometry_and_cell_areas_exact=True,
                                         cross_version_restart_bitwise_equal=True,
                                         constant_field_max_abs_error=constant_error,
                                         constant_field_actual=actual.tolist(),
                                         constant_field_expected=expected.tolist()))
            print('PASS', label, 'fields/history, diagnostics, areas, old restart, constant fields', flush=True)
        report['analytical_statistics'] = verify_statistics(new)
        print('PASS analytical full/half-window statistics, RMS Re, relative differences', flush=True)
        report['status'] = 'passed'
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = str(error)
        raise
    finally:
        result = args.report
        result.write_text(json.dumps(report, indent=2), encoding='utf-8')
        print('Report:', result, flush=True)
        print('BUILD', v.BUILD, flush=True)


if __name__ == '__main__':
    main()
