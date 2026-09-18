"""Compare ordinary-array storage with an explicitly supplied pre-change solver.

Usage: python verify_plain_arrays.py PATH_TO_PRE_CHANGE_SOURCE
All generated sources, checkpoints and outputs stay in a temporary directory.
"""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import verify_multiblock as v


def checkpoint(folder):
    return folder / (folder / 'reloadFile2DOpenaccMultiblock-latest.meta').read_text().strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('--report-name', default='plain_arrays_verification.json',
                        help='Report filename within the solver directory')
    args = parser.parse_args()
    old = args.baseline.read_text(encoding='utf-8-sig')
    new = v.SOURCE.read_text(encoding='utf-8-sig')
    report = {
        'source_sha256': hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
        'baseline_sha256': hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
        'baseline_path': str(args.baseline.resolve()),
        'device': 'gfortran OpenACC host',
        'build_directory': str(v.BUILD),
        'checks': [],
    }
    # Exercise field history, both output formats, diagnostics and steady errors.
    diagnostics = '''
    call calNuRe()
    call check()
    call output_SnapshotFile()
    call output_Tecplot()
'''
    output_driver = v.RESTART_DRIVER.replace('    call output_ReloadFile()', diagnostics+'    call output_ReloadFile()')
    dims = dict(n=96, ny=80, walls=(16,17,16,17), ra=10000)
    for ratio, side, steady in [(1,False,False),(2,False,False),(4,True,False),(8,True,True)]:
        label = f'r{ratio}_side{int(side)}_steady{int(steady)}'
        kw = dict(ratio=ratio, side=side, steady=steady, legacy=True)
        case_dims = dict(n=128, ny=112, walls=(32,33,32,33), ra=10000) if ratio == 8 else dims
        before, be = v.compile_source(label+'_before', v.variant(old, **kw), output_driver, **case_dims)
        after, ae = v.compile_source(label+'_after', v.variant(new, **kw), output_driver, **case_dims)
        v.run([be,'40'], before)
        v.run([ae,'40'], after)
        state = np.fromfile(after/'allstate.bin', dtype='<f8')
        assert np.all(np.isfinite(state))
        assert (before/'allstate.bin').read_bytes() == (after/'allstate.bin').read_bytes()
        assert checkpoint(before).read_bytes() == checkpoint(after).read_bytes()
        files = ['NuRe_2DOpenaccMultiblock.dat',
                 'buoyancyCavity2DOpenaccMultiblockSnapshot-0000000001.bin',
                 'buoyancyCavity2DOpenaccMultiblockTecplot-0000000001.dat']
        if steady:
            files.append('Convergence_2DOpenaccMultiblock.dat')
        for filename in files:
            assert (before/filename).read_bytes() == (after/filename).read_bytes(), filename
        report['checks'].append({'case':label, 'fine_steps':40,
                                 'fields_packets_checkpoint_outputs_bitwise_equal':True})
        print('Fields, packet history, diagnostics and outputs identical:', label, flush=True)

    # Real main loops: identical scheduled output and bidirectional v10 resume.
    import re
    def settings(source, steady, restart):
        src = v.variant(source, ratio=2, side=True, steady=steady, legacy=True, restart=restart)
        src = re.sub(r'(::\s*(?:outputSnapshotInterval|reloadFileInterval|outputPltFileInterval)\s*=)\s*[^\n!]+',
                     r'\1 0.05d0 ', src)
        if steady:
            src = re.sub(r'(::\s*itc_max\s*=)\s*\d+',r'\1 '+('320' if restart else '160'),src)
            src = re.sub(r'(::\s*eps[UT]\s*=)\s*[^\n]+',r'\1 0.0d0',src)
            src = src.replace('mod(itc, 2000)', 'mod(itc, 20)')
        else:
            src = re.sub(r'(::\s*unsteadyRunDuration\s*=)\s*[^\n!]+',r'\1 '+('0.2d0 ' if restart else '0.1d0 '),src)
        return src
    for steady in [False, True]:
        label = 'steady' if steady else 'unsteady'
        outputs = []
        initial_checkpoint = None
        for source_name, source, other in [('before',old,new),('after',new,old)]:
            folder, exe = v.compile_source(f'main_{label}_{source_name}',settings(source,steady,False), **dims)
            v.run([exe],folder)
            if initial_checkpoint is None:
                initial_checkpoint = checkpoint(folder).read_bytes()
            else:
                assert initial_checkpoint == checkpoint(folder).read_bytes()
            _, rexe = v.compile_source(f'main_{label}_{source_name}_cross_resume',settings(other,steady,True), **dims)
            v.run([rexe],folder)
            outputs.append(folder)
        assert checkpoint(outputs[0]).read_bytes() == checkpoint(outputs[1]).read_bytes()
        files = ['NuRe_2DOpenaccMultiblock.dat']
        if steady:
            files.append('Convergence_2DOpenaccMultiblock.dat')
        files += [f.name for pattern in ['buoyancy*Snapshot-*.bin','buoyancy*Tecplot-*.dat']
                  for f in outputs[0].glob(pattern)]
        for filename in files:
            assert (outputs[0]/filename).read_bytes() == (outputs[1]/filename).read_bytes(), filename
        report['checks'].append({'actual_main_bidirectional_v10_restart':label,
                                 'checkpoints_counters_histories_outputs_bitwise_equal':True})
        print('Scheduled main-loop bidirectional v10 restart identical:', label, flush=True)
    report['status'] = 'passed'
    if Path(args.report_name).name != args.report_name:
        raise ValueError('report-name must be a filename, not a path')
    (v.HERE/args.report_name).write_text(json.dumps(report, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
