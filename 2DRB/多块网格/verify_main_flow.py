"""Small OpenACC-host checks for final diagnostics, outputs and restart counters.

Run: python verify_main_flow.py
All edited sources and output files stay in the temporary build directory.
"""
from pathlib import Path
import re
import struct

import numpy as np

import verify_compact as c
import verify_multiblock as v


def settings(source, end, restart=False, converge=False):
    source = v.variant(source, ratio=2, side=True, steady=True, restart=restart)
    source = re.sub(r'(::\s*itc_max\s*=)\s*\d+', rf'\g<1> {end}', source)
    tolerance = '200.0d0' if converge else '0.0d0'
    source = re.sub(r'(::\s*eps[UT]\s*=)\s*[^\n]+', rf'\g<1> {tolerance}', source)
    for name, steps in [('outputSnapshotInterval', 5.5), ('reloadFileInterval', 9.5),
                        ('outputPltFileInterval', 7.5)]:
        source = re.sub(r'(::\s*'+name+r'\s*=)\s*[^\n!]+',
                        rf'\g<1> {steps}d0 / timeUnit ', source)
    source = re.sub(r'(^[ \t]*checkIntervalItc[ \t]*=[ \t]*)[^\n]+',
                    r'\g<1> 10', source, flags=re.M)
    # Prove the real main program releases host arrays after device cleanup.
    source = source.replace('    call deallocate_grid_arrays()', '''    call deallocate_grid_arrays()
    if( allocated(f_coarse) .OR. allocated(f_left) .OR. allocated(f_right) .OR. &
        allocated(f_bottom) .OR. allocated(f_top) .OR. allocated(rhoReceive) .OR. &
        allocated(coarseToLeftReceiveIndexX) .OR. allocated(up_coarse) ) error stop 'Host arrays still allocated' ''')
    return source


def run_case(source, label, end, restart=False, folder=None, converge=False):
    build, exe = v.compile_source(label, settings(source, end, restart, converge),
                                 n=96, ny=80, walls=(16, 17, 16, 17), ra=10000)
    folder = folder or build
    v.run([exe], folder)
    return folder


def snapshot_steps(folder):
    steps = []
    for path in folder.glob('*Snapshot-*.bin'):
        c.validate_snapshot(path)
        steps.append(struct.unpack_from('<i', path.read_bytes(), 12)[0])
    assert len(steps) == len(set(steps)), 'Repeated snapshot at the same step'
    return sorted(steps)


def main():
    source = v.SOURCE.read_text(encoding='utf-8-sig')
    split = run_case(source, 'off_sample_half', 22)
    assert snapshot_steps(split) == [6, 12, 18, 22]
    history = np.loadtxt(split/'NuRe_2DOpenaccMultiblock.dat')
    final = np.loadtxt(split/'FinalNuRe_2DOpenaccMultiblock.dat')
    assert len(history) == 3 and final[0] > history[-1, 0]
    # Header still stores next sample number = completed sample count + 1.
    itc, sample_index = struct.unpack_from('<2i', c.checkpoint(split).read_bytes(), 212)
    assert (itc, sample_index) == (22, 4)

    continuous = run_case(source, 'off_sample_full', 40)
    run_case(source, 'off_sample_resume', 40, restart=True, folder=split)
    assert c.checkpoint(continuous).read_bytes()[256:] == c.checkpoint(split).read_bytes()[256:]
    for name in ('NuRe_2DOpenaccMultiblock.dat', 'SteadyMonitor_2DOpenaccMultiblock.dat',
                 'Convergence_2DOpenaccMultiblock.dat', 'FinalNuRe_2DOpenaccMultiblock.dat'):
        assert (continuous/name).read_bytes() == (split/name).read_bytes(), name
    assert snapshot_steps(split) == [6, 12, 18, 22, 24, 30, 36, 40]
    log = (split/'SimulationSettings2DOpenaccMultiblock.txt').read_text()
    work = [float(x) for x in re.findall(r'Lattice updates in this run\s*=\s*(\S+)', log)]
    per_cycle = 37*29 + 2*(18*80 + 19*80 + 59*18 + 59*19)
    assert work == [11*per_cycle, 9*per_cycle], work
    assert log.count('End:') == 2
    assert log.count('steady convergence not achieved') == 2
    print('PASS off-sample final diagnostics, exact restart, host cleanup and run-only work count', flush=True)

    aligned = run_case(source, 'aligned_final', 24)
    assert snapshot_steps(aligned) == [6, 12, 18, 24]
    plots = list(aligned.rglob('*Tecplot*.dat'))
    assert len(plots) == 3, plots
    print('PASS final snapshot/Tecplot do not duplicate periodic outputs', flush=True)

    stopped = run_case(source, 'zero_step_converged', 24, converge=True)
    log = (stopped/'SimulationSettings2DOpenaccMultiblock.txt').read_text()
    assert 'Steady calculation converged.' in log
    assert re.search(r'Lattice updates in this run\s*=\s*0\.0', log)
    assert snapshot_steps(stopped) == [0]
    assert '# converged = T' in (stopped/'FinalNuRe_2DOpenaccMultiblock.dat').read_text()
    print('PASS zero-step termination and convergence status', flush=True)
    print('BUILD', v.BUILD, flush=True)


if __name__ == '__main__':
    main()
