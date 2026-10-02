"""Check cumulative-tff filenames and old/new checkpoint selection on OpenACC host.

Usage: python verify_output_names.py PATH_TO_PRE_CHANGE_SOURCE
All compiled sources and solver outputs stay in a temporary directory.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import struct

import verify_compact as c
import verify_main_flow as m
import verify_multiblock as v


META = 'reloadFile2DOpenaccMultiblock-latest.meta'
DIMS = dict(n=96, ny=80, walls=(16, 17, 16, 17), ra=10000)


def manual_name(source, name):
    source, count = re.subn(r'(::\s*reloadFileName\s*=)"[^"]*"',
                            lambda match: match[1] + '"' + name + '"', source)
    assert count == 1
    return source


def settings(source, end, steady=True, restart=False):
    source = m.settings(source, end, restart)
    if not steady:
        source = re.sub(r'\.OR\.\s*(?:&\s*)?allocated\(up_coarse\)', '', source)
        source = re.sub(r'^#define steadyFlow\s*$', '!#define steadyFlow', source, flags=re.M)
        source = re.sub(r'^!#define unsteadyFlow\s*$', '#define unsteadyFlow', source, flags=re.M)
        source = re.sub(r'(::\s*unsteadyRunDuration\s*=)\s*[^\n!]+',
                        rf'\g<1> {end}.0d0/timeUnit ', source)
    # Expose the actual scaling for checking filenames, including Tecplot.
    return source.replace('    call initial()', '''    call initial()
    open(unit=79,file='filename_time_unit.dat',status='replace')
    write(79,*) timeUnit
    close(79)''', 1)


def build(source, label, end, steady=True, restart=False, name=None):
    if name is not None:
        source = manual_name(source, name)
    return v.compile_source(label, settings(source, end, steady, restart), **DIMS)


def check_names(folder):
    time_unit = float((folder/'filename_time_unit.dat').read_text())
    snapshots = list(folder.glob('*Snapshot-*tff.bin'))
    checkpoints = list(folder.glob('reloadFile*-*tff.bin'))
    plots = list(folder.glob('*Tecplot-*tff.dat'))
    assert snapshots and checkpoints and plots
    for paths in (snapshots, checkpoints, plots):
        times = []
        for path in paths:
            match = re.search(r'-(\d+(?:\.\d+)?)tff\.(?:bin|dat)$', path.name)
            assert match, path.name
            value = float(match[1])
            assert '.' not in match[1] or not match[1].endswith('0'), path.name
            times.append(value)
            if 'Snapshot-' in path.name:
                data = path.read_bytes()
                step = struct.unpack_from('<i', data, 12)[0]
                recorded_time = struct.unpack_from('<d', data, 16)[0]
                assert abs(value-recorded_time) <= 5.1e-13
                c.validate_snapshot(path)
            elif path.name.startswith('reloadFile'):
                data = path.read_bytes()
                step = struct.unpack_from('<i', data, 212)[0]
                assert struct.unpack_from('<3i', data, 228) == (0, 0, 0)
            else:
                step = round(value*time_unit)
            assert abs(value-step/time_unit) <= 5.1e-13, path.name
            assert step % 2 == 0
        assert len(times) == len(set(times))
    assert c.checkpoint(folder).name.endswith('tff.bin')
    return dict(snapshots=len(snapshots), checkpoints=len(checkpoints), plots=len(plots))


def compare_final(full, resumed):
    # The three retired naming counters occupy bytes 228:240; all other bytes agree.
    a, b = c.checkpoint(full).read_bytes(), c.checkpoint(resumed).read_bytes()
    assert a[:228] == b[:228] and a[240:] == b[240:]
    assert (full/'NuRe_2DOpenaccMultiblock.dat').read_bytes() == (
        resumed/'NuRe_2DOpenaccMultiblock.dat').read_bytes()
    assert c.checkpoint(full).name == c.checkpoint(resumed).name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    args = parser.parse_args()
    source = v.SOURCE.read_text(encoding='utf-8-sig')
    old = args.baseline.read_text(encoding='utf-8-sig')
    report = dict(source_sha256=hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
                  baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
                  compiler=v.COMPILER, device='OpenACC host', build_directory=str(v.BUILD), checks=[])
    for steady in (True, False):
        mode = 'steady' if steady else 'unsteady'
        full, exe = build(source, mode+'_full', 40, steady, name='does-not-exist.bin')
        v.run([exe], full)  # Fresh runs ignore manual restart settings.
        split, exe = build(source, mode+'_split', 22, steady)
        v.run([exe], split)
        prior = c.checkpoint(split)
        saved_name, saved_data = prior.name, prior.read_bytes()
        manual = v.BUILD/(mode+'_manual')
        shutil.copytree(split, manual)
        (manual/META).unlink()
        resume_end = 18 if steady else 40
        _, resume = build(source, mode+'_resume', resume_end, steady, True, saved_name)
        v.run([resume], split)
        v.run([resume], manual)
        for folder in (split, manual):
            compare_final(full, folder)
            assert (folder/saved_name).read_bytes() == saved_data
            check_names(folder)
        report['checks'].append(dict(mode=mode, continuous_and_resumed_exact=True,
                                     manual_filename_exact=True, previous_checkpoint_preserved=True,
                                     files=check_names(full)))
        print('PASS cumulative tff, continuous/resumed, manual filename:', mode, flush=True)

        # The old numbered checkpoint must remain readable with either selection path.
        legacy, exe = build(old, mode+'_legacy', 22, steady)
        v.run([exe], legacy)
        legacy_name = c.checkpoint(legacy).name
        assert not legacy_name.endswith('tff.bin')
        legacy_manual = v.BUILD/(mode+'_legacy_manual')
        shutil.copytree(legacy, legacy_manual)
        (legacy_manual/META).unlink()
        # A deliberately missing manual filename proves that latest.meta takes priority.
        _, automatic = build(source, mode+'_meta_priority', resume_end, steady, True, 'missing.bin')
        v.run([automatic], legacy)
        _, manual_exe = build(source, mode+'_legacy_resume', resume_end, steady, True, legacy_name)
        v.run([manual_exe], legacy_manual)
        for folder in (legacy, legacy_manual):
            compare_final(full, folder)
        report['checks'].append(dict(mode=mode, old_numbered_checkpoint_exact=True,
                                     old_manual_filename_exact=True, meta_priority=True))
        print('PASS old numbered checkpoints and latest.meta priority:', mode, flush=True)

    zero, exe = build(source, 'zero_time', 0)
    v.run([exe], zero)
    check_names(zero)
    assert c.checkpoint(zero).name.endswith('-0tff.bin')
    assert all('-0tff.' in p.name for p in zero.glob('*Snapshot-*.bin'))
    report['checks'].append(dict(zero_time=True))
    report_path = v.HERE/'output_names_verification.json'
    report_path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('PASS zero-time output; report:', report_path, flush=True)
    print('BUILD', v.BUILD, flush=True)


if __name__ == '__main__':
    main()
