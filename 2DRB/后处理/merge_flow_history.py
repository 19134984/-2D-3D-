#!/usr/bin/env python3
"""Non-destructively merge ordered FLOW-TEST continuation histories.

Later parts take precedence at overlapping sample times.  The program merges
Nu, Re, dissipation and binary temperature-profile histories, then copies the
latest result metadata and hard-links (or copies) the restart file referenced
by ``*-latest.meta``.  Existing output directories are never overwritten.
"""

from __future__ import print_function

import argparse
import glob
import hashlib
import json
import math
import os
import shutil
import struct
import sys


TIME_DIGITS = 8


def _time_key(value):
    return round(value, TIME_DIGITS)


def _first_match(directory, pattern):
    matches = sorted(glob.glob(os.path.join(directory, pattern)))
    if not matches:
        raise RuntimeError("missing {} in {}".format(pattern, directory))
    return matches[0]


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _read_text_history(path):
    headers = []
    records = {}
    with open(path, "r") as handle:
        for raw in handle:
            text = raw.strip()
            if not text or text.startswith("#"):
                headers.append(raw)
                continue
            fields = text.split()
            try:
                values = [float(field.replace("D", "E").replace("d", "e"))
                          for field in fields]
            except ValueError:
                raise RuntimeError("malformed history row in {}: {}".format(path, text))
            if not values or not all(math.isfinite(value) for value in values):
                raise RuntimeError("non-finite history row in {}: {}".format(path, text))
            records[_time_key(values[0])] = raw if raw.endswith("\n") else raw + "\n"
    return headers, records


def _merge_text(parts, pattern, output):
    merged = {}
    first_headers = []
    sources = []
    for index, directory in enumerate(parts):
        path = _first_match(directory, pattern)
        headers, records = _read_text_history(path)
        if index == 0:
            first_headers = headers
        before = len(merged)
        overlap = len(set(merged).intersection(records))
        merged.update(records)
        sources.append({
            "path": path,
            "sha256": _sha256(path),
            "rows": len(records),
            "overlap_with_previous": overlap,
            "new_rows": len(merged) - before,
        })
    times = sorted(merged)
    if len(times) > 1:
        deltas = [right - left for left, right in zip(times[:-1], times[1:])]
        positive = sorted(delta for delta in deltas if delta > 0.0)
        median = positive[len(positive) // 2]
        gaps = [(left, right) for left, right in zip(times[:-1], times[1:])
                if right - left > 1.5 * median]
    else:
        median = None
        gaps = []
    if gaps:
        raise RuntimeError("merged {} still has gaps: {}".format(pattern, gaps))
    with open(output, "w") as handle:
        handle.writelines(first_headers)
        for time_value in times:
            handle.write(merged[time_value])
    return {
        "output": output,
        "sha256": _sha256(output),
        "rows": len(times),
        "first": times[0] if times else None,
        "last": times[-1] if times else None,
        "median_dt": median,
        "gaps": gaps,
        "sources": sources,
    }


def _read_profile_history(path, ny):
    values_per_record = 1 + 2 * ny
    record_bytes = 8 * values_per_record
    size = os.path.getsize(path)
    if size % record_bytes:
        raise RuntimeError("profile byte count is not a whole record: {}".format(path))
    unpack_time = struct.Struct("<d")
    records = {}
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(record_bytes)
            if not chunk:
                break
            time_value = unpack_time.unpack(chunk[:8])[0]
            if not math.isfinite(time_value):
                raise RuntimeError("non-finite profile time in {}".format(path))
            records[_time_key(time_value)] = chunk
    return record_bytes, records


def _merge_profiles(parts, ny, output):
    merged = {}
    sources = []
    record_bytes = None
    for directory in parts:
        path = _first_match(directory, "TemperatureProfileHistory_*.bin")
        current_bytes, records = _read_profile_history(path, ny)
        if record_bytes is None:
            record_bytes = current_bytes
        elif record_bytes != current_bytes:
            raise RuntimeError("profile record sizes differ")
        before = len(merged)
        overlap = len(set(merged).intersection(records))
        merged.update(records)
        sources.append({
            "path": path,
            "sha256": _sha256(path),
            "rows": len(records),
            "overlap_with_previous": overlap,
            "new_rows": len(merged) - before,
        })
    times = sorted(merged)
    if len(times) > 1:
        deltas = [right - left for left, right in zip(times[:-1], times[1:])]
        positive = sorted(delta for delta in deltas if delta > 0.0)
        median = positive[len(positive) // 2]
        gaps = [(left, right) for left, right in zip(times[:-1], times[1:])
                if right - left > 1.5 * median]
    else:
        median = None
        gaps = []
    if gaps:
        raise RuntimeError("merged profile history still has gaps: {}".format(gaps))
    with open(output, "wb") as handle:
        for time_value in times:
            handle.write(merged[time_value])
    return {
        "output": output,
        "sha256": _sha256(output),
        "rows": len(times),
        "first": times[0] if times else None,
        "last": times[-1] if times else None,
        "median_dt": median,
        "gaps": gaps,
        "record_bytes": record_bytes,
        "sources": sources,
    }


def _copy_latest_artifacts(latest_part, output):
    copied = []
    auxiliary_patterns = [
        "*-latest.meta",
        "NuRe_TimeAverage_*.txt",
        "NuReDissStatistics_*.dat",
        "SimulationSettings*.txt",
        "runtime_source.sha256",
        "run.status",
        "compile.status",
        "timing.txt",
    ]
    for pattern in auxiliary_patterns:
        for source in sorted(glob.glob(os.path.join(latest_part, pattern))):
            target = os.path.join(output, os.path.basename(source))
            shutil.copy2(source, target)
            copied.append(target)

    meta_paths = sorted(glob.glob(os.path.join(latest_part, "*-latest.meta")))
    if meta_paths:
        restart_name = None
        with open(meta_paths[0], "r") as handle:
            for raw in handle:
                fields = raw.split(None, 1)
                if len(fields) == 2 and fields[0] == "reloadFileName":
                    restart_name = fields[1].strip()
                    break
        if restart_name:
            matches = sorted(glob.glob(os.path.join(
                latest_part, "reloadFile*" + restart_name + ".bin")))
            if not matches:
                raise RuntimeError("latest.meta points to a missing restart file")
            source = matches[0]
            target = os.path.join(output, os.path.basename(source))
            try:
                os.link(source, target)
            except OSError:
                shutil.copy2(source, target)
            copied.append(target)
    return copied


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="new merged result directory")
    parser.add_argument("--ny", required=True, type=int, help="number of y nodes")
    parser.add_argument("--parts", required=True, nargs="+",
                        help="ordered result directories; later parts win overlaps")
    args = parser.parse_args(argv)

    output = os.path.abspath(args.output)
    parts = [os.path.abspath(path) for path in args.parts]
    if os.path.exists(output):
        print("refusing to overwrite existing output: {}".format(output), file=sys.stderr)
        return 2
    for directory in parts:
        if not os.path.isdir(directory):
            print("missing input directory: {}".format(directory), file=sys.stderr)
            return 2

    os.makedirs(output)
    try:
        report = {
            "precedence": "later parts replace earlier records at equal times",
            "parts": parts,
            "ny": args.ny,
            "histories": {},
        }
        for label, pattern in (
                ("Nu", "Nu_VolAvg_*.dat"),
                ("Re", "Re_VolAvg_*.dat"),
                ("Dissipation", "DissipationHistory_*.dat")):
            name = os.path.basename(_first_match(parts[-1], pattern))
            report["histories"][label] = _merge_text(
                parts, pattern, os.path.join(output, name))
        profile_name = os.path.basename(
            _first_match(parts[-1], "TemperatureProfileHistory_*.bin"))
        report["histories"]["TemperatureProfile"] = _merge_profiles(
            parts, args.ny, os.path.join(output, profile_name))
        report["latest_artifacts"] = _copy_latest_artifacts(parts[-1], output)
        manifest = os.path.join(output, "merge_manifest.json")
        with open(manifest, "w") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print("MERGED {} parts into {}".format(len(parts), output))
        for label, detail in sorted(report["histories"].items()):
            print("{} rows {} first {} last {} gaps {}".format(
                label, detail["rows"], detail["first"], detail["last"], detail["gaps"]))
    except Exception:
        shutil.rmtree(output)
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
