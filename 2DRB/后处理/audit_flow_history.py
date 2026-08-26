#!/usr/bin/env python3
"""Audit FLOW-TEST history files without changing any result data.

The solver writes three formatted histories (Nu, Re and dissipation), a
stream-format temperature-profile history, and a restart ``*-latest.meta``
ledger.  This program checks that their sample clocks, sample counts and
finite values agree.  It is deliberately read-only: redirect stdout to a
``history_audit.txt`` file when a persistent audit record is wanted.

Examples
--------
Audit one completed result directory::

    python3 audit_flow_history.py --case /data2/XLLi/LBMCDE/FLOW-TEST/ORIGINAL/Ra1e7/chinu0.5/repair_from_t1400/results \
      --strict > /data2/XLLi/LBMCDE/FLOW-TEST/ORIGINAL/Ra1e7/chinu0.5/repair_from_t1400/results/history_audit.txt

Audit every result directory below FLOW-TEST::

    python3 audit_flow_history.py --root /data2/XLLi/LBMCDE/FLOW-TEST --strict \
      > /data2/XLLi/LBMCDE/FLOW-TEST/history_audit_all.txt

``--strict`` returns non-zero only for structural errors (missing/corrupt
files, NaN/Inf, mismatched time axes or internal gaps).  A one-sample tail
after the latest restart metadata is reported separately: it is a valid final
output record, but it is not a fully committed restart ledger and must not be
continued blindly.
"""

from __future__ import print_function

import argparse
import glob
import json
import math
import os
import statistics
import struct
import sys


FLOAT_TOLERANCE = 1.0e-8


def _as_float(token):
    """Read Fortran E/D exponent notation."""
    return float(token.replace("D", "E").replace("d", "e"))


def _first_match(directory, patterns):
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.join(directory, pattern)))
        if matches:
            return matches[0]
    return None


def _parse_meta(path):
    values = {}
    if path is None:
        return values
    with open(path, "r") as handle:
        for raw in handle:
            fields = raw.split(None, 1)
            if len(fields) == 2:
                values[fields[0]] = fields[1].strip()
    return values


def _read_text_history(path, minimum_columns):
    """Return (times, nonfinite_count, malformed_count)."""
    times = []
    nonfinite = 0
    malformed = 0
    with open(path, "r") as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) < minimum_columns:
                malformed += 1
                continue
            try:
                values = [_as_float(field) for field in fields]
            except ValueError:
                malformed += 1
                continue
            if not all(math.isfinite(value) for value in values):
                nonfinite += 1
            times.append(values[0])
    return times, nonfinite, malformed


def _sampling_description(times):
    """Describe monotonicity and internal gaps using a robust median clock."""
    result = {
        "first": None,
        "last": None,
        "median_dt": None,
        "gaps": [],
        "nonmonotone": 0,
    }
    if not times:
        return result
    result["first"] = times[0]
    result["last"] = times[-1]
    if len(times) < 2:
        return result
    deltas = [right - left for left, right in zip(times[:-1], times[1:])]
    positive = [delta for delta in deltas if delta > FLOAT_TOLERANCE]
    result["nonmonotone"] = sum(delta <= FLOAT_TOLERANCE for delta in deltas)
    if not positive:
        return result
    median_dt = statistics.median(positive)
    result["median_dt"] = median_dt
    if median_dt <= FLOAT_TOLERANCE:
        return result
    for left, right, delta in zip(times[:-1], times[1:], deltas):
        # 1.5 handles harmless rounding while exposing a missing whole sample.
        if delta > 1.5 * median_dt:
            result["gaps"].append((left, right))
    return result


def _read_profile_history(path, ny):
    """Read stream-format records: t, <T>_x(:), <T^2>_x(:)."""
    record_values = 1 + 2 * ny
    record_bytes = 8 * record_values
    size = os.path.getsize(path)
    if size % record_bytes:
        return [], 0, 1, size, record_bytes

    times = []
    nonfinite = 0
    malformed = 0
    unpacker = struct.Struct("<" + "d" * record_values)
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(record_bytes)
            if not chunk:
                break
            if len(chunk) != record_bytes:
                malformed += 1
                break
            values = unpacker.unpack(chunk)
            if not all(math.isfinite(value) for value in values):
                nonfinite += 1
            times.append(values[0])
    return times, nonfinite, malformed, size, record_bytes


def _axis_delta(reference, candidate):
    if len(reference) != len(candidate):
        return None
    if not reference:
        return 0.0
    return max(abs(left - right) for left, right in zip(reference, candidate))


def _format_gaps(gaps):
    if not gaps:
        return "[]"
    return "[" + ", ".join("{:.12g}->{:.12g}".format(left, right)
                           for left, right in gaps) + "]"


def _format_float(value):
    return "NA" if value is None else "{:.12g}".format(value)


class Reporter(object):
    def __init__(self):
        self.errors = []
        self.warnings = []

    def error(self, text):
        self.errors.append(text)

    def warning(self, text):
        self.warnings.append(text)


def _report_history(label, times, nonfinite, malformed, reporter):
    detail = _sampling_description(times)
    print("{} rows {} first {} last {} median_dt {} gaps {} nonfinite {} malformed {}".format(
        label, len(times), _format_float(detail["first"]), _format_float(detail["last"]),
        _format_float(detail["median_dt"]), _format_gaps(detail["gaps"]), nonfinite, malformed))
    if not times:
        reporter.error(label + " has no data rows")
    if nonfinite:
        reporter.error(label + " has non-finite values")
    if malformed:
        reporter.error(label + " has malformed records")
    if detail["nonmonotone"]:
        reporter.error(label + " is not strictly increasing in time")
    if detail["gaps"]:
        reporter.error(label + " has internal time gaps")
    return detail


def _audit_case(directory):
    reporter = Reporter()
    print("CASE {}".format(os.path.abspath(directory)))

    merge_manifest = None
    merge_manifest_path = os.path.join(directory, "merge_manifest.json")
    if os.path.isfile(merge_manifest_path):
        with open(merge_manifest_path, "r") as handle:
            merge_manifest = json.load(handle)
        print("MERGE_MANIFEST {}".format(merge_manifest_path))

    meta_path = _first_match(directory, ["*-latest.meta"])
    meta = _parse_meta(meta_path)
    print("META file {}".format(meta_path if meta_path else "MISSING"))
    if not meta_path:
        reporter.error("latest.meta is missing")
    else:
        for key in ("reload_meta_version", "nx", "ny", "time_tf",
                    "cumulativeStatisticSampleCount", "reloadFileName"):
            print("META {} {}".format(key, meta.get(key, "MISSING")))

    nu_path = _first_match(directory, ["Nu_VolAvg_*.dat"])
    re_path = _first_match(directory, ["Re_VolAvg_*.dat"])
    diss_path = _first_match(directory, ["DissipationHistory_*.dat"])
    profile_path = _first_match(directory, ["TemperatureProfileHistory_*.bin"])
    sources = []
    for label, path, columns in (("NU", nu_path, 2), ("RE", re_path, 2),
                                 ("DISS", diss_path, 10)):
        if path is None:
            print("{} file MISSING".format(label))
            reporter.error(label + " history is missing")
            sources.append((label, []))
        else:
            times, nonfinite, malformed = _read_text_history(path, columns)
            _report_history(label, times, nonfinite, malformed, reporter)
            sources.append((label, times))

    try:
        ny = int(meta.get("ny", ""))
    except ValueError:
        ny = 0
    profile_times = []
    if profile_path is None:
        print("PROFILE file MISSING")
        reporter.error("temperature profile history is missing")
    elif ny <= 0:
        print("PROFILE file {} cannot parse because META ny is invalid".format(profile_path))
        reporter.error("temperature profile cannot be checked without valid META ny")
    else:
        profile_times, nonfinite, malformed, byte_count, record_bytes = _read_profile_history(profile_path, ny)
        print("PROFILE file {} bytes {} record_bytes {}".format(profile_path, byte_count, record_bytes))
        _report_history("PROFILE", profile_times, nonfinite, malformed, reporter)
        sources.append(("PROFILE", profile_times))

    populated = [(label, times) for label, times in sources if times]
    if populated:
        reference_name, reference_times = populated[0]
        for label, times in populated[1:]:
            delta = _axis_delta(reference_times, times)
            if delta is None:
                print("AXIS {}_vs_{} length_mismatch {} {}".format(
                    reference_name, label, len(reference_times), len(times)))
                reporter.error("{} and {} have different sample counts".format(reference_name, label))
            else:
                print("AXIS {}_vs_{} max_abs_time_difference {:.12g}".format(
                    reference_name, label, delta))
                if delta > FLOAT_TOLERANCE:
                    reporter.error("{} and {} have different time axes".format(reference_name, label))

    if meta_path and populated:
        try:
            committed = int(meta.get("cumulativeStatisticSampleCount", ""))
        except ValueError:
            committed = None
        if committed is None or committed < 0:
            reporter.error("META cumulativeStatisticSampleCount is invalid")
        else:
            row_count = len(populated[0][1])
            ledger_rows = row_count
            prepended_rows = 0
            if merge_manifest:
                try:
                    latest_rows = int(merge_manifest["histories"]["Nu"]["sources"][-1]["rows"])
                    ledger_rows = latest_rows
                    prepended_rows = row_count - latest_rows
                except (KeyError, IndexError, TypeError, ValueError):
                    reporter.error("merge manifest does not describe the latest Nu source rows")
                print("META_COMMIT merged_history_rows {} latest_part_rows {} prepended_rows {}".format(
                    row_count, ledger_rows, prepended_rows))
            tail = ledger_rows - committed
            print("META_COMMIT committed_samples {} ledger_history_rows {} tail_rows {}".format(
                committed, ledger_rows, tail))
            if committed > ledger_rows:
                reporter.error("META claims more committed samples than the histories contain")
            elif tail == 0:
                print("META_COMMIT_STATUS exact")
            elif tail == 1:
                print("META_COMMIT_STATUS one_uncommitted_final_tail")
                reporter.warning("one history sample follows latest.meta; do not continue from this meta without rollback")
            else:
                print("META_COMMIT_STATUS mismatch")
                reporter.error("latest.meta and histories differ by more than one sample")

            if "time_tf" in meta and reference_times:
                try:
                    meta_time = _as_float(meta["time_tf"])
                    last_time = reference_times[-1]
                    clock = _sampling_description(reference_times)["median_dt"] or 0.0
                    print("META_TIME meta_tf {:.12g} history_last_tf {:.12g} abs_delta {:.12g}".format(
                        meta_time, last_time, abs(meta_time - last_time)))
                    if abs(meta_time - last_time) > max(0.55 * clock, FLOAT_TOLERANCE):
                        reporter.warning("latest.meta time and last history time are more than half a sample apart")
                except ValueError:
                    reporter.error("META time_tf is invalid")

    if meta_path and meta.get("reloadFileName"):
        restart_name = meta["reloadFileName"]
        restart_matches = glob.glob(os.path.join(directory, "*" + restart_name + ".bin"))
        print("RESTART latest_bin_matches {}".format(len(restart_matches)))
        if not restart_matches:
            reporter.error("latest.meta points to a missing restart binary")

    if reporter.errors:
        verdict = "FAIL"
    elif reporter.warnings:
        verdict = "PASS_WITH_WARNING"
    else:
        verdict = "PASS"
    print("VERDICT {}".format(verdict))
    for text in reporter.errors:
        print("ERROR {}".format(text))
    for text in reporter.warnings:
        print("WARNING {}".format(text))
    print()
    return reporter.errors


def _discover_cases(root):
    cases = set()
    for pattern in ("**/*-latest.meta", "**/Nu_VolAvg_*.dat", "**/DissipationHistory_*.dat"):
        for path in glob.glob(os.path.join(root, pattern), recursive=True):
            cases.add(os.path.dirname(path))
    return sorted(cases)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--case", help="one result directory (or a case directory containing results/)")
    group.add_argument("--root", help="recursively audit every result directory under this root")
    parser.add_argument("--strict", action="store_true",
                        help="return non-zero when a structural consistency error is found")
    args = parser.parse_args(argv)

    if args.case:
        candidate = os.path.abspath(args.case)
        nested_results = os.path.join(candidate, "results")
        cases = [nested_results if os.path.isdir(nested_results) else candidate]
    else:
        cases = _discover_cases(os.path.abspath(args.root))
    if not cases:
        print("No result directories found.", file=sys.stderr)
        return 2

    all_errors = []
    for directory in cases:
        all_errors.extend(_audit_case(directory))
    print("SUMMARY cases {} structural_errors {}".format(len(cases), len(all_errors)))
    return 1 if args.strict and all_errors else 0


if __name__ == "__main__":
    sys.exit(main())
