#!/usr/bin/env python3
"""Collect the six chinu=0 FLOW-TEST postprocessing outputs into one TSV."""

from __future__ import print_function

import argparse
import os
from pathlib import Path


CASES = (
    ("ORIGINAL", "Ra1e6", 129, "1e6"),
    ("ORIGINAL", "Ra1e7", 257, "1e7"),
    ("ORIGINAL", "Ra1e8", 513, "1e8"),
    ("EFFECTIVE", "Ra1e6", 129, "1e6"),
    ("EFFECTIVE", "Ra1e7", 257, "1e7"),
    ("EFFECTIVE", "Ra1e8", 513, "1e8"),
)


def read_key_values(path):
    values = {}
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("# convergence_status "):
                values["convergence_status"] = stripped.split()[-1]
                continue
            if stripped.startswith("#"):
                continue
            fields = stripped.split()
            values[fields[0]] = [value.rstrip("%") for value in fields[1:]]
    return values


def read_nure(path):
    with path.open("r", encoding="utf-8") as stream:
        rows = [line.strip() for line in stream if line.strip() and not line.lstrip().startswith("#")]
    if len(rows) != 1:
        raise ValueError("{}: expected one data row, got {}".format(path, len(rows)))
    fields = [value.rstrip("%") for value in rows[0].split()]
    if len(fields) != 14:
        raise ValueError("{}: expected 14 columns, got {}".format(path, len(fields)))
    return fields


def one(values, key):
    return values[key][0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()

    header = [
        "branch", "Ra", "grid", "convergence_status",
        "Nu_whole_500_1000", "Re_whole_500_1000",
        "Nu_first_500_750", "Re_first_500_750",
        "Nu_final_750_1000", "Re_final_750_1000",
        "Nu_first_vs_final_percent", "Re_first_vs_final_percent",
        "Nu_wall", "Nu_wall_diff_percent", "Nu_kinetic", "Nu_kinetic_diff_percent",
        "Nu_thermal", "Nu_thermal_diff_percent", "eps_u", "eps_u_exact", "R_u",
        "eps_T", "eps_T_exact", "R_T", "N_BL_global", "N_BL_rms",
        "grid_over_eta", "grid_over_etaB", "dt_over_tauEta", "maximum_lattice_CFL",
    ]
    rows = []
    for branch, ra_dir, grid, ra_label in CASES:
        results = root / branch / ra_dir / "chinu0.0" / "results"
        manifest = results / "postprocess_750_1000.manifest.txt"
        if not manifest.is_file():
            raise FileNotFoundError(manifest)
        stats = read_key_values(results / "NuReDissStatistics_2DOpenaccLBMCDE.dat")
        nure = read_nure(results / "NuRe_TimeAverage_2DOpenaccLBMCDE_D2Q5.txt")
        scales = stats["grid_over_eta_grid_over_etaB_dt_over_tauEta"]
        rows.append([
            branch, ra_label, "{}x{}".format(grid, grid), stats["convergence_status"],
            nure[6], nure[7], nure[8], nure[9], nure[10], nure[11], nure[12], nure[13],
            one(stats, "Nu_wall_from_mean_wall_gradient"),
            one(stats, "Nu_wall_relative_difference_percent"),
            one(stats, "Nu_kinetic_from_mean_dissipation"),
            one(stats, "Nu_kinetic_relative_difference_percent"),
            one(stats, "Nu_thermal_from_mean_dissipation"),
            one(stats, "Nu_thermal_relative_difference_percent"),
            one(stats, "eps_u_calculated_space_time_mean"), one(stats, "eps_u_exact_from_Nu"),
            one(stats, "R_u_calculated_over_exact"),
            one(stats, "eps_T_calculated_space_time_mean"), one(stats, "eps_T_exact_from_Nu"),
            one(stats, "R_T_calculated_over_exact"), one(stats, "N_BL_global_estimate"),
            one(stats, "N_BL_from_temperature_rms_peak"), scales[0], scales[1], scales[2],
            one(stats, "maximum_lattice_CFL"),
        ])

    output = root / "analysis" / "chinu0.0_reprocessed_750_1000.tsv"
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("\t".join(header) + "\n")
        for row in rows:
            stream.write("\t".join(row) + "\n")
    os.replace(str(temporary), str(output))
    print(output)
    for row in rows:
        print("{} {} {} Nu={} Re={} errors={}%/{}%".format(
            row[0], row[1], row[3], row[8], row[9], row[10], row[11]
        ))


if __name__ == "__main__":
    main()
