#!/usr/bin/env python3
"""Rebuild D2Q5 FLOW-TEST statistics from complete 500--1000 t_ff histories."""

import argparse
import array
import hashlib
import math
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Tuple


START_TF = 500.0
MID_TF = 750.0
END_TF = 1000.0
EXPECTED_COMPLETE = 501
EXPECTED_FIRST = 250
EXPECTED_FINAL = 251
PRANDTL = 0.7
MACH = 0.1
CS2 = 1.0 / 3.0
THOT = 0.5
TCOLD = -0.5
DELTA_T = THOT - TCOLD


def read_table(path: Path, columns: int) -> List[List[float]]:
    rows = []  # type: List[List[float]]
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != columns:
                raise ValueError(f"{path}:{line_number}: expected {columns} columns, got {len(fields)}")
            rows.append([float(value) for value in fields])
    return rows


def assert_close(left: float, right: float, label: str, tolerance: float = 1.0e-10) -> None:
    if abs(left - right) > tolerance * max(1.0, abs(left), abs(right)):
        raise ValueError(f"{label}: {left:.17e} != {right:.17e}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, content: str) -> None:
    temporary = path.with_name(path.name + ".tmp-750-1000")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(content)
    os.replace(temporary, path)


def format_real(value: float) -> str:
    return f"{value: .16E}"


def parse_profiles(path: Path, ny: int, expected_times: List[float]) -> Tuple[List[List[float]], List[List[float]]]:
    values = array.array("d")
    with path.open("rb") as stream:
        values.fromfile(stream, path.stat().st_size // values.itemsize)
    if sys.byteorder != "little":
        values.byteswap()
    record_size = 1 + 2 * ny
    expected_values = len(expected_times) * record_size
    if len(values) != expected_values:
        raise ValueError(f"{path}: expected {expected_values} doubles, got {len(values)}")
    temperature = []  # type: List[List[float]]
    temperature_squared = []  # type: List[List[float]]
    for record, expected_time in enumerate(expected_times):
        offset = record * record_size
        profile_time = values[offset]
        assert_close(profile_time, expected_time, f"temperature profile time at record {record + 1}")
        temperature.append(list(values[offset + 1 : offset + 1 + ny]))
        temperature_squared.append(list(values[offset + 1 + ny : offset + 1 + 2 * ny]))
    return temperature, temperature_squared


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--rayleigh", type=float, required=True)
    parser.add_argument("--branch", choices=("ORIGINAL", "EFFECTIVE"), required=True)
    args = parser.parse_args()

    results = args.results.resolve()
    nx = args.nx
    ny = args.ny
    rayleigh = args.rayleigh
    length_unit = float(ny)

    paths = {
        "nu": results / "Nu_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",
        "re": results / "Re_VolAvg_2DOpenaccLBMCDE_D2Q5.dat",
        "diss": results / "DissipationHistory_2DOpenaccLBMCDE.dat",
        "profiles": results / "TemperatureProfileHistory_2DOpenaccLBMCDE.bin",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    nu_rows = read_table(paths["nu"], 2)
    re_rows = read_table(paths["re"], 2)
    diss_rows = read_table(paths["diss"], 10)
    if not (len(nu_rows) == len(re_rows) == len(diss_rows) == EXPECTED_COMPLETE):
        raise ValueError(
            f"history row counts must all be {EXPECTED_COMPLETE}: "
            f"Nu={len(nu_rows)}, Re={len(re_rows)}, dissipation={len(diss_rows)}"
        )

    times = []  # type: List[float]
    for index, (nu_row, re_row, diss_row) in enumerate(zip(nu_rows, re_rows, diss_rows)):
        expected_time = START_TF + float(index)
        assert_close(nu_row[0], expected_time, f"Nu time at row {index + 1}")
        assert_close(re_row[0], expected_time, f"Re time at row {index + 1}")
        assert_close(diss_row[0], expected_time, f"dissipation time at row {index + 1}")
        assert_close(nu_row[0], re_row[0], f"Nu/Re time at row {index + 1}")
        assert_close(nu_row[1], diss_row[1], f"Nu/dissipation value at row {index + 1}")
        assert_close(re_row[1], diss_row[2], f"Re/dissipation value at row {index + 1}")
        times.append(expected_time)

    first_indices = [i for i, value in enumerate(times) if START_TF <= value < MID_TF]
    final_indices = [i for i, value in enumerate(times) if MID_TF <= value <= END_TF]
    if len(first_indices) != EXPECTED_FIRST or len(final_indices) != EXPECTED_FINAL:
        raise ValueError(f"window counts are first={len(first_indices)}, final={len(final_indices)}")

    temperature, temperature_squared = parse_profiles(paths["profiles"], ny, times)

    nu_values = [row[1] for row in nu_rows]
    re_values = [row[1] for row in re_rows]

    def mean(values: List[float], indices: List[int]) -> float:
        return math.fsum(values[index] for index in indices) / float(len(indices))

    all_indices = list(range(EXPECTED_COMPLETE))
    nu_whole = mean(nu_values, all_indices)
    nu_first = mean(nu_values, first_indices)
    nu_final = mean(nu_values, final_indices)
    re_whole = math.sqrt(mean([value * value for value in re_values], all_indices))
    re_first = math.sqrt(mean([value * value for value in re_values], first_indices))
    re_final = math.sqrt(mean([value * value for value in re_values], final_indices))
    nu_error_percent = 100.0 * abs(nu_first - nu_final) / max(abs(nu_final), sys.float_info.min)
    re_error_percent = 100.0 * abs(re_first - re_final) / max(abs(re_final), sys.float_info.min)
    converged = nu_error_percent < 1.0 and re_error_percent < 1.0

    time_average_lines = [
        "# 2D OpenACC Nu mean and literature-defined RMS Re in the averaging window",
        "# final Nu/Re use the second half; relative-error columns are percentages (%), target < 1.0",
        "# start_tf mid_tf end_tf whole_count first_count second_count Nu_mean_whole Re_rms_whole "
        "Nu_mean_first Re_rms_first Nu_final_second Re_final_second "
        "Nu_first_vs_final_percent Re_first_vs_final_percent",
        " ".join(
            [
                format_real(START_TF), format_real(MID_TF), format_real(END_TF),
                str(EXPECTED_COMPLETE), str(EXPECTED_FIRST), str(EXPECTED_FINAL),
                format_real(nu_whole), format_real(re_whole), format_real(nu_first), format_real(re_first),
                format_real(nu_final), format_real(re_final),
                format_real(nu_error_percent) + "%", format_real(re_error_percent) + "%",
            ]
        ),
        "# 达到统计收敛了。" if converged else "# 未达到统计收敛，当前计算时间不足，还需要继续往后算。",
    ]

    instantaneous_lines = [
        'TITLE = "2D OpenACC instantaneous Nu/Re in statistics window"',
        'VARIABLES = "time_over_tff" "NuVolumeInstantaneous" "ReInstantaneousRms"',
        'ZONE T="NuReInstantaneous", F=POINT',
    ]
    running_lines = [
        'TITLE = "2D OpenACC running statistics in averaging window"',
        'VARIABLES = "time_over_tff" "NuVolumeRunningMean" "ReRunningRms"',
        'ZONE T="NuReRunningStatistics", F=POINT',
    ]
    nu_running_sum = 0.0
    re_squared_running_sum = 0.0
    for count, (time_tf, nu_value, re_value) in enumerate(zip(times, nu_values, re_values), 1):
        instantaneous_lines.append(f"{format_real(time_tf)} {format_real(nu_value)} {format_real(re_value)}")
        nu_running_sum += nu_value
        re_squared_running_sum += re_value * re_value
        running_lines.append(
            f"{format_real(time_tf)} {format_real(nu_running_sum / count)} "
            f"{format_real(math.sqrt(max(0.0, re_squared_running_sum / count)))}"
        )

    final_diss = [diss_rows[index] for index in final_indices]
    nu_volume = math.fsum(row[1] for row in final_diss) / EXPECTED_FINAL
    re_rms = math.sqrt(math.fsum(row[2] * row[2] for row in final_diss) / EXPECTED_FINAL)
    eps_u = math.fsum(row[3] for row in final_diss) / EXPECTED_FINAL
    eps_t = math.fsum(row[4] for row in final_diss) / EXPECTED_FINAL
    rho_rms = math.sqrt(math.fsum(row[5] * row[5] for row in final_diss) / EXPECTED_FINAL)
    div_rms = math.sqrt(math.fsum(row[6] * row[6] for row in final_diss) / EXPECTED_FINAL)
    maximum_cfl = max(row[7] * math.sqrt(CS2) for row in final_diss)
    assert_close(nu_volume, nu_final, "final Nu from Nu and dissipation histories")
    assert_close(re_rms, re_final, "final Re from Re and dissipation histories")

    mean_temperature = [
        math.fsum(temperature[index][j] for index in final_indices) / EXPECTED_FINAL for j in range(ny)
    ]
    mean_temperature_squared = [
        math.fsum(temperature_squared[index][j] for index in final_indices) / EXPECTED_FINAL for j in range(ny)
    ]
    theta_rms = [
        math.sqrt(max(0.0, mean_temperature_squared[j] - mean_temperature[j] * mean_temperature[j]))
        for j in range(ny)
    ]

    viscosity = MACH * length_unit * math.sqrt(PRANDTL / (3.0 * rayleigh))
    diffusivity = viscosity / PRANDTL
    g_beta = rayleigh * viscosity * diffusivity / (DELTA_T * length_unit**3)
    time_unit = math.sqrt(length_unit / (g_beta * DELTA_T))
    wall_spacing = 1.0 / length_unit
    nu_bottom = (8.0 * THOT - 9.0 * mean_temperature[0] + mean_temperature[1]) / (
        3.0 * wall_spacing * DELTA_T
    )
    nu_top = (-8.0 * TCOLD + 9.0 * mean_temperature[-1] - mean_temperature[-2]) / (
        3.0 * wall_spacing * DELTA_T
    )
    nu_wall = 0.5 * (nu_bottom + nu_top)
    nu_kinetic = 1.0 + eps_u * length_unit**4 * PRANDTL**2 / (viscosity**3 * rayleigh)
    nu_thermal = eps_t * length_unit**2 / (diffusivity * DELTA_T**2)
    wall_difference = 100.0 * abs(nu_wall - nu_volume) / abs(nu_volume)
    kinetic_difference = 100.0 * abs(nu_kinetic - nu_volume) / abs(nu_volume)
    thermal_difference = 100.0 * abs(nu_thermal - nu_volume) / abs(nu_volume)
    eps_u_exact = viscosity**3 / length_unit**4 * (nu_volume - 1.0) * rayleigh / PRANDTL**2
    eps_t_exact = diffusivity * DELTA_T**2 / length_unit**2 * nu_volume
    ratio_u = eps_u / eps_u_exact
    ratio_t = eps_t / eps_t_exact

    lower_peak = max(range(max(1, ny // 2)), key=lambda index: theta_rms[index])
    upper_peak = max(range(min(ny, ny // 2 + 1) - 1, ny), key=lambda index: theta_rms[index])
    lower_peak_fortran = lower_peak + 1
    upper_peak_fortran = upper_peak + 1
    delta_rms = 0.5 * (
        (float(lower_peak_fortran) - 0.5) + (float(ny - upper_peak_fortran) + 0.5)
    ) / length_unit
    n_bl_rms = max(1, math.ceil(delta_rms * length_unit))
    delta_global = 1.0 / (2.0 * nu_volume)
    n_bl_global = max(1, math.ceil(delta_global * length_unit))
    delta_difference = 100.0 * abs(delta_rms - delta_global) / delta_global
    eta_over_h = math.sqrt(PRANDTL) / (rayleigh * (nu_volume - 1.0)) ** 0.25
    eta_b_over_h = eta_over_h / math.sqrt(PRANDTL)
    grid_over_eta = (1.0 / length_unit) / eta_over_h
    grid_over_eta_b = (1.0 / length_unit) / eta_b_over_h
    dt_over_tau_eta = (1.0 / time_unit) / math.sqrt(PRANDTL / (nu_volume - 1.0))

    profile_lines = [
        "# z_over_H theta_rms",
        f"# delta_theta_global_estimate_over_H {format_real(delta_global)}",
        f"# delta_theta_rms_over_H {format_real(delta_rms)}",
        f"# N_BL_global_estimate {n_bl_global}",
        f"# N_BL_from_temperature_rms_peak {n_bl_rms}",
        f"# delta_theta_relative_difference_percent {format_real(delta_difference)}%",
    ]
    for j, value in enumerate(theta_rms, 1):
        profile_lines.append(f"{format_real((float(j) - 0.5) / length_unit)} {format_real(value)}")

    nint = lambda value: max(1, int(math.floor(value + 0.5)))
    statistics_lines = [
        "# Xu/Zhang final statistics over the second half of the configured history window",
        "# Nu/Re convergence still compares the first and second halves of the complete history",
        "# final statistics are trustworthy only when both Nu/Re relative errors are below 1%",
        f"# convergence_status {'CONVERGED' if converged else 'NOT_CONVERGED'}",
        f"# complete_history_samples {EXPECTED_COMPLETE}",
        f"# final_statistics_samples {EXPECTED_FINAL}",
        f"# history_start_mid_end_sample_itc {nint(START_TF * time_unit)} "
        f"{nint(MID_TF * time_unit)} {nint(END_TF * time_unit)}",
        f"# final_start_end_time_tf {format_real(MID_TF)} {format_real(END_TF)}",
        "# Xu et al. Table 4: Nu_volume is the reference for all relative differences",
        f"Nu_volume_from_heat_flux {format_real(nu_volume)}",
        f"Re_rms_space_time {format_real(re_rms)}",
        f"rho_fluctuation_rms_space_time {format_real(rho_rms)}",
        f"velocity_divergence_rms_space_time {format_real(div_rms)}",
        f"Nu_wall_from_mean_wall_gradient {format_real(nu_wall)}",
        f"Nu_wall_relative_difference_percent {format_real(wall_difference)}%",
        f"Nu_kinetic_from_mean_dissipation {format_real(nu_kinetic)}",
        f"Nu_kinetic_relative_difference_percent {format_real(kinetic_difference)}%",
        f"Nu_thermal_from_mean_dissipation {format_real(nu_thermal)}",
        f"Nu_thermal_relative_difference_percent {format_real(thermal_difference)}%",
        "# Zhang et al.: directly calculated dissipation over exact-relation dissipation",
        f"eps_u_calculated_space_time_mean {format_real(eps_u)}",
        f"eps_u_exact_from_Nu {format_real(eps_u_exact)}",
        f"R_u_calculated_over_exact {format_real(ratio_u)}",
        f"eps_T_calculated_space_time_mean {format_real(eps_t)}",
        f"eps_T_exact_from_Nu {format_real(eps_t_exact)}",
        f"R_T_calculated_over_exact {format_real(ratio_t)}",
        "# Nu/Re half-window convergence is rebuilt once from the complete Nu/Re history",
        "# See NuRe_TimeAverage_2DOpenaccLBMCDE_D2Q5.txt for first/second-half results",
        "# Zhang Table 1: N_BL uses the global estimate delta_theta/H=1/(2*Nu_volume)",
        f"delta_theta_global_estimate_over_H {format_real(delta_global)}",
        f"N_BL_global_estimate {n_bl_global}",
        "# Retained diagnostic: thermal boundary layer from the T_rms peak position",
        f"delta_theta_rms_over_H {format_real(delta_rms)}",
        f"N_BL_from_temperature_rms_peak {n_bl_rms}",
        f"delta_theta_relative_difference_percent {format_real(delta_difference)}%",
        "temperature_rms_profile_file TemperatureRmsProfile_2DOpenaccLBMCDE_D2Q5.dat",
        f"grid_over_eta_grid_over_etaB_dt_over_tauEta {format_real(grid_over_eta)} "
        f"{format_real(grid_over_eta_b)} {format_real(dt_over_tau_eta)}",
        f"maximum_lattice_CFL {format_real(maximum_cfl)}",
    ]

    output_names = [
        "NuReDissStatistics_2DOpenaccLBMCDE.dat",
        "NuRe_TimeAverage_2DOpenaccLBMCDE_D2Q5.txt",
        "NuRe_InstantaneousVolAvg_2DOpenaccLBMCDE_D2Q5.plt",
        "NuRe_RunningStatistics_2DOpenaccLBMCDE_D2Q5.plt",
        "TemperatureRmsProfile_2DOpenaccLBMCDE_D2Q5.dat",
    ]
    backup = results / "pre_750_1000_postprocess_20260817"
    backup.mkdir(exist_ok=True)
    for name in output_names:
        source = results / name
        target = backup / name
        if source.exists() and not target.exists():
            shutil.copy2(source, target)

    atomic_write_text(results / output_names[0], "\n".join(statistics_lines) + "\n")
    atomic_write_text(results / output_names[1], "\n".join(time_average_lines) + "\n")
    atomic_write_text(results / output_names[2], "\n".join(instantaneous_lines) + "\n")
    atomic_write_text(results / output_names[3], "\n".join(running_lines) + "\n")
    atomic_write_text(results / output_names[4], "\n".join(profile_lines) + "\n")

    manifest_lines = [
        "postprocess_version 750_1000_v1",
        f"generated_utc {datetime.now(timezone.utc).isoformat()}",
        f"branch {args.branch}",
        f"rayleigh {rayleigh:.1e}",
        f"nx {nx}",
        f"ny {ny}",
        f"complete_history_samples {EXPECTED_COMPLETE}",
        f"first_window_samples {EXPECTED_FIRST}",
        f"final_window_samples {EXPECTED_FINAL}",
        f"convergence_status {'CONVERGED' if converged else 'NOT_CONVERGED'}",
        f"Nu_first_vs_final_percent {nu_error_percent:.16e}",
        f"Re_first_vs_final_percent {re_error_percent:.16e}",
    ]
    for label, path in paths.items():
        manifest_lines.append(f"input_sha256_{label} {sha256(path)}")
    for name in output_names:
        manifest_lines.append(f"output_sha256_{name} {sha256(results / name)}")
    atomic_write_text(results / "postprocess_750_1000.manifest.txt", "\n".join(manifest_lines) + "\n")

    print(
        f"{args.branch} Ra={rayleigh:.1e}: Nu={nu_final:.12g}, Re={re_final:.12g}, "
        f"errors={nu_error_percent:.6g}%/{re_error_percent:.6g}%, "
        f"{'CONVERGED' if converged else 'NOT_CONVERGED'}"
    )


if __name__ == "__main__":
    main()
