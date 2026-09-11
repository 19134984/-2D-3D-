"""Check local-gradient and wall-flux diagnostics in a 2-D LBM-CDE Tecplot file.

The script intentionally uses only the Python standard library so the same
post-processing can run on the legacy P100 nodes without installing NumPy.
"""

import argparse
from pathlib import Path



def parse_record(record):
    """Parse a whitespace-delimited or legacy fixed-width Tecplot record."""

    fields = record.split()
    if len(fields) == 8:
        return [float(value) for value in fields]
    if len(record) < 192:
        raise ValueError(f"cannot parse Tecplot record: {record!r}")
    return [float(record[offset : offset + 24]) for offset in range(0, 192, 24)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tecplot", type=Path)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--hot", type=float, default=0.5)
    parser.add_argument("--cold", type=float, default=-0.5)
    args = parser.parse_args()

    temperature = []
    grad_x_local = []
    grad_y_local = []
    with args.tecplot.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream):
            if line_number < 3:
                continue
            record = line.rstrip("\r\n")
            if not record:
                continue
            row = parse_record(record)
            temperature.append(row[4])
            grad_x_local.append(row[6])
            grad_y_local.append(row[7])

    expected_rows = args.nx * args.ny
    if len(temperature) != expected_rows:
        raise ValueError(f"expected {expected_rows} rows, got {len(temperature)}")

    local_squared = 0.0
    fd_squared = 0.0
    local_fd_dot = 0.0
    for j in range(2, args.ny - 2):
        row_offset = j * args.nx
        for i in range(2, args.nx - 2):
            index = row_offset + i
            grad_x_fd = 0.5 * (temperature[index + 1] - temperature[index - 1])
            grad_y_fd = 0.5 * (
                temperature[index + args.nx] - temperature[index - args.nx]
            )
            local_squared += (
                grad_x_local[index] * grad_x_local[index]
                + grad_y_local[index] * grad_y_local[index]
            )
            fd_squared += grad_x_fd * grad_x_fd + grad_y_fd * grad_y_fd
            local_fd_dot += (
                grad_x_local[index] * grad_x_fd
                + grad_y_local[index] * grad_y_fd
            )

    local_norm = local_squared**0.5
    fd_norm = fd_squared**0.5
    correlation = local_fd_dot / (local_norm * fd_norm)

    delta_temperature = args.hot - args.cold
    bottom = temperature[: args.nx]
    top = temperature[(args.ny - 1) * args.nx :]
    nu_hot_halfway = (
        2.0
        * args.ny
        * sum(args.hot - value for value in bottom)
        / (args.nx * delta_temperature)
    )
    nu_cold_halfway = (
        2.0
        * args.ny
        * sum(value - args.cold for value in top)
        / (args.nx * delta_temperature)
    )

    print(
        {
            "file": str(args.tecplot),
            "interior_local_to_fd_gradient_norm": local_norm / fd_norm,
            "interior_gradient_correlation": correlation,
            "Nu_hot_halfway_temperature": nu_hot_halfway,
            "Nu_cold_halfway_temperature": nu_cold_halfway,
            "temperature_min": min(temperature),
            "temperature_max": max(temperature),
        }
    )


if __name__ == "__main__":
    main()
