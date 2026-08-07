"""Symbolic checks for the D2Q5 LBM-CDE derivation.

The script uses the exact moment ordering adopted by
``均匀网格/2DRBOpenmp.F90``:

    n = [T, jx, jy, e, nu]^T

It verifies:

1. the local-gradient source collapses to an effective flux relaxation;
2. the fourth-order pure-diffusion coefficients reduce to Dubois and
   Lallemand (2009), Eqs. (40)-(42);
3. the constant-advection hydrodynamic eigenvalue through fourth order.

Run with the system Python that provides SymPy:

    D:\\App\\python\\python.exe tools\\derive_d2q5_cde.py

Use ``--scan`` to repeat the coarse von Neumann scan documented in the
derivation note.  The scan is illustrative, not a substitute for the target
case's exact velocity and parameter range.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import sympy as sp


I = sp.I


@dataclass
class Expansion:
    """Coefficients of log(lambda_h(epsilon*k)) through fourth order."""

    z1: sp.Expr
    z2: sp.Expr
    z3: sp.Expr
    z4: sp.Expr


def d2q5_matrices():
    """Return the D2Q5 moment matrix and lattice velocities."""

    moment = sp.Matrix(
        [
            [1, 1, 1, 1, 1],
            [0, 1, 0, -1, 0],
            [0, 0, 1, 0, -1],
            [-4, 1, 1, 1, 1],
            [0, 1, -1, 1, -1],
        ]
    )
    velocities = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)]
    return moment, velocities


def hydrodynamic_log_expansion() -> tuple[Expansion, dict[str, sp.Symbol]]:
    """Derive the hydrodynamic Fourier eigenvalue by matrix perturbation."""

    kx, ky = sp.symbols("kx ky", real=True)
    ux, uy, a = sp.symbols("ux uy a", real=True)
    q1, q3, q4 = sp.symbols("q1 q3 q4", nonzero=True, real=True)

    moment, velocities = d2q5_matrices()
    inv_moment = moment.inv()

    # Collision map n* = C n.  The equilibrium vector is
    # [T, ux*T, uy*T, a*T, 0]^T.
    collision = sp.zeros(5)
    collision[0, 0] = 1
    collision[1, 0] = q1 * ux
    collision[1, 1] = 1 - q1
    collision[2, 0] = q1 * uy
    collision[2, 2] = 1 - q1
    collision[3, 0] = q3 * a
    collision[3, 3] = 1 - q3
    collision[4, 4] = 1 - q4

    # B_r is the epsilon^r coefficient of
    # M diag(exp(-i epsilon k.c_i)) M^{-1} C.
    streaming_coefficients = []
    for order in range(5):
        diagonal = sp.diag(
            *[
                (-I * (cx * kx + cy * ky)) ** order / sp.factorial(order)
                for cx, cy in velocities
            ]
        )
        streaming_coefficients.append(
            moment * diagonal * inv_moment * collision
        )

    eigenvectors = [sp.Matrix([1, ux, uy, a, 0])]
    eigenvalues = [sp.Integer(1)]

    # Normalization: the conserved component of every correction is zero.
    # The fifth unknown at each order is the eigenvalue coefficient.
    solve_matrix = (
        (streaming_coefficients[0] - sp.eye(5))[:, 1:5]
        .row_join(-eigenvectors[0])
    )
    inverse_solve_matrix = solve_matrix.inv()

    for order in range(1, 5):
        rhs = sp.zeros(5, 1)
        for r in range(1, order):
            rhs += eigenvalues[r] * eigenvectors[order - r]
        for r in range(1, order + 1):
            rhs -= streaming_coefficients[r] * eigenvectors[order - r]

        solution = inverse_solve_matrix * rhs
        eigenvectors.append(sp.Matrix([0, *solution[:4]]))
        eigenvalues.append(sp.cancel(solution[4]))

    l1, l2, l3, l4 = eigenvalues[1:]
    z1 = l1
    z2 = sp.expand(l2 - l1**2 / 2)
    z3 = sp.expand(l3 - l1 * l2 + l1**3 / 3)
    z4 = sp.expand(
        l4 - l1 * l3 - l2**2 / 2 + l1**2 * l2 - l1**4 / 4
    )

    return (
        Expansion(z1=z1, z2=z2, z3=z3, z4=z4),
        {
            "kx": kx,
            "ky": ky,
            "ux": ux,
            "uy": uy,
            "a": a,
            "q1": q1,
            "q3": q3,
            "q4": q4,
        },
    )


def verify_source_collapse() -> None:
    """Check q_eff and its Hénon-parameter mapping."""

    q0, chi = sp.symbols("q0 chi", nonzero=True)
    q_eff = 2 * q0 / (2 * (1 - chi) + chi * q0)
    sigma0 = 1 / q0 - sp.Rational(1, 2)
    sigma_eff = 1 / q_eff - sp.Rational(1, 2)
    assert sp.simplify(sigma_eff - (1 - chi) * sigma0) == 0

    # R_alpha = -chi*q_eff*j_alpha^neq makes the source-augmented
    # flux collision exactly j* = (1-q_eff)j + q_eff*j_eq.
    source_ratio = -chi * q_eff
    q_from_collision = q0 - (1 - q0 / 2) * source_ratio
    assert sp.simplify(q_from_collision - q_eff) == 0


def verify_pure_diffusion(expansion: Expansion, symbols) -> None:
    """Check the Dubois-Lallemand pure-diffusion coefficients."""

    kx = symbols["kx"]
    ky = symbols["ky"]
    ux = symbols["ux"]
    uy = symbols["uy"]
    a = symbols["a"]
    q1 = symbols["q1"]
    q3 = symbols["q3"]
    q4 = symbols["q4"]
    h1, h3, h4 = sp.symbols("h1 h3 h4", real=True)

    relaxation_to_henon = {
        q1: 1 / (h1 + sp.Rational(1, 2)),
        q3: 1 / (h3 + sp.Rational(1, 2)),
        q4: 1 / (h4 + sp.Rational(1, 2)),
    }
    z2 = sp.cancel(
        expansion.z2.subs(relaxation_to_henon).subs({ux: 0, uy: 0})
    )
    expected_z2 = -h1 * (a + 4) * (kx**2 + ky**2) / 10
    assert sp.simplify(z2 - expected_z2) == 0

    k40 = (
        8
        - 3 * a
        + 12 * (a + 4) * h1**2
        - 12 * (1 - a) * h1 * h3
        - 60 * h1 * h4
    )
    k22 = (
        -6 * (a + 4)
        + 24 * (a + 4) * h1**2
        - 24 * (1 - a) * h1 * h3
        + 120 * h1 * h4
    )
    expected_z4 = -h1 * (a + 4) * (
        k40 * (kx**4 + ky**4) + k22 * kx**2 * ky**2
    ) / 1200
    z4 = sp.cancel(
        expansion.z4.subs(relaxation_to_henon).subs({ux: 0, uy: 0})
    )
    assert sp.simplify(z4 - expected_z4) == 0

    # Isotropy and exact cancellation conditions.
    assert sp.simplify(
        k22 - 2 * k40 - 40 * (6 * h1 * h4 - 1)
    ) == 0
    h4_quartic = 1 / (6 * h1)
    h3_quartic = (
        h1 * (a + 4) / (1 - a)
        - (2 + 3 * a) / (12 * h1 * (1 - a))
    )
    assert sp.simplify(k40.subs({h4: h4_quartic, h3: h3_quartic})) == 0
    assert sp.simplify(k22.subs({h4: h4_quartic, h3: h3_quartic})) == 0


def print_henon_coefficients(expansion: Expansion, symbols) -> None:
    """Print the constant-advection coefficients in compact polynomial form."""

    kx = symbols["kx"]
    ky = symbols["ky"]
    q1 = symbols["q1"]
    q3 = symbols["q3"]
    q4 = symbols["q4"]
    h1, h3, h4 = sp.symbols("h r t", real=True)
    substitutions = {
        q1: 1 / (h1 + sp.Rational(1, 2)),
        q3: 1 / (h3 + sp.Rational(1, 2)),
        q4: 1 / (h4 + sp.Rational(1, 2)),
    }

    print("log(lambda_h) coefficients in Hénon parameters")
    for name, expression, degree, imaginary in (
        ("z1", expansion.z1, 1, True),
        ("z2", expansion.z2, 2, False),
        ("z3/i", expansion.z3 / I, 3, False),
        ("z4", expansion.z4, 4, False),
    ):
        converted = sp.cancel(expression.subs(substitutions))
        polynomial = sp.Poly(sp.expand(converted), kx, ky)
        print(f"\n{name}:")
        for power_x in range(degree, -1, -1):
            power_y = degree - power_x
            monomial = kx**power_x * ky**power_y
            coefficient = sp.factor(polynomial.coeff_monomial(monomial))
            if coefficient != 0:
                print(f"  kx^{power_x} ky^{power_y}: {coefficient}")


def run_scan(
    a_values: list[float],
    speed_values: list[float],
    sigma_values: list[float],
    wave_points: int,
    q_points: int,
) -> None:
    """Run a configurable coarse von Neumann scan for target-case screening."""

    import numpy as np

    moment = np.array(
        [
            [1, 1, 1, 1, 1],
            [0, 1, 0, -1, 0],
            [0, 0, 1, 0, -1],
            [-4, 1, 1, 1, 1],
            [0, 1, -1, 1, -1],
        ],
        dtype=float,
    )
    inv_moment = np.linalg.inv(moment)
    velocities = np.array(
        [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)],
        dtype=float,
    )
    wave_numbers = np.linspace(0, math.pi, wave_points)
    wave_vectors = np.array(
        [(kx, ky) for kx in wave_numbers for ky in wave_numbers]
    )
    streaming = np.exp(-1j * (wave_vectors @ velocities.T))
    q_grid = np.linspace(0.05, 1.95, q_points)

    def spectral_radius(a, ux, uy, sigma_eff, q3, q4):
        q1 = 1 / (sigma_eff + 0.5)
        collision = np.zeros((5, 5))
        collision[0, 0] = 1
        collision[1, 0] = q1 * ux
        collision[1, 1] = 1 - q1
        collision[2, 0] = q1 * uy
        collision[2, 2] = 1 - q1
        collision[3, 0] = q3 * a
        collision[3, 3] = 1 - q3
        collision[4, 4] = 1 - q4
        distribution_collision = inv_moment @ collision @ moment
        amplification = (
            streaming[:, :, None] * distribution_collision[None, :, :]
        )
        return float(
            np.max(np.abs(np.linalg.eigvals(amplification)))
        )

    for a in a_values:
        for speed in speed_values:
            for sigma_eff in sigma_values:
                stable = []
                for q3 in q_grid:
                    for q4 in q_grid:
                        axial = spectral_radius(
                            a, speed, 0, sigma_eff, q3, q4
                        )
                        diagonal = spectral_radius(
                            a,
                            speed / math.sqrt(2),
                            speed / math.sqrt(2),
                            sigma_eff,
                            q3,
                            q4,
                        )
                        if max(axial, diagonal) <= 1 + 2.0e-10:
                            stable.append((q3, q4))

                candidates = {}
                for label, q3, q4 in (
                    ("q3=q4=1", 1.0, 1.0),
                    ("q3=q4=1.25", 1.25, 1.25),
                    ("q3=q4=1.5", 1.5, 1.5),
                    ("q3=q4=1.8", 1.8, 1.8),
                ):
                    candidates[label] = max(
                        spectral_radius(a, speed, 0, sigma_eff, q3, q4),
                        spectral_radius(
                            a,
                            speed / math.sqrt(2),
                            speed / math.sqrt(2),
                            sigma_eff,
                            q3,
                            q4,
                        ),
                    )

                sigma4_quartic = 1.0 / (6.0 * sigma_eff)
                q4_quartic = 1.0 / (0.5 + sigma4_quartic)
                if abs(1.0 - a) > 1.0e-14:
                    sigma3_quartic = (
                        sigma_eff * (a + 4.0) / (1.0 - a)
                        - (2.0 + 3.0 * a)
                        / (12.0 * sigma_eff * (1.0 - a))
                    )
                    q3_quartic = 1.0 / (0.5 + sigma3_quartic)
                else:
                    q3_quartic = math.nan

                print(
                    {
                        "a": a,
                        "|u|": speed,
                        "sigma_eff": sigma_eff,
                        "q_eff": 1 / (sigma_eff + 0.5),
                        "stable_count": len(stable),
                        "total": len(q_grid) ** 2,
                        "stable_q3_envelope": (
                            min(x[0] for x in stable),
                            max(x[0] for x in stable),
                        )
                        if stable
                        else None,
                        "stable_q4_envelope": (
                            min(x[1] for x in stable),
                            max(x[1] for x in stable),
                        )
                        if stable
                        else None,
                        "candidate_spectral_radius": candidates,
                        "source_free_quartic_q3_q4": (
                            q3_quartic,
                            q4_quartic,
                        ),
                    }
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scan",
        action="store_true",
        help="also run the coarse von Neumann parameter scan",
    )
    parser.add_argument(
        "--scan-a",
        type=float,
        action="append",
        help="D2Q5 a value; repeat for multiple values",
    )
    parser.add_argument(
        "--scan-speed",
        type=float,
        action="append",
        help="constant advection speed magnitude; repeat as needed",
    )
    parser.add_argument(
        "--scan-sigma",
        type=float,
        action="append",
        help="effective flux Henon parameter; repeat as needed",
    )
    parser.add_argument("--wave-points", type=int, default=25)
    parser.add_argument("--q-points", type=int, default=39)
    parser.add_argument(
        "--quiet-coefficients",
        action="store_true",
        help="verify the symbolic identities without printing every expansion coefficient",
    )
    args = parser.parse_args()

    verify_source_collapse()
    expansion, symbols = hydrodynamic_log_expansion()
    verify_pure_diffusion(expansion, symbols)
    print("All symbolic identity checks passed.")
    if not args.quiet_coefficients:
        print_henon_coefficients(expansion, symbols)
    if args.scan:
        run_scan(
            a_values=args.scan_a or [0.0],
            speed_values=args.scan_speed or [0.1],
            sigma_values=args.scan_sigma or [0.1, 0.01, 0.001],
            wave_points=args.wave_points,
            q_points=args.q_points,
        )


if __name__ == "__main__":
    main()
