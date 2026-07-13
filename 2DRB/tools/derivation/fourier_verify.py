"""High-precision directional checks for D2Q9 temperature amplification."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from mpmath import mp
from sympy import Expr, Matrix, N, Rational, sqrt, sympify

from tools.derivation.d2q9_temperature import (
    build_d2q9_temperature_model,
    d2q9_amplification,
)


@dataclass(frozen=True)
class DirectionalFourierFit:
    """Shrinking-wave-number estimates from direct matrix eigenvalues."""

    precision: int
    wave_numbers: tuple
    diffusion: object
    axis_gamma: tuple
    equal_diagonal_gamma: tuple
    axis_quartic: tuple
    equal_diagonal_quartic: tuple
    axis_orders: tuple
    equal_diagonal_orders: tuple


def high_precision_directional_fit(
    *,
    case: str,
    pi: Expr,
    chi_kappa: Expr,
    sigma_flux: Expr,
    sigma_odd_ghost: Expr,
    sigma_even: Expr,
    dt: Expr = 1,
    precision: int = 80,
    wave_numbers: Sequence[Expr] = (
        Rational(1, 50),
        Rational(1, 100),
        Rational(1, 200),
        Rational(1, 400),
    ),
) -> DirectionalFourierFit:
    """Fit axis/equal-magnitude-diagonal residual orders at high precision."""

    if precision < 80:
        raise ValueError("at least 80 digits are required")
    if len(wave_numbers) < 3:
        raise ValueError("at least three shrinking wave numbers are required")
    exact_dt = sympify(dt)
    exact_sigma_flux = sympify(sigma_flux)
    model = build_d2q9_temperature_model(
        pi=pi,
        chi_kappa=chi_kappa,
        sigma_odd=sigma_odd_ghost,
        sigma_even=sigma_even,
    )
    diffusivity_factor = {
        "baseline": model.lattice.cs2,
        "external": model.b,
        "feedback": model.a,
    }.get(case)
    if diffusivity_factor is None:
        raise ValueError("case must be 'baseline', 'external', or 'feedback'")
    exact_diffusion = diffusivity_factor * exact_sigma_flux * exact_dt

    with mp.workdps(precision):
        q_values = tuple(_to_mpmath(sympify(q), precision).real for q in wave_numbers)
        if any(q <= 0 for q in q_values):
            raise ValueError("wave numbers must be positive")
        if any(q_values[index + 1] >= q_values[index] for index in range(len(q_values) - 1)):
            raise ValueError("wave numbers must form a strictly shrinking sequence")
        diffusion = _to_mpmath(exact_diffusion, precision).real
        dt_value = _to_mpmath(exact_dt, precision).real
        axis_gamma = tuple(
            _hydrodynamic_decay(
                case=case,
                pi=pi,
                chi_kappa=chi_kappa,
                sigma_flux=sigma_flux,
                sigma_odd_ghost=sigma_odd_ghost,
                sigma_even=sigma_even,
                dt=exact_dt,
                kx=sympify(q),
                ky=0,
                precision=precision,
                dt_value=dt_value,
            )
            for q in wave_numbers
        )
        equal_diagonal_gamma = tuple(
            _hydrodynamic_decay(
                case=case,
                pi=pi,
                chi_kappa=chi_kappa,
                sigma_flux=sigma_flux,
                sigma_odd_ghost=sigma_odd_ghost,
                sigma_even=sigma_even,
                dt=exact_dt,
                kx=sympify(q) / sqrt(2),
                ky=sympify(q) / sqrt(2),
                precision=precision,
                dt_value=dt_value,
            )
            for q in wave_numbers
        )
        axis_residuals = tuple(
            gamma - diffusion * q**2
            for gamma, q in zip(axis_gamma, q_values, strict=True)
        )
        equal_diagonal_residuals = tuple(
            gamma - diffusion * q**2
            for gamma, q in zip(
                equal_diagonal_gamma, q_values, strict=True
            )
        )
        axis_quartic = tuple(
            residual / q**4
            for residual, q in zip(axis_residuals, q_values, strict=True)
        )
        equal_diagonal_quartic = tuple(
            residual / q**4
            for residual, q in zip(
                equal_diagonal_residuals, q_values, strict=True
            )
        )
        axis_orders = _observed_orders(axis_residuals, q_values)
        equal_diagonal_orders = _observed_orders(
            equal_diagonal_residuals,
            q_values,
        )
        return DirectionalFourierFit(
            precision=precision,
            wave_numbers=q_values,
            diffusion=diffusion,
            axis_gamma=axis_gamma,
            equal_diagonal_gamma=equal_diagonal_gamma,
            axis_quartic=axis_quartic,
            equal_diagonal_quartic=equal_diagonal_quartic,
            axis_orders=axis_orders,
            equal_diagonal_orders=equal_diagonal_orders,
        )


def _hydrodynamic_decay(
    *,
    case: str,
    pi: Expr,
    chi_kappa: Expr,
    sigma_flux: Expr,
    sigma_odd_ghost: Expr,
    sigma_even: Expr,
    dt: Expr,
    kx: Expr,
    ky: Expr,
    precision: int,
    dt_value,
):
    exact_matrix = d2q9_amplification(
        case=case,
        pi=pi,
        chi_kappa=chi_kappa,
        sigma_flux=sigma_flux,
        sigma_odd_ghost=sigma_odd_ghost,
        sigma_even=sigma_even,
        dt=dt,
        kx=kx,
        ky=ky,
    )
    numeric_matrix = _matrix_to_mpmath(exact_matrix, precision)
    eigenvalues = mp.eig(numeric_matrix, left=False, right=False)
    hydrodynamic_root = min(eigenvalues, key=lambda value: abs(value - 1))
    return -mp.log(hydrodynamic_root) / dt_value


def _matrix_to_mpmath(matrix: Matrix, precision: int):
    return mp.matrix(
        [
            [_to_mpmath(matrix[row, column], precision) for column in range(matrix.cols)]
            for row in range(matrix.rows)
        ]
    )


def _to_mpmath(expression: Expr, precision: int):
    real, imaginary = N(expression, precision).as_real_imag()
    return mp.mpc(str(real), str(imaginary))


def _observed_orders(residuals: tuple, wave_numbers: tuple) -> tuple:
    return tuple(
        mp.log(abs(residuals[index] / residuals[index + 1]))
        / mp.log(wave_numbers[index] / wave_numbers[index + 1])
        for index in range(len(residuals) - 1)
    )
