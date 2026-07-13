"""Exact transformed-distribution TRT and BGK collision operators."""

from __future__ import annotations

from typing import Sequence

from sympy import Expr, Matrix, simplify, sympify

from tools.derivation.lattice import Lattice, parity_projectors, raw_moment


def trt_collision(
    lattice: Lattice,
    *,
    h_tilde: Sequence[Expr],
    h_eq: Sequence[Expr],
    source_plus: Sequence[Expr],
    source_minus: Sequence[Expr],
    s_plus: Expr,
    s_minus: Expr,
    dt: Expr = 1,
) -> Matrix:
    """Apply the operator-trapezoidal TRT collision before streaming."""

    h_tilde = _column(h_tilde)
    h_eq = _column(h_eq)
    source_plus = _column(source_plus)
    source_minus = _column(source_minus)
    s_plus = sympify(s_plus)
    s_minus = sympify(s_minus)
    dt = sympify(dt)
    p_plus, p_minus = parity_projectors(lattice)
    nonequilibrium = h_tilde - h_eq
    return (
        h_tilde
        - s_plus * p_plus * nonequilibrium
        - s_minus * p_minus * nonequilibrium
        + dt
        * (
            (1 - s_plus / 2) * source_plus
            + (1 - s_minus / 2) * source_minus
        )
    ).applyfunc(simplify)


def bgk_collision(
    *,
    h_tilde: Sequence[Expr],
    h_eq: Sequence[Expr],
    source: Sequence[Expr],
    rate: Expr,
    dt: Expr = 1,
) -> Matrix:
    """Apply the equal-rate transformed BGK collision componentwise."""

    h_tilde = _column(h_tilde)
    h_eq = _column(h_eq)
    source = _column(source)
    rate = sympify(rate)
    dt = sympify(dt)
    return (
        h_tilde
        - rate * (h_tilde - h_eq)
        + dt * (1 - rate / 2) * source
    ).applyfunc(simplify)


def reconstruct_momentum(
    lattice: Lattice,
    *,
    f_tilde: Sequence[Expr],
    force: Sequence[Expr],
    dt: Expr = 1,
) -> Matrix:
    """Return ``rho0*u`` from transformed populations and the half force."""

    f_tilde = _column(f_tilde)
    force = _column(force)
    dt = sympify(dt)
    kinetic_momentum = Matrix(
        [
            raw_moment(f_tilde, lattice.velocities, (1, 0)),
            raw_moment(f_tilde, lattice.velocities, (0, 1)),
        ]
    )
    return (kinetic_momentum + dt * force / 2).applyfunc(simplify)


def reconstruct_scalar(
    *,
    g_tilde: Sequence[Expr],
    heat_source: Expr,
    dt: Expr = 1,
) -> Expr:
    """Return the scalar from transformed populations and the half source."""

    g_tilde = _column(g_tilde)
    return simplify(sum(g_tilde) + sympify(dt) * sympify(heat_source) / 2)


def _column(values: Sequence[Expr]) -> Matrix:
    column = Matrix(values)
    if column.cols != 1:
        column = column.reshape(len(column), 1)
    return column.applyfunc(sympify)
