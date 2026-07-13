"""Exact D2Q lattice definitions and moment operators."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from sympy import Expr, Integer, Matrix, Rational, eye, simplify, sympify, zeros


Velocity = tuple[int, int]


@dataclass(frozen=True)
class Lattice:
    """A two-dimensional velocity lattice with exact symbolic constants."""

    velocities: tuple[Velocity, ...]
    weights: tuple[Expr, ...]
    lambda_t: tuple[Expr, ...]
    opposite: tuple[int, ...]
    cs2: Expr


def d2q9() -> Lattice:
    """Return the standard D2Q9 lattice used by the LBM-CDE model."""

    return Lattice(
        velocities=(
            (0, 0),
            (1, 0),
            (0, 1),
            (-1, 0),
            (0, -1),
            (1, 1),
            (-1, 1),
            (-1, -1),
            (1, -1),
        ),
        weights=(
            Rational(4, 9),
            Rational(1, 9),
            Rational(1, 9),
            Rational(1, 9),
            Rational(1, 9),
            Rational(1, 36),
            Rational(1, 36),
            Rational(1, 36),
            Rational(1, 36),
        ),
        lambda_t=(
            Rational(-5, 9),
            Rational(1, 9),
            Rational(1, 9),
            Rational(1, 9),
            Rational(1, 9),
            Rational(1, 36),
            Rational(1, 36),
            Rational(1, 36),
            Rational(1, 36),
        ),
        opposite=(0, 3, 4, 1, 2, 7, 8, 5, 6),
        cs2=Rational(1, 3),
    )


def d2q5() -> Lattice:
    """Return the standard isotropic-second-moment D2Q5 verifier lattice."""

    return Lattice(
        velocities=((0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)),
        weights=(
            Rational(1, 3),
            Rational(1, 6),
            Rational(1, 6),
            Rational(1, 6),
            Rational(1, 6),
        ),
        lambda_t=(),
        opposite=(0, 3, 4, 1, 2),
        cs2=Rational(1, 3),
    )


def parity_projectors(lattice: Lattice) -> tuple[Matrix, Matrix]:
    """Return exact projectors onto opposite-pair even and odd populations."""

    size = len(lattice.velocities)
    reversal = zeros(size)
    for index, opposite in enumerate(lattice.opposite):
        reversal[index, opposite] = 1
    return (eye(size) + reversal) / 2, (eye(size) - reversal) / 2


def raw_moment(
    coefficients: Sequence[Expr],
    velocities: Sequence[Velocity],
    powers: tuple[int, int],
) -> Expr:
    """Return ``sum_i coefficients_i c_ix**px c_iy**py`` exactly."""

    power_x, power_y = powers
    return simplify(
        sum(
            (
                sympify(coefficient)
                * Integer(velocity_x) ** power_x
                * Integer(velocity_y) ** power_y
            )
            for coefficient, (velocity_x, velocity_y) in zip(
                coefficients, velocities, strict=True
            )
        )
    )


def hermite_moment(
    coefficients: Sequence[Expr],
    velocities: Sequence[Velocity],
    axes: tuple[int, ...],
    cs2: Expr,
) -> Expr:
    """Return a rank-zero through rank-three Cartesian Hermite moment."""

    return simplify(
        sum(
            sympify(coefficient)
            * _hermite_polynomial(velocity, axes, sympify(cs2))
            for coefficient, velocity in zip(
                coefficients, velocities, strict=True
            )
        )
    )


def _hermite_polynomial(
    velocity: Iterable[int], axes: tuple[int, ...], cs2: Expr
) -> Expr:
    components = tuple(Integer(component) for component in velocity)
    rank = len(axes)
    if rank == 0:
        return Integer(1)
    if rank == 1:
        return components[axes[0]]
    if rank == 2:
        alpha, beta = axes
        return components[alpha] * components[beta] - cs2 * Integer(alpha == beta)
    if rank == 3:
        alpha, beta, gamma = axes
        return (
            components[alpha] * components[beta] * components[gamma]
            - cs2
            * (
                components[alpha] * Integer(beta == gamma)
                + components[beta] * Integer(alpha == gamma)
                + components[gamma] * Integer(alpha == beta)
            )
        )
    raise ValueError("Hermite moments are implemented only through rank three")
