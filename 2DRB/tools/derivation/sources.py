"""Exact D2Q source terms and their parity-resolved raw moments."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from sympy import Expr, Matrix, eye, simplify, sympify

from tools.derivation.lattice import Lattice, parity_projectors, raw_moment


Power = tuple[int, int]
MomentTable = Mapping[Power, Expr]


@dataclass(frozen=True)
class SourceDecomposition:
    """A source vector, its exact parity projections, and raw moments."""

    raw: Matrix
    plus: Matrix
    minus: Matrix
    moment_tables: Mapping[str, MomentTable]


def source_moment_table(
    coefficients: Sequence[Expr],
    lattice: Lattice,
    max_order: int = 4,
) -> dict[Power, Expr]:
    """Return every two-dimensional raw moment through ``max_order``."""

    return {
        (power_x, total_order - power_x): raw_moment(
            coefficients,
            lattice.velocities,
            (power_x, total_order - power_x),
        )
        for total_order in range(max_order + 1)
        for power_x in range(total_order, -1, -1)
    }


def flow_source(
    lattice: Lattice,
    *,
    velocity: Sequence[Expr],
    force: Sequence[Expr],
    rho0: Expr,
    strain: Matrix,
    chi_s: Expr,
    chi_b: Expr,
) -> SourceDecomposition:
    """Return the Eq. (13) flow source split into even and odd parts.

    The second Hermite tensor is dimensionless,
    ``H2_i = c_i c_i / cs2 - I``.  Consequently the positive contraction
    ``rho0 * H2_i:A`` uses
    ``A = chi_s*S + (chi_b-chi_s)*tr(S)*I/D`` without another ``1/cs2``.
    The caller supplies the symmetric strain tensor from Eq. (2).
    """

    velocity_vector = _column(velocity)
    force_vector = _column(force)
    strain_tensor = Matrix(strain).applyfunc(sympify)
    dimension = len(velocity_vector)
    if dimension != 2 or strain_tensor.shape != (2, 2):
        raise ValueError("flow_source requires two-dimensional vectors and strain")

    rho0 = sympify(rho0)
    chi_s = sympify(chi_s)
    chi_b = sympify(chi_b)
    cs2 = sympify(lattice.cs2)
    correction_tensor = (
        chi_s * strain_tensor
        + (chi_b - chi_s)
        * strain_tensor.trace()
        * eye(dimension)
        / dimension
    )
    velocity_force = (velocity_vector.T * force_vector)[0]

    populations = []
    for weight, velocity_i in zip(
        lattice.weights, lattice.velocities, strict=True
    ):
        discrete_velocity = Matrix(velocity_i)
        discrete_force = (discrete_velocity.T * force_vector)[0]
        discrete_flow = (discrete_velocity.T * velocity_vector)[0]
        hermite_two = discrete_velocity * discrete_velocity.T / cs2 - eye(dimension)
        contraction = sum(
            hermite_two[alpha, beta] * correction_tensor[alpha, beta]
            for alpha in range(dimension)
            for beta in range(dimension)
        )
        populations.append(
            simplify(
                weight
                * (
                    discrete_force / cs2
                    + discrete_flow * discrete_force / cs2**2
                    - velocity_force / cs2
                    + rho0 * contraction
                )
            )
        )

    return _decompose(Matrix(populations), lattice)


def scalar_source(
    lattice: Lattice,
    *,
    velocity: Sequence[Expr],
    force: Sequence[Expr],
    pressure: Expr,
    temperature: Expr,
    grad_temperature: Sequence[Expr],
    rho0: Expr,
    heat_source: Expr,
    chi_kappa: Expr,
) -> SourceDecomposition:
    """Return the Eq. (24) scalar source split into even and odd parts."""

    velocity_vector = _column(velocity)
    force_vector = _column(force)
    gradient = _column(grad_temperature)
    if len(velocity_vector) != 2 or len(force_vector) != 2 or len(gradient) != 2:
        raise ValueError("scalar_source requires two-dimensional vectors")

    pressure = sympify(pressure)
    temperature = sympify(temperature)
    rho0 = sympify(rho0)
    heat_source = sympify(heat_source)
    chi_kappa = sympify(chi_kappa)
    cs2 = sympify(lattice.cs2)
    coupling = pressure * gradient + temperature * force_vector

    populations = []
    for weight, velocity_i in zip(
        lattice.weights, lattice.velocities, strict=True
    ):
        discrete_velocity = Matrix(velocity_i)
        discrete_coupling = (discrete_velocity.T * coupling)[0]
        discrete_flow = (discrete_velocity.T * velocity_vector)[0]
        discrete_gradient = (discrete_velocity.T * gradient)[0]
        populations.append(
            simplify(
                weight
                * (
                    discrete_coupling / (rho0 * cs2)
                    + heat_source * (1 + discrete_flow / cs2)
                    + chi_kappa * discrete_gradient
                )
            )
        )

    return _decompose(Matrix(populations), lattice)


def _decompose(raw: Matrix, lattice: Lattice) -> SourceDecomposition:
    p_plus, p_minus = parity_projectors(lattice)
    plus = (p_plus * raw).applyfunc(simplify)
    minus = (p_minus * raw).applyfunc(simplify)
    raw = raw.applyfunc(simplify)
    return SourceDecomposition(
        raw=raw,
        plus=plus,
        minus=minus,
        moment_tables={
            "raw": source_moment_table(raw, lattice),
            "plus": source_moment_table(plus, lattice),
            "minus": source_moment_table(minus, lattice),
        },
    )


def _column(values: Sequence[Expr]) -> Matrix:
    column = Matrix(values)
    if column.cols != 1:
        column = column.reshape(len(column), 1)
    return column.applyfunc(sympify)
