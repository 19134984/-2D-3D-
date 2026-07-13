"""Parameterized D2Q5 pure-diffusion reference generator."""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial

from sympy import (
    Expr,
    I,
    Matrix,
    Rational,
    diag,
    exp,
    eye,
    expand,
    simplify,
    symbols,
    sympify,
    zeros,
)

from tools.derivation.series import truncate_in_scale


@dataclass(frozen=True)
class ParameterizedD2Q5:
    """Dubois Appendix-B D2Q5 model with a free equilibrium parameter."""

    alpha: Expr
    s1: Expr
    s3: Expr
    s4: Expr
    lattice_speed: Expr
    velocities: tuple[tuple[int, int], ...]
    moment_matrix: Matrix
    inverse_moment_matrix: Matrix
    collision_matrix: Matrix
    equilibrium_moments: Matrix
    equilibrium_weights: tuple[Expr, ...]
    equilibrium_second_moment: Expr


@dataclass(frozen=True)
class D2Q5EquivalentCoefficients:
    """Second- and optional fourth-order coefficients of the pure thermal mode."""

    diffusion: Expr
    cubic: Expr | None = None
    fourth_axis: Expr | None = None
    fourth_mixed: Expr | None = None
    kappa40: Expr | None = None
    kappa22: Expr | None = None


def parameterized_d2q5(
    alpha: Expr,
    s1: Expr,
    s3: Expr,
    s4: Expr,
    *,
    lattice_speed: Expr = 1,
) -> ParameterizedD2Q5:
    """Return the parameterized D2Q5 reference model."""

    alpha = sympify(alpha)
    s1 = sympify(s1)
    s3 = sympify(s3)
    s4 = sympify(s4)
    lattice_speed = sympify(lattice_speed)
    moment_matrix = Matrix(
        [
            [1, 1, 1, 1, 1],
            [0, lattice_speed, 0, -lattice_speed, 0],
            [0, 0, lattice_speed, 0, -lattice_speed],
            [-4, 1, 1, 1, 1],
            [0, 1, -1, 1, -1],
        ]
    )
    inverse_moment_matrix = moment_matrix.inv()
    collision_matrix = Matrix(
        [
            [1, 0, 0, 0, 0],
            [0, 1 - s1, 0, 0, 0],
            [0, 0, 1 - s1, 0, 0],
            [alpha * s3, 0, 0, 1 - s3, 0],
            [0, 0, 0, 0, 1 - s4],
        ]
    )
    equilibrium_moments = Matrix([1, 0, 0, alpha, 0])
    equilibrium_populations = inverse_moment_matrix * equilibrium_moments
    equilibrium_weights = tuple(
        simplify(value) for value in equilibrium_populations
    )

    return ParameterizedD2Q5(
        alpha=alpha,
        s1=s1,
        s3=s3,
        s4=s4,
        lattice_speed=lattice_speed,
        velocities=((0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)),
        moment_matrix=moment_matrix,
        inverse_moment_matrix=inverse_moment_matrix,
        collision_matrix=collision_matrix,
        equilibrium_moments=equilibrium_moments,
        equilibrium_weights=equilibrium_weights,
        equilibrium_second_moment=simplify((4 + alpha) / 10),
    )


def amplification_matrix(
    model: ParameterizedD2Q5,
    kx: Expr,
    ky: Expr,
    *,
    time_step: Expr = 1,
) -> Matrix:
    """Return Dubois' positive-phase population amplification matrix."""

    kx = sympify(kx)
    ky = sympify(ky)
    time_step = sympify(time_step)
    phase_x = kx * model.lattice_speed * time_step
    phase_y = ky * model.lattice_speed * time_step
    streaming = diag(
        1,
        exp(I * phase_x),
        exp(I * phase_y),
        exp(-I * phase_x),
        exp(-I * phase_y),
    )
    return (
        streaming
        * model.inverse_moment_matrix
        * model.collision_matrix
        * model.moment_matrix
    )


def amplification_route(
    alpha: Expr,
    sigma1: Expr,
    sigma3: Expr,
    sigma4: Expr,
    *,
    lattice_speed: Expr = 1,
    time_step: Expr = 1,
    order: int = 2,
) -> D2Q5EquivalentCoefficients:
    """Generate coefficients by a formal hydrodynamic amplification series."""

    if order not in {2, 4}:
        raise ValueError("the amplification route supports order=2 or order=4")
    alpha = sympify(alpha)
    sigma1 = sympify(sigma1)
    lattice_speed = sympify(lattice_speed)
    time_step = sympify(time_step)
    model = _model_from_shifts(
        alpha,
        sigma1,
        sigma3,
        sigma4,
        lattice_speed=1,
    )
    kx, ky = symbols("kx ky")
    phase_directions = (
        0,
        kx,
        ky,
        -kx,
        -ky,
    )
    population_collision = (
        model.inverse_moment_matrix
        * model.collision_matrix
        * model.moment_matrix
    )
    operator_terms = tuple(
        diag(
            *(
                (I * direction) ** degree / factorial(degree)
                for direction in phase_directions
            )
        )
        * population_collision
        for degree in range(order + 1)
    )
    base_vector = Matrix(model.equilibrium_weights)
    conserved_row = Matrix([[1, 1, 1, 1, 1]])
    augmented = (
        (operator_terms[0] - eye(5)).row_join(-base_vector)
    ).col_join(conserved_row.row_join(Matrix([[0]])))
    augmented_inverse = augmented.inv()
    vector_terms = [base_vector]
    value_terms: list[Expr] = [sympify(1)]

    for degree in range(1, order + 1):
        right_hand_side = -sum(
            (
                operator_terms[index] * vector_terms[degree - index]
                for index in range(1, degree + 1)
            ),
            zeros(5, 1),
        )
        right_hand_side += sum(
            (
                value_terms[index] * vector_terms[degree - index]
                for index in range(1, degree)
            ),
            zeros(5, 1),
        )
        solution = augmented_inverse * right_hand_side.col_join(Matrix([0]))
        vector_terms.append(solution[:5, :].applyfunc(simplify))
        value_terms.append(simplify(solution[5]))

    decay_degree_two = simplify(
        -(value_terms[2] - value_terms[1] ** 2 / 2)
    )
    diffusion_dimensionless = simplify(expand(decay_degree_two).coeff(kx, 2))
    diffusion = simplify(
        diffusion_dimensionless * lattice_speed**2 * time_step
    )
    if order == 2:
        return D2Q5EquivalentCoefficients(diffusion=diffusion)

    decay_degree_three = simplify(
        -(
            value_terms[3]
            - value_terms[1] * value_terms[2]
            + value_terms[1] ** 3 / 3
        )
    )
    decay_degree_four = simplify(
        -(
            value_terms[4]
            - value_terms[1] * value_terms[3]
            - value_terms[2] ** 2 / 2
            + value_terms[1] ** 2 * value_terms[2]
            - value_terms[1] ** 4 / 4
        )
    )
    fourth_axis_dimensionless = simplify(expand(decay_degree_four).coeff(kx, 4))
    fourth_mixed_dimensionless = simplify(
        expand(decay_degree_four).coeff(kx, 2).coeff(ky, 2)
    )
    fourth_axis = simplify(
        fourth_axis_dimensionless * lattice_speed**4 * time_step**3
    )
    fourth_mixed = simplify(
        fourth_mixed_dimensionless * lattice_speed**4 * time_step**3
    )
    return D2Q5EquivalentCoefficients(
        diffusion=diffusion,
        cubic=simplify(
            decay_degree_three * lattice_speed**3 * time_step**2
        ),
        fourth_axis=fourth_axis,
        fourth_mixed=fourth_mixed,
        kappa40=simplify(
            1200 * fourth_axis_dimensionless / (sigma1 * (4 + alpha))
        ),
        kappa22=simplify(
            1200 * fourth_mixed_dimensionless / (sigma1 * (4 + alpha))
        ),
    )


def taylor_moment_route(
    alpha: Expr,
    sigma1: Expr,
    sigma3: Expr,
    sigma4: Expr,
    *,
    lattice_speed: Expr = 1,
    time_step: Expr = 1,
    order: int = 2,
) -> D2Q5EquivalentCoefficients:
    """Generate coefficients by physical-space Taylor/moment elimination."""

    if order not in {2, 4}:
        raise ValueError("the Taylor/moment route supports order=2 or order=4")
    alpha = sympify(alpha)
    sigma1 = sympify(sigma1)
    lattice_speed = sympify(lattice_speed)
    time_step = sympify(time_step)
    model = _model_from_shifts(
        alpha,
        sigma1,
        sigma3,
        sigma4,
        lattice_speed=1,
    )
    epsilon, derivative_x, derivative_y = symbols(
        "epsilon derivative_x derivative_y"
    )
    recursion_step = sympify(1)
    base_moments = model.equilibrium_moments
    augmented = (
        (eye(5) - model.collision_matrix).row_join(
            recursion_step * base_moments
        )
    ).col_join(Matrix([[1, 0, 0, 0, 0, 0]]))
    augmented_inverse = augmented.inv()
    moment_terms = [base_moments]
    time_terms: list[Expr] = []
    time_operator = sympify(0)

    for degree in range(1, order + 1):
        residual = _streaming_taylor_residual(
            model,
            moment_terms,
            time_operator,
            epsilon,
            derivative_x,
            derivative_y,
            recursion_step,
            degree,
        )
        coefficient = residual.applyfunc(
            lambda value: expand(value).coeff(epsilon, degree)
        )
        solution = augmented_inverse * (-coefficient).col_join(Matrix([0]))
        moment_terms.append(solution[:5, :].applyfunc(simplify))
        time_term = simplify(solution[5])
        time_terms.append(time_term)
        time_operator += epsilon ** (degree - 1) * time_term

    diffusion_dimensionless = simplify(
        expand(time_terms[1]).coeff(derivative_x, 2)
    )
    diffusion = simplify(
        diffusion_dimensionless * lattice_speed**2 * time_step
    )
    if order == 2:
        return D2Q5EquivalentCoefficients(diffusion=diffusion)

    fourth_axis_dimensionless = simplify(
        -expand(time_terms[3]).coeff(derivative_x, 4)
    )
    fourth_mixed_dimensionless = simplify(
        -expand(time_terms[3]).coeff(derivative_x, 2).coeff(derivative_y, 2)
    )
    fourth_axis = simplify(
        fourth_axis_dimensionless * lattice_speed**4 * time_step**3
    )
    fourth_mixed = simplify(
        fourth_mixed_dimensionless * lattice_speed**4 * time_step**3
    )
    return D2Q5EquivalentCoefficients(
        diffusion=diffusion,
        cubic=simplify(time_terms[2] * lattice_speed**3 * time_step**2),
        fourth_axis=fourth_axis,
        fourth_mixed=fourth_mixed,
        kappa40=simplify(
            1200 * fourth_axis_dimensionless / (sigma1 * (4 + alpha))
        ),
        kappa22=simplify(
            1200 * fourth_mixed_dimensionless / (sigma1 * (4 + alpha))
        ),
    )


def _model_from_shifts(
    alpha: Expr,
    sigma1: Expr,
    sigma3: Expr,
    sigma4: Expr,
    *,
    lattice_speed: Expr,
) -> ParameterizedD2Q5:
    rates = tuple(
        simplify(1 / (sympify(shift) + Rational(1, 2)))
        for shift in (sigma1, sigma3, sigma4)
    )
    return parameterized_d2q5(
        alpha,
        *rates,
        lattice_speed=lattice_speed,
    )


def _streaming_taylor_residual(
    model: ParameterizedD2Q5,
    moment_terms: list[Matrix],
    time_operator: Expr,
    epsilon: Expr,
    derivative_x: Expr,
    derivative_y: Expr,
    time_step: Expr,
    degree: int,
) -> Matrix:
    moments = sum(
        (
            epsilon**index * moment
            for index, moment in enumerate(moment_terms)
        ),
        zeros(5, 1),
    )
    populations = model.inverse_moment_matrix * moments
    postcollision = (
        model.inverse_moment_matrix * model.collision_matrix * moments
    )
    temporal_taylor = sum(
        (
            (epsilon * time_step * time_operator) ** index
            / factorial(index)
            for index in range(degree + 1)
        ),
        sympify(0),
    )
    directions = tuple(
        model.lattice_speed
        * (velocity_x * derivative_x + velocity_y * derivative_y)
        for velocity_x, velocity_y in model.velocities
    )
    residual_populations = Matrix(
        [
            populations[index] * temporal_taylor
            - postcollision[index]
            * sum(
                (
                    (epsilon * time_step * directions[index]) ** power
                    / factorial(power)
                    for power in range(degree + 1)
                ),
                sympify(0),
            )
            for index in range(5)
        ]
    )
    residual_moments = model.moment_matrix * residual_populations
    return residual_moments.applyfunc(
        lambda value: truncate_in_scale(value, epsilon, degree)
    )
