"""Exact frozen-coefficient D2Q9 temperature modified-equation tools."""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial

from sympy import (
    Eq,
    Expr,
    I,
    Matrix,
    Rational,
    diag,
    exp,
    expand,
    eye,
    simplify,
    solve,
    symbols,
    sympify,
    zeros,
)

from tools.derivation.lattice import Lattice, d2q9, parity_projectors
from tools.derivation.series import truncate_in_scale


@dataclass(frozen=True)
class D2Q9TemperatureModel:
    """Raw frozen D2Q9 operators for the three linear temperature schemes."""

    lattice: Lattice
    ell: Matrix
    e: Matrix
    G: Matrix
    p_plus: Matrix
    p_minus: Matrix
    S: Matrix
    H: Matrix
    J: Matrix
    a: Expr
    b: Expr
    d: Expr
    s_odd: Expr
    s_even: Expr
    sigma_odd: Expr
    sigma_even: Expr
    dt: Expr
    kx: Expr
    ky: Expr
    C0: Matrix
    C_ext: Matrix
    C_fb: Matrix | None


@dataclass(frozen=True)
class D2Q9EquivalentCoefficients:
    """Generated second- and fourth-order data under the fixed convention."""

    diffusion: Expr
    gamma_axis3: Expr | None = None
    gamma_diagonal3: Expr | None = None
    gamma_qq3: Expr | None = None
    gamma_axis4: Expr | None = None
    gamma_diagonal4: Expr | None = None
    gamma_qq4: Expr | None = None
    gamma_mixed4: Expr | None = None
    pde_c40: Expr | None = None
    pde_c22: Expr | None = None
    pde_k_axis: Expr | None = None
    pde_k_equal_diagonal: Expr | None = None
    isotropy_residual: Expr | None = None
    cancellation_residual: tuple[Expr, Expr] | None = None


@dataclass(frozen=True)
class QuarticConditionSystem:
    """A generated complete-cancellation system and its exact solutions."""

    c40: Expr
    c22: Expr
    isotropy_residual: Expr
    equations: tuple[Eq, Eq]
    solve_for: tuple[Expr, ...]
    solutions: tuple[dict[Expr, Expr], ...]
    status: str


@dataclass(frozen=True)
class PrintedDuboisCoefficients:
    """The D2Q9 expressions printed after Dubois--Lallemand Eq. (44)."""

    kappa40: Expr
    kappa22: Expr


def build_d2q9_temperature_model(
    *,
    pi: Expr,
    chi_kappa: Expr,
    sigma_odd: Expr,
    sigma_even: Expr,
    dt: Expr = 1,
    kx: Expr = 0,
    ky: Expr = 0,
    require_feedback: bool = True,
) -> D2Q9TemperatureModel:
    """Build the raw matrices in the frozen linear model.

    ``require_feedback=False`` keeps the baseline/external operators available
    when only the local-feedback closure denominator is singular.  In that
    case ``C_fb`` is ``None``.
    """

    lattice = d2q9()
    pi = sympify(pi)
    chi_kappa = sympify(chi_kappa)
    sigma_odd = sympify(sigma_odd)
    sigma_even = sympify(sigma_even)
    dt = sympify(dt)
    kx = sympify(kx)
    ky = sympify(ky)
    cs2 = lattice.cs2

    a = simplify(cs2 + pi)
    b = simplify((1 - chi_kappa) * cs2)
    d = simplify(a - b)
    if a == 0:
        raise ValueError("a=0 is a singular equilibrium closure")
    feedback_denominator = simplify(a + 2 * b * sigma_odd)
    if require_feedback and feedback_denominator == 0:
        raise ValueError("a+2*b*sigma_odd=0 is a singular feedback closure")
    s_odd = simplify(1 / (sigma_odd + Rational(1, 2)))
    s_even = simplify(1 / (sigma_even + Rational(1, 2)))

    ell = Matrix([[1] * 9])
    e = Matrix(
        [
            simplify(weight + pi * correction / cs2)
            for weight, correction in zip(
                lattice.weights, lattice.lambda_t, strict=True
            )
        ]
    )
    G = e * ell
    p_plus, p_minus = parity_projectors(lattice)
    S = s_even * p_plus + s_odd * p_minus
    H = Matrix(
        [
            [
                simplify(weight * velocity[axis] * d / cs2)
                for axis in range(2)
            ]
            for weight, velocity in zip(
                lattice.weights, lattice.velocities, strict=True
            )
        ]
    )
    J = Matrix(
        [
            [velocity[axis] for velocity in lattice.velocities]
            for axis in range(2)
        ]
    )

    identity = eye(9)
    nonequilibrium = identity - G
    C0 = identity - S * nonequilibrium
    wavevector = Matrix([kx, ky])
    C_ext = (
        C0
        + I * dt * (1 - s_odd / 2) * H * wavevector * ell
    )
    C_fb = None
    if feedback_denominator != 0:
        C_fb = (
            C0
            - 2
            * (1 - s_odd / 2)
            / feedback_denominator
            * H
            * J
            * nonequilibrium
        )

    return D2Q9TemperatureModel(
        lattice=lattice,
        ell=ell,
        e=e,
        G=G,
        p_plus=p_plus,
        p_minus=p_minus,
        S=S,
        H=H,
        J=J,
        a=a,
        b=b,
        d=d,
        s_odd=s_odd,
        s_even=s_even,
        sigma_odd=sigma_odd,
        sigma_even=sigma_even,
        dt=dt,
        kx=kx,
        ky=ky,
        C0=C0.applyfunc(simplify),
        C_ext=C_ext.applyfunc(simplify),
        C_fb=None if C_fb is None else C_fb.applyfunc(simplify),
    )


def d2q9_amplification(
    *,
    case: str,
    pi: Expr,
    chi_kappa: Expr,
    sigma_flux: Expr,
    sigma_odd_ghost: Expr,
    sigma_even: Expr,
    dt: Expr = 1,
    kx: Expr = 0,
    ky: Expr = 0,
) -> Matrix:
    """Return ``E(k) C(k)`` under the fixed negative streaming phase."""

    if case not in {"baseline", "external", "feedback"}:
        raise ValueError("case must be 'baseline', 'external', or 'feedback'")
    pi = sympify(pi)
    chi_kappa = sympify(chi_kappa)
    sigma_flux = sympify(sigma_flux)
    sigma_odd_ghost = sympify(sigma_odd_ghost)
    sigma_even = sympify(sigma_even)
    dt = sympify(dt)
    kx = sympify(kx)
    ky = sympify(ky)
    if case == "baseline" and (pi != 0 or chi_kappa != 0):
        raise ValueError("true baseline requires pi=chi_kappa=0")

    model = build_d2q9_temperature_model(
        pi=pi,
        chi_kappa=chi_kappa,
        sigma_odd=sigma_odd_ghost,
        sigma_even=sigma_even,
        dt=dt,
        kx=kx,
        ky=ky,
        require_feedback=case == "feedback",
    )
    expected_sigma_flux = (
        sigma_odd_ghost
        if case in {"baseline", "external"}
        else simplify(model.b * sigma_odd_ghost / model.a)
    )
    if simplify(sigma_flux - expected_sigma_flux) != 0:
        raise ValueError(
            "sigma_flux does not match the selected actual scheme"
        )
    if case == "feedback":
        assert model.C_fb is not None
        collision = model.C_fb
    else:
        collision = {"baseline": model.C0, "external": model.C_ext}[case]
    return _streaming_matrix(model) * collision


def three_block_amplification(
    *,
    block_case: str,
    pi: Expr,
    chi_kappa: Expr,
    sigma_flux: Expr,
    sigma_odd_ghost: Expr,
    sigma_even: Expr,
    dt: Expr = 1,
    kx: Expr = 0,
    ky: Expr = 0,
) -> Matrix:
    """Return the explicit three-shift block generator before specialization."""

    if block_case not in {"homogeneous", "external"}:
        raise ValueError("block_case must be 'homogeneous' or 'external'")
    model = build_d2q9_temperature_model(
        pi=pi,
        chi_kappa=chi_kappa,
        sigma_odd=sigma_odd_ghost,
        sigma_even=sigma_even,
        dt=dt,
        kx=kx,
        ky=ky,
        require_feedback=False,
    )
    sigma_flux = sympify(sigma_flux)
    p_flux, p_odd_ghost = _flux_projectors(model)
    s_flux = simplify(1 / (sigma_flux + Rational(1, 2)))
    s_odd_ghost = simplify(
        1 / (sympify(sigma_odd_ghost) + Rational(1, 2))
    )
    s_even = simplify(1 / (sympify(sigma_even) + Rational(1, 2)))
    relaxation = (
        s_flux * p_flux
        + s_odd_ghost * p_odd_ghost
        + s_even * model.p_plus
    )
    collision = eye(9) - relaxation * (eye(9) - model.G)
    if block_case == "external":
        collision += (
            I
            * sympify(dt)
            * (1 - s_flux / 2)
            * model.H
            * Matrix([sympify(kx), sympify(ky)])
            * model.ell
        )
    return _streaming_matrix(model) * collision


def amplification_route(
    *,
    case: str,
    pi: Expr,
    chi_kappa: Expr,
    sigma_flux: Expr,
    sigma_odd_ghost: Expr,
    sigma_even: Expr,
    dt: Expr = 1,
    order: int = 2,
) -> D2Q9EquivalentCoefficients:
    """Generate the hydrodynamic root by a formal amplification series."""

    if order not in {2, 4}:
        raise ValueError("the amplification route supports order=2 or order=4")
    context = _block_context(
        case=case,
        pi=pi,
        chi_kappa=chi_kappa,
        sigma_flux=sigma_flux,
        sigma_odd_ghost=sigma_odd_ghost,
        sigma_even=sigma_even,
    )
    axis_terms = _amplification_direction(context, (1, 0), order)
    diffusion = simplify(axis_terms[2] * sympify(dt))
    if order == 2:
        return D2Q9EquivalentCoefficients(diffusion=diffusion)
    diagonal_terms = _amplification_direction(context, (1, 1), order)
    return _coefficients_from_spectral_data(
        diffusion=diffusion,
        axis3=axis_terms[3],
        diagonal3=diagonal_terms[3],
        axis4=axis_terms[4],
        diagonal4=diagonal_terms[4],
        dt=sympify(dt),
    )


def taylor_moment_route(
    *,
    case: str,
    pi: Expr,
    chi_kappa: Expr,
    sigma_flux: Expr,
    sigma_odd_ghost: Expr,
    sigma_even: Expr,
    dt: Expr = 1,
    order: int = 2,
) -> D2Q9EquivalentCoefficients:
    """Generate the PDE by physical-space Taylor/block elimination."""

    if order not in {2, 4}:
        raise ValueError("the Taylor/moment route supports order=2 or order=4")
    context = _block_context(
        case=case,
        pi=pi,
        chi_kappa=chi_kappa,
        sigma_flux=sigma_flux,
        sigma_odd_ghost=sigma_odd_ghost,
        sigma_even=sigma_even,
    )
    axis_terms = _directional_taylor_coefficients(context, (1, 0), order)
    diffusion = simplify(axis_terms[1] * sympify(dt))
    if order == 2:
        return D2Q9EquivalentCoefficients(diffusion=diffusion)
    diagonal_terms = _directional_taylor_coefficients(context, (1, 1), order)
    return _coefficients_from_pde_data(
        diffusion=diffusion,
        axis3=axis_terms[2],
        diagonal3=diagonal_terms[2],
        c40=axis_terms[3],
        c22=diagonal_terms[3] - 2 * axis_terms[3],
        dt=sympify(dt),
    )


def printed_dubois_coefficients(
    *,
    sigma1: Expr,
    sigma3: Expr,
    sigma5: Expr,
    sigma7: Expr,
    sigma8: Expr,
    xi: Expr,
    a4: Expr,
) -> PrintedDuboisCoefficients:
    """Encode the two printed D2Q9 quartic expressions without correction."""

    sigma1 = sympify(sigma1)
    sigma3 = sympify(sigma3)
    sigma5 = sympify(sigma5)
    sigma7 = sympify(sigma7)
    sigma8 = sympify(sigma8)
    xi = sympify(xi)
    a4 = sympify(a4)
    kappa40 = sigma1 * (
        2 * sigma5 * (sigma7 - sigma3) * (a4 - 4)
        + 6
        * xi
        * (
            1
            - sigma1 * sigma7
            - 5 * sigma1 * sigma3
            + 2 * sigma5 * (sigma7 - sigma3)
        )
    )
    kappa22 = (
        2
        * (
            sigma1
            + sigma5
            - 2 * sigma1 * sigma5 * (sigma3 + sigma7 + 4 * sigma8)
        )
        * (a4 - 4)
        + 12
        * xi
        * (
            sigma5
            + 3 * sigma1
            - 2 * sigma1 * sigma5 * (sigma3 + sigma7)
            - 2 * sigma1 * sigma3 * sigma5
            - 8 * sigma1 * sigma8 * (sigma1 + sigma5)
            + sigma1**2 * sigma7
        )
    )
    return PrintedDuboisCoefficients(
        kappa40=simplify(kappa40),
        kappa22=simplify(kappa22),
    )


def quartic_condition_system(
    coefficients: D2Q9EquivalentCoefficients,
    *,
    solve_for: tuple[Expr, ...],
) -> QuarticConditionSystem:
    """Generate and solve ``C40=C22=0`` from any route-produced result."""

    if coefficients.pde_c40 is None or coefficients.pde_c22 is None:
        raise ValueError("fourth-order coefficients are required")
    variables = tuple(sympify(variable) for variable in solve_for)
    if not variables:
        raise ValueError("at least one solve variable is required")
    c40 = simplify(coefficients.pde_c40)
    c22 = simplify(coefficients.pde_c22)
    raw_solutions = solve((c40, c22), variables, dict=True)
    solutions = tuple(
        {
            variable: simplify(value)
            for variable, value in solution.items()
        }
        for solution in raw_solutions
    )
    if c40 == 0 and c22 == 0:
        status = "identically_satisfied"
    elif solutions:
        status = "solved"
    else:
        status = "incompatible"
    return QuarticConditionSystem(
        c40=c40,
        c22=c22,
        isotropy_residual=simplify(c22 - 2 * c40),
        equations=(Eq(c40, 0), Eq(c22, 0)),
        solve_for=variables,
        solutions=solutions,
        status=status,
    )


def _flux_projectors(model: D2Q9TemperatureModel) -> tuple[Matrix, Matrix]:
    """Return physical-flux and odd-ghost projectors without dividing by d."""

    cs2 = model.lattice.cs2
    flux_injection = Matrix(
        [
            [weight * velocity[axis] / cs2 for axis in range(2)]
            for weight, velocity in zip(
                model.lattice.weights, model.lattice.velocities, strict=True
            )
        ]
    )
    p_flux = flux_injection * model.J
    return p_flux, model.p_minus - p_flux


def _streaming_matrix(model: D2Q9TemperatureModel) -> Matrix:
    return diag(
        *(
            exp(
                -I
                * model.dt
                * (model.kx * velocity_x + model.ky * velocity_y)
            )
            for velocity_x, velocity_y in model.lattice.velocities
        )
    )


@dataclass(frozen=True)
class _BlockContext:
    model: D2Q9TemperatureModel
    collision: Matrix
    relaxation_inverse: Matrix
    source_gradient: Matrix


def _block_context(
    *,
    case: str,
    pi: Expr,
    chi_kappa: Expr,
    sigma_flux: Expr,
    sigma_odd_ghost: Expr,
    sigma_even: Expr,
) -> _BlockContext:
    """Specialize the independent three-block algebra to an actual scheme."""

    if case not in {"baseline", "external", "feedback"}:
        raise ValueError("case must be 'baseline', 'external', or 'feedback'")
    pi = sympify(pi)
    chi_kappa = sympify(chi_kappa)
    sigma_flux = sympify(sigma_flux)
    sigma_odd_ghost = sympify(sigma_odd_ghost)
    sigma_even = sympify(sigma_even)
    if case == "baseline" and (pi != 0 or chi_kappa != 0):
        raise ValueError("true baseline requires pi=chi_kappa=0")
    model = build_d2q9_temperature_model(
        pi=pi,
        chi_kappa=chi_kappa,
        sigma_odd=sigma_odd_ghost,
        sigma_even=sigma_even,
        require_feedback=case == "feedback",
    )
    expected_sigma_flux = (
        sigma_odd_ghost
        if case in {"baseline", "external"}
        else simplify(model.b * sigma_odd_ghost / model.a)
    )
    if simplify(sigma_flux - expected_sigma_flux) != 0:
        raise ValueError(
            "sigma_flux does not match the selected actual scheme"
        )

    p_flux, p_odd_ghost = _flux_projectors(model)
    s_flux = simplify(1 / (sigma_flux + Rational(1, 2)))
    s_odd_ghost = simplify(1 / (sigma_odd_ghost + Rational(1, 2)))
    s_even = simplify(1 / (sigma_even + Rational(1, 2)))
    relaxation = (
        s_flux * p_flux
        + s_odd_ghost * p_odd_ghost
        + s_even * model.p_plus
    )
    relaxation_inverse = (
        (sigma_flux + Rational(1, 2)) * p_flux
        + (sigma_odd_ghost + Rational(1, 2)) * p_odd_ghost
        + (sigma_even + Rational(1, 2)) * model.p_plus
    )
    collision = eye(9) - relaxation * (eye(9) - model.G)
    source_gradient = (
        (1 - s_flux / 2) * model.H
        if case == "external"
        else zeros(9, 2)
    )
    return _BlockContext(
        model=model,
        collision=collision,
        relaxation_inverse=relaxation_inverse,
        source_gradient=source_gradient,
    )


def _amplification_direction(
    context: _BlockContext,
    direction: tuple[int, int],
    order: int,
) -> list[Expr]:
    """Return coefficients of minus log of the hydrodynamic root."""

    model = context.model
    direction_vector = Matrix(direction)
    collision_terms = [context.collision]
    collision_terms.append(
        I * context.source_gradient * direction_vector * model.ell
    )
    collision_terms.extend(zeros(9) for _ in range(order - 1))
    phase_directions = tuple(
        velocity_x * direction[0] + velocity_y * direction[1]
        for velocity_x, velocity_y in model.lattice.velocities
    )
    streaming_terms = [
        diag(
            *(
                (-I * phase) ** degree / factorial(degree)
                for phase in phase_directions
            )
        )
        for degree in range(order + 1)
    ]
    operator_terms = [
        sum(
            (
                streaming_terms[index] * collision_terms[degree - index]
                for index in range(degree + 1)
            ),
            zeros(9),
        )
        for degree in range(order + 1)
    ]

    vector_terms = [model.e]
    value_terms: list[Expr] = [sympify(1)]
    nonequilibrium = eye(9) - model.G
    for degree in range(1, order + 1):
        right_hand_side = -sum(
            (
                operator_terms[index] * vector_terms[degree - index]
                for index in range(1, degree + 1)
            ),
            zeros(9, 1),
        )
        right_hand_side += sum(
            (
                value_terms[index] * vector_terms[degree - index]
                for index in range(1, degree)
            ),
            zeros(9, 1),
        )
        value = simplify(-(model.ell * right_hand_side)[0])
        vector = (
            -context.relaxation_inverse
            * nonequilibrium
            * right_hand_side
        ).applyfunc(simplify)
        value_terms.append(value)
        vector_terms.append(vector)

    decay_terms = [sympify(0)] * (order + 1)
    if order >= 1:
        decay_terms[1] = simplify(-value_terms[1])
    if order >= 2:
        decay_terms[2] = simplify(
            -(value_terms[2] - value_terms[1] ** 2 / 2)
        )
    if order >= 3:
        decay_terms[3] = simplify(
            -(
                value_terms[3]
                - value_terms[1] * value_terms[2]
                + value_terms[1] ** 3 / 3
            )
        )
    if order >= 4:
        decay_terms[4] = simplify(
            -(
                value_terms[4]
                - value_terms[1] * value_terms[3]
                - value_terms[2] ** 2 / 2
                + value_terms[1] ** 2 * value_terms[2]
                - value_terms[1] ** 4 / 4
            )
        )
    return decay_terms


def _directional_taylor_coefficients(
    context: _BlockContext,
    direction: tuple[int, int],
    order: int,
) -> list[Expr]:
    epsilon = symbols("epsilon")
    model = context.model
    population_terms = [model.e]
    time_terms: list[Expr] = []
    time_operator = sympify(0)
    nonequilibrium = eye(9) - model.G

    for degree in range(1, order + 1):
        residual = _directional_taylor_residual(
            context,
            population_terms,
            time_operator,
            epsilon,
            direction,
            degree,
        )
        known = residual.applyfunc(
            lambda value: expand(value).coeff(epsilon, degree)
        )
        time_term = simplify(-(model.ell * known)[0])
        population_term = (
            -context.relaxation_inverse * nonequilibrium * known
        ).applyfunc(simplify)
        population_terms.append(population_term)
        time_terms.append(time_term)
        time_operator += epsilon ** (degree - 1) * time_term
    return time_terms


def _directional_taylor_residual(
    context: _BlockContext,
    population_terms: list[Matrix],
    time_operator: Expr,
    epsilon: Expr,
    direction: tuple[int, int],
    degree: int,
) -> Matrix:
    model = context.model
    populations = sum(
        (
            epsilon**index * population
            for index, population in enumerate(population_terms)
        ),
        zeros(9, 1),
    )
    source_matrix = (
        context.source_gradient * Matrix(direction) * model.ell
    )
    postcollision = (
        context.collision * populations
        + epsilon * source_matrix * populations
    )
    temporal_taylor = sum(
        (
            (epsilon * time_operator) ** power / factorial(power)
            for power in range(degree + 1)
        ),
        sympify(0),
    )
    spatial_directions = tuple(
        velocity_x * direction[0] + velocity_y * direction[1]
        for velocity_x, velocity_y in model.lattice.velocities
    )
    residual = Matrix(
        [
            populations[index] * temporal_taylor
            - postcollision[index]
            * sum(
                (
                    (-epsilon * spatial_directions[index]) ** power
                    / factorial(power)
                    for power in range(degree + 1)
                ),
                sympify(0),
            )
            for index in range(9)
        ]
    )
    return residual.applyfunc(
        lambda value: truncate_in_scale(value, epsilon, degree)
    )


def _coefficients_from_spectral_data(
    *,
    diffusion: Expr,
    axis3: Expr,
    diagonal3: Expr,
    axis4: Expr,
    diagonal4: Expr,
    dt: Expr,
) -> D2Q9EquivalentCoefficients:
    gamma_axis4 = simplify(axis4 * dt**3)
    gamma_diagonal4 = simplify(diagonal4 * dt**3)
    gamma_mixed4 = simplify(gamma_diagonal4 - 2 * gamma_axis4)
    c40 = simplify(-axis4)
    c22 = simplify(-(diagonal4 - 2 * axis4))
    return D2Q9EquivalentCoefficients(
        diffusion=simplify(diffusion),
        gamma_axis3=simplify(axis3 * dt**2),
        gamma_diagonal3=simplify(diagonal3 * dt**2),
        gamma_qq3=simplify(diagonal3 * dt**2),
        gamma_axis4=gamma_axis4,
        gamma_diagonal4=gamma_diagonal4,
        gamma_qq4=gamma_diagonal4,
        gamma_mixed4=gamma_mixed4,
        pde_c40=c40,
        pde_c22=c22,
        pde_k_axis=c40,
        pde_k_equal_diagonal=simplify(c40 / 2 + c22 / 4),
        isotropy_residual=simplify(c22 - 2 * c40),
        cancellation_residual=(c40, c22),
    )


def _coefficients_from_pde_data(
    *,
    diffusion: Expr,
    axis3: Expr,
    diagonal3: Expr,
    c40: Expr,
    c22: Expr,
    dt: Expr,
) -> D2Q9EquivalentCoefficients:
    c40 = simplify(c40)
    c22 = simplify(c22)
    gamma_axis4 = simplify(-dt**3 * c40)
    gamma_mixed4 = simplify(-dt**3 * c22)
    return D2Q9EquivalentCoefficients(
        diffusion=simplify(diffusion),
        gamma_axis3=simplify(I * axis3 * dt**2),
        gamma_diagonal3=simplify(I * diagonal3 * dt**2),
        gamma_qq3=simplify(I * diagonal3 * dt**2),
        gamma_axis4=gamma_axis4,
        gamma_diagonal4=simplify(2 * gamma_axis4 + gamma_mixed4),
        gamma_qq4=simplify(2 * gamma_axis4 + gamma_mixed4),
        gamma_mixed4=gamma_mixed4,
        pde_c40=c40,
        pde_c22=c22,
        pde_k_axis=c40,
        pde_k_equal_diagonal=simplify(c40 / 2 + c22 / 4),
        isotropy_residual=simplify(c22 - 2 * c40),
        cancellation_residual=(c40, c22),
    )
