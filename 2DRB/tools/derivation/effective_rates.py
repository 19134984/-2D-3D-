"""Exact local-feedback elimination and second-order CE residual algebra.

The functions in this module deliberately distinguish nominal TRT rates from
the effective rates of the physical shear, trace, and scalar-flux modes.  They
operate on exact SymPy expressions and do not modify a solver implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from sympy import Expr, Integer, Rational, diff, simplify, symbols, sympify


@dataclass(frozen=True)
class EffectiveModeRate:
    """A nominal Hénon shift and its feedback-modified physical-mode data."""

    mode: str
    sigma_nominal: Expr
    sigma_effective: Expr
    rate_effective: Expr
    transport_coefficient: Expr
    pressure_ratio: Expr | None = None
    cs2: Expr | None = None
    dt: Expr | None = None


@dataclass(frozen=True)
class NominalGhostRate:
    """An unchanged ghost block on the nominal TRT skeleton."""

    mode: str
    parity: str
    sigma_nominal: Expr
    rate_nominal: Expr


@dataclass(frozen=True)
class EquilibriumMomentConstraints:
    """Named coefficients from LBM-CDE equilibrium moments, Eqs. (6) and (17)."""

    flow_second_convection: Expr = Integer(1)
    scalar_first_advection: Expr = Integer(1)
    scalar_second_pressure: Expr = Integer(1)


@dataclass(frozen=True)
class SourceMomentConstraints:
    """Named coefficients from the inverse-designed source moments."""

    flow_first_force: Expr = Integer(1)
    flow_second_u_force: Expr = Integer(1)
    scalar_zeroth_heat: Expr = Integer(1)
    scalar_first_p_grad_t: Expr = Integer(1)
    scalar_first_t_force: Expr = Integer(1)
    scalar_first_q_velocity: Expr = Integer(1)


@dataclass(frozen=True)
class TrapezoidalFactors:
    """Parity-specific explicit source multipliers in the transformed LBE."""

    flow_even: Expr
    flow_odd: Expr
    scalar_even: Expr
    scalar_odd: Expr

    @classmethod
    def from_rates(
        cls,
        *,
        s_f_plus: Expr,
        s_f_minus: Expr,
        s_g_plus: Expr,
        s_g_minus: Expr,
    ) -> "TrapezoidalFactors":
        """Return the exact ``1-s/2`` factors for all four parity blocks."""

        return cls(
            flow_even=simplify(1 - sympify(s_f_plus) / 2),
            flow_odd=simplify(1 - sympify(s_f_minus) / 2),
            scalar_even=simplify(1 - sympify(s_g_plus) / 2),
            scalar_odd=simplify(1 - sympify(s_g_minus) / 2),
        )


def _rate_from_shift(sigma: Expr) -> Expr:
    return simplify(1 / (sympify(sigma) + Rational(1, 2)))


def scalar_flux_effective_rate(
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    *,
    pressure_ratio: Expr = 0,
    cs2: Expr = Rational(1, 3),
    dt: Expr = 1,
) -> EffectiveModeRate:
    """Eliminate Eq. (35) in the odd scalar-flux collision moment.

    ``pressure_ratio`` is the locally frozen value ``p/rho0``.  The returned
    effective rate is a modal identity; it is not a replacement nominal input
    for the TRT collision skeleton.
    """

    sigma_g_minus = sympify(sigma_g_minus)
    chi_kappa = sympify(chi_kappa)
    pressure_ratio = sympify(pressure_ratio)
    cs2 = sympify(cs2)
    dt = sympify(dt)
    sigma_effective = simplify(
        (1 - chi_kappa)
        * cs2
        * sigma_g_minus
        / (cs2 + pressure_ratio)
    )
    transport = simplify(
        (cs2 + pressure_ratio) * dt * sigma_effective
    )
    return EffectiveModeRate(
        mode="scalar_flux",
        sigma_nominal=sigma_g_minus,
        sigma_effective=sigma_effective,
        rate_effective=_rate_from_shift(sigma_effective),
        transport_coefficient=transport,
        pressure_ratio=pressure_ratio,
        cs2=cs2,
        dt=dt,
    )


def shear_effective_rate(
    sigma_f_plus: Expr,
    chi_s: Expr,
    *,
    cs2: Expr = Rational(1, 3),
    dt: Expr = 1,
    mode: str = "off_diagonal",
) -> EffectiveModeRate:
    """Return the feedback-modified shear or deviatoric-diagonal mode rate."""

    if mode not in {"off_diagonal", "deviatoric_diagonal"}:
        raise ValueError(
            "shear mode must be 'off_diagonal' or 'deviatoric_diagonal'"
        )
    sigma_f_plus = sympify(sigma_f_plus)
    chi_s = sympify(chi_s)
    sigma_effective = simplify((1 - chi_s) * sigma_f_plus)
    transport = simplify(sympify(cs2) * sympify(dt) * sigma_effective)
    return EffectiveModeRate(
        mode=mode,
        sigma_nominal=sigma_f_plus,
        sigma_effective=sigma_effective,
        rate_effective=_rate_from_shift(sigma_effective),
        transport_coefficient=transport,
        cs2=sympify(cs2),
        dt=sympify(dt),
    )


def bulk_effective_rate(
    sigma_f_plus: Expr,
    chi_b: Expr,
    *,
    dimension: Expr = 2,
    cs2: Expr = Rational(1, 3),
    dt: Expr = 1,
) -> EffectiveModeRate:
    """Return the trace-block rate and the paper's D-dimensional bulk value."""

    sigma_f_plus = sympify(sigma_f_plus)
    chi_b = sympify(chi_b)
    dimension = sympify(dimension)
    sigma_effective = simplify((1 - chi_b) * sigma_f_plus)
    transport = simplify(
        2
        * sympify(cs2)
        * sympify(dt)
        * sigma_effective
        / dimension
    )
    return EffectiveModeRate(
        mode="trace_bulk",
        sigma_nominal=sigma_f_plus,
        sigma_effective=sigma_effective,
        rate_effective=_rate_from_shift(sigma_effective),
        transport_coefficient=transport,
        cs2=sympify(cs2),
        dt=sympify(dt),
    )


def effective_operator_blocks(
    *,
    sigma_f_plus: Expr,
    sigma_f_minus: Expr,
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_s: Expr,
    chi_b: Expr,
    chi_kappa: Expr,
    pressure_ratio: Expr = 0,
    dimension: Expr = 2,
    cs2: Expr = Rational(1, 3),
    dt: Expr = 1,
) -> Mapping[str, EffectiveModeRate | NominalGhostRate]:
    """Expose the restricted block-MRT operator induced inside TRT parity.

    The input collision still has four nominal TRT shifts.  Feedback changes
    only the physical shear, trace, and scalar-flux blocks; ghosts retain the
    nominal shifts of their containing parity subspace.
    """

    sigma_f_plus = sympify(sigma_f_plus)
    sigma_f_minus = sympify(sigma_f_minus)
    sigma_g_plus = sympify(sigma_g_plus)
    sigma_g_minus = sympify(sigma_g_minus)
    return {
        "flow_shear": shear_effective_rate(
            sigma_f_plus,
            chi_s,
            cs2=cs2,
            dt=dt,
            mode="off_diagonal",
        ),
        "flow_deviatoric": shear_effective_rate(
            sigma_f_plus,
            chi_s,
            cs2=cs2,
            dt=dt,
            mode="deviatoric_diagonal",
        ),
        "flow_bulk": bulk_effective_rate(
            sigma_f_plus,
            chi_b,
            dimension=dimension,
            cs2=cs2,
            dt=dt,
        ),
        "flow_even_ghost": _ghost("flow_even_ghost", "even", sigma_f_plus),
        "flow_odd_ghost": _ghost("flow_odd_ghost", "odd", sigma_f_minus),
        "scalar_flux": scalar_flux_effective_rate(
            sigma_g_minus,
            chi_kappa,
            pressure_ratio=pressure_ratio,
            cs2=cs2,
            dt=dt,
        ),
        "scalar_even_ghost": _ghost(
            "scalar_even_ghost", "even", sigma_g_plus
        ),
        "scalar_odd_ghost": _ghost(
            "scalar_odd_ghost", "odd", sigma_g_minus
        ),
    }


def scalar_variable_pressure_residual(
    frozen_result: EffectiveModeRate,
    *,
    pressure_perturbation: Expr,
    grad_pressure_ratio_dot_grad_temperature: Expr,
    laplacian_temperature: Expr,
) -> Expr:
    """Return the first product-derivative residual beyond a frozen-p mode.

    The rate is frozen at the reference pressure while ``p/rho0`` is perturbed.
    This intentionally does not claim that the local modal identity is a
    variable-coefficient modified equation.
    """

    if frozen_result.mode != "scalar_flux":
        raise ValueError("a frozen scalar-flux result is required")
    delta_pi = sympify(pressure_perturbation)
    grad_product = sympify(grad_pressure_ratio_dot_grad_temperature)
    laplacian = sympify(laplacian_temperature)
    return simplify(
        frozen_result.sigma_effective
        * _transport_dt(frozen_result)
        * (grad_product + delta_pi * laplacian)
    )


def second_order_residual_table(
    *,
    s_f_plus: Expr | None = None,
    s_f_minus: Expr | None = None,
    s_g_plus: Expr | None = None,
    s_g_minus: Expr | None = None,
    chi_s: Expr = 0,
    chi_b: Expr = 0,
    chi_kappa: Expr = 0,
    cs2: Expr = Rational(1, 3),
    dt: Expr = 1,
    dimension: Expr = 2,
    pressure_ratio: Expr = 0,
    equilibrium: EquilibriumMomentConstraints | None = None,
    sources: SourceMomentConstraints | None = None,
    trapezoidal: TrapezoidalFactors | None = None,
) -> Mapping[str, Expr | int | str]:
    """Generate second-order residuals from named moments and CE transfers.

    In the transformed-LBE CE hierarchy, a nonconserved source block with
    explicit factor ``b`` enters a first nonequilibrium moment through ``b/s``.
    Division by ``sigma=1/s-1/2`` gives its normalized constitutive transfer.
    A conserved force/heat source also receives ``s/2`` from relaxing the
    half-source nonequilibrium moment, so its Euler balance coefficient is
    ``b+s/2``.  The correct ``b=1-s/2`` makes both routes one; perturbing a
    moment or factor makes the corresponding generated residual nonzero.
    """

    default_rates = symbols(
        "s_f_plus s_f_minus s_g_plus s_g_minus", nonzero=True
    )
    rate_values = (s_f_plus, s_f_minus, s_g_plus, s_g_minus)
    s_f_plus, s_f_minus, s_g_plus, s_g_minus = tuple(
        default if value is None else sympify(value)
        for value, default in zip(rate_values, default_rates, strict=True)
    )
    equilibrium = equilibrium or EquilibriumMomentConstraints()
    sources = sources or SourceMomentConstraints()
    trapezoidal = trapezoidal or TrapezoidalFactors.from_rates(
        s_f_plus=s_f_plus,
        s_f_minus=s_f_minus,
        s_g_plus=s_g_plus,
        s_g_minus=s_g_minus,
    )

    sigma_f_plus = simplify(1 / s_f_plus - Rational(1, 2))
    sigma_f_minus = simplify(1 / s_f_minus - Rational(1, 2))
    sigma_g_plus = simplify(1 / s_g_plus - Rational(1, 2))
    sigma_g_minus = simplify(1 / s_g_minus - Rational(1, 2))

    flow_stress_transfer = _normalized_source_transfer(
        trapezoidal.flow_even,
        s_f_plus,
        sigma_f_plus,
    )
    flow_force_balance = simplify(
        sympify(trapezoidal.flow_odd) + s_f_minus / 2
    )
    scalar_heat_balance = simplify(
        sympify(trapezoidal.scalar_even) + s_g_plus / 2
    )
    scalar_flux_transfer = _normalized_source_transfer(
        trapezoidal.scalar_odd,
        s_g_minus,
        sigma_g_minus,
    )

    nu_shear = simplify(
        (1 - sympify(chi_s)) * sympify(cs2) * sympify(dt) * sigma_f_plus
    )
    nu_bulk = simplify(
        2
        * (1 - sympify(chi_b))
        * sympify(cs2)
        * sympify(dt)
        * sigma_f_plus
        / sympify(dimension)
    )
    kappa = simplify(
        (1 - sympify(chi_kappa))
        * sympify(cs2)
        * sympify(dt)
        * sigma_g_minus
    )

    residuals: dict[str, Expr | int | str] = {
        "p_grad_T": simplify(
            sympify(equilibrium.scalar_second_pressure)
            - scalar_flux_transfer
            * sympify(sources.scalar_first_p_grad_t)
        ),
        "T_F": simplify(
            sympify(equilibrium.scalar_first_advection)
            * flow_force_balance
            * sympify(sources.flow_first_force)
            - scalar_flux_transfer
            * sympify(sources.scalar_first_t_force)
        ),
        "Q_u": simplify(
            sympify(equilibrium.scalar_first_advection)
            * scalar_heat_balance
            * sympify(sources.scalar_zeroth_heat)
            - scalar_flux_transfer
            * sympify(sources.scalar_first_q_velocity)
        ),
        "u_F": simplify(
            sympify(equilibrium.flow_second_convection)
            * flow_force_balance
            * sympify(sources.flow_first_force)
            - flow_stress_transfer
            * sympify(sources.flow_second_u_force)
        ),
        "d_t_F": simplify(
            sympify(trapezoidal.flow_odd) / s_f_minus - sigma_f_minus
        ),
        "d_t_Q": simplify(
            sympify(trapezoidal.scalar_even) / s_g_plus - sigma_g_plus
        ),
        "s_f_minus_transport": simplify(diff(nu_shear, s_f_minus)),
        "s_g_plus_transport": simplify(diff(kappa, s_g_plus)),
        "nu_shear": nu_shear,
        "nu_bulk": nu_bulk,
        "kappa": kappa,
        "sigma_scalar_flux_frozen": simplify(
            (1 - sympify(chi_kappa))
            * sympify(cs2)
            * sigma_g_minus
            / (sympify(cs2) + sympify(pressure_ratio))
        ),
        "coefficient_variation_epsilon_order": 2,
        "coefficient_assumption": "frozen_through_ce2",
        "constant_coefficient_first_omitted_epsilon_order": 3,
        "first_omitted_epsilon_order": 3,
        "first_omitted_mach_order": 3,
    }
    return residuals


def _normalized_source_transfer(factor: Expr, rate: Expr, sigma: Expr) -> Expr:
    return simplify(sympify(factor) / (sympify(rate) * sympify(sigma)))


def _ghost(mode: str, parity: str, sigma: Expr) -> NominalGhostRate:
    sigma = sympify(sigma)
    return NominalGhostRate(
        mode=mode,
        parity=parity,
        sigma_nominal=sigma,
        rate_nominal=_rate_from_shift(sigma),
    )


def _transport_dt(result: EffectiveModeRate) -> Expr:
    """Recover ``dt`` multiplying a frozen scalar modal flux."""

    if result.dt is None:
        raise ValueError("time step is unavailable")
    return sympify(result.dt)
