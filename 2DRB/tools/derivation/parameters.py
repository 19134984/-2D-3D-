"""Exact compatibility reports assembled from reviewed bulk and wall APIs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from sympy import (
    And,
    Expr,
    Rational,
    Q,
    limit,
    simplify,
    solve,
    sqrt,
    symbols,
    sympify,
    together,
)

from tools.derivation.boundary import (
    BoundaryResidual,
    MagicClassification,
    classify_magic,
    temperature_abb_residual,
    velocity_bb_residual,
)
from tools.derivation.d2q9_temperature import (
    D2Q9EquivalentCoefficients,
    QuarticConditionSystem,
    amplification_route,
    canonical_quartic_condition,
    quartic_condition_system,
)


PARAMETER_STATUSES = {
    "feasible_exact",
    "feasible_restricted",
    "no_feasible_solution",
    "degenerate_branch",
    "boundary_correction_required",
    "mrt_extension_required",
}


@dataclass(frozen=True)
class ConsumedEvidence:
    """Actual reviewed objects used to construct one compatibility report."""

    bulk_route: D2Q9EquivalentCoefficients | None
    bulk_conditions: QuarticConditionSystem | None
    boundary_residual: BoundaryResidual | None = None
    boundary_classification: MagicClassification | None = None
    bulk_specialization: Mapping[str, Expr] | None = None
    bulk_provenance: str | None = None


@dataclass(frozen=True)
class ParameterReport:
    """Structured exact solver result; never only a tuple of rates."""

    status: str
    case: str
    exact_substitutions: Mapping[str, Expr]
    collision_rates: Mapping[str, Expr]
    open_interval_checks: Mapping[str, object]
    violated_constraints: tuple[str, ...]
    assumptions: tuple[str, ...]
    minimal_extension: str | None
    consumed_evidence: ConsumedEvidence
    is_conditional: bool = False
    feasibility_conditions: Mapping[str, object] | None = None


def recover_baseline_quartic_family(*, sigma_odd: Expr) -> ParameterReport:
    """Generate the baseline family from Task 5's fourth-order API."""

    sigma_odd = sympify(sigma_odd)
    sigma_even = symbols("sigma_e")
    route = amplification_route(
        case="baseline",
        pi=0,
        chi_kappa=0,
        sigma_flux=sigma_odd,
        sigma_odd_ghost=sigma_odd,
        sigma_even=sigma_even,
        order=4,
    )
    conditions = quartic_condition_system(route, solve_for=(sigma_even,))
    if conditions.status != "solved" or len(conditions.solutions) != 1:
        raise ValueError("Task 5 did not return one baseline quartic family")
    solved_even = simplify(conditions.solutions[0][sigma_even])
    nominal_odd_rate = simplify(1 / (sigma_odd + sympify(1) / 2))
    nominal_even_rate = simplify(1 / (solved_even + sympify(1) / 2))
    physical_flux_rate = nominal_odd_rate
    return ParameterReport(
        status="feasible_exact",
        case="baseline",
        exact_substitutions={
            "sigma_odd": sigma_odd,
            "sigma_even": solved_even,
            "sigma_flux": sigma_odd,
        },
        collision_rates={
            "nominal_odd": nominal_odd_rate,
            "nominal_even": nominal_even_rate,
            "physical_flux": physical_flux_rate,
        },
        open_interval_checks={
            "nominal_odd_rate": _open_rate_check(nominal_odd_rate),
            "nominal_even_rate": _open_rate_check(nominal_even_rate),
            "physical_flux_rate": _open_rate_check(physical_flux_rate),
        },
        violated_constraints=(),
        assumptions=(
            "task5_frozen_baseline_d2q9",
            "complete_bulk_quartic_cancellation",
            "positive_henon_branch_when_numerically_specialized",
        ),
        minimal_extension=None,
        consumed_evidence=ConsumedEvidence(
            bulk_route=route,
            bulk_conditions=conditions,
        ),
    )


def derive_scalar_compatibility(
    *,
    case: str,
    pressure_ratio: Expr,
    chi_kappa: Expr,
    cs2: Expr = Rational(1, 3),
    dt: Expr = 1,
) -> ParameterReport:
    """Eliminate Task 5 bulk and Task 6 restricted ABB constraints exactly."""

    if case not in {"external", "feedback"}:
        raise ValueError("scalar compatibility case must be external or feedback")
    pressure_ratio, chi_kappa, cs2, dt = map(
        sympify, (pressure_ratio, chi_kappa, cs2, dt)
    )
    if cs2 != Rational(1, 3):
        raise ValueError("Task 5's reviewed D2Q9 route requires cs2=1/3")
    a = simplify(cs2 + pressure_ratio)
    b = simplify((1 - chi_kappa) * cs2)
    input_conditions, input_violations, _ = _condition_partition(
        {
            "a_positive": _positive_check(a),
            "b_positive": _positive_check(b),
        }
    )
    if input_violations:
        return ParameterReport(
            status="no_feasible_solution",
            case=case,
            exact_substitutions={"a": a, "b": b},
            collision_rates={},
            open_interval_checks={},
            violated_constraints=input_violations,
            assumptions=("positive_physical_scalar_transport_required",),
            minimal_extension="choose_positive_equilibrium_and_transport_blocks",
            consumed_evidence=ConsumedEvidence(None, None),
            feasibility_conditions=input_conditions,
        )

    sigma_odd, sigma_even = symbols("sigma_o sigma_e")
    canonical = canonical_quartic_condition(
        case=case,
        pi=pressure_ratio,
        chi_kappa=chi_kappa,
        sigma_odd_ghost=sigma_odd,
        sigma_even=sigma_even,
    )
    route = canonical.coefficients
    conditions = canonical.conditions

    regime = (
        "steady_1d_quadratic"
        if case == "external"
        else "steady_1d_quadratic_local_feedback"
    )
    boundary = temperature_abb_residual(
        sigma_g_plus=sigma_even,
        sigma_g_minus=sigma_odd,
        chi_kappa=chi_kappa,
        cs2=cs2,
        pressure_ratio=pressure_ratio,
        regime=regime,
        normalization="transformed_cde_chain",
        scalar_gauge="half_source",
        geometry="flat_grid_aligned_halfway",
    )
    boundary_classification = classify_magic(boundary)
    if boundary_classification.status != "restricted_calibration":
        raise ValueError("Task 6 did not return the restricted ABB calibration")
    boundary_product = boundary.calibration_parameters[0]
    boundary_conditions = boundary.solve_zero_conditions()
    if boundary_product not in boundary_conditions:
        raise ValueError("Task 6 did not solve the restricted ABB product")

    transport_ratio = symbols("kappa_over_dt")
    transport_odd = simplify(transport_ratio / b)
    boundary_even_roots = solve(
        boundary_product.subs(sigma_odd, transport_odd)
        - boundary_conditions[boundary_product],
        sigma_even,
        dict=False,
    )
    if len(boundary_even_roots) != 1:
        raise ValueError("restricted ABB did not determine one nominal even shift")
    boundary_even = simplify(boundary_even_roots[0])
    bulk_even = None
    if conditions.status == "solved" and len(conditions.solutions) == 1:
        bulk_even = conditions.solutions[0].get(sigma_even)
    compatibility = together(
        canonical.undivided_polynomial.subs(
            {
                sigma_odd: transport_odd,
                sigma_even: boundary_even,
            }
        )
    )
    numerator = simplify(compatibility.as_numer_denom()[0])
    squared_roots = solve(numerator, transport_ratio**2, dict=False)
    if len(squared_roots) != 1:
        raise ValueError(
            "canonical bulk plus ABB did not determine one squared transport ratio"
        )
    required_squared = simplify(squared_roots[0])
    branch_assumptions = (
        ("a_equals_one_direct_condition", "task5_canonical_undivided_condition")
        if a == 1
        else ("a_not_equal_one_main_branch", "task5_canonical_undivided_condition")
    )

    evidence = ConsumedEvidence(
        bulk_route=route,
        bulk_conditions=conditions,
        boundary_residual=boundary,
        boundary_classification=boundary_classification,
        bulk_provenance=canonical.provenance,
    )
    substitutions: dict[str, Expr] = {
        "a": a,
        "b": b,
        "boundary_sigma_even_family": boundary_even,
        "required_transport_ratio_squared": required_squared,
    }
    if bulk_even is not None:
        substitutions["bulk_sigma_even_family"] = bulk_even
    if required_squared == 0 or required_squared.is_negative is True:
        violation = (
            "zero_compatibility_radicand"
            if required_squared == 0
            else "negative_compatibility_radicand"
        )
        return ParameterReport(
            status="no_feasible_solution",
            case=case,
            exact_substitutions=substitutions,
            collision_rates={},
            open_interval_checks={},
            violated_constraints=(violation,),
            assumptions=(
                "task5_frozen_coefficient_complete_quartic_cancellation",
                "positive_physical_main_branch",
            )
            + branch_assumptions
            + boundary.assumptions,
            minimal_extension="explicit_restricted_abb_boundary_correction",
            consumed_evidence=evidence,
        )

    required_ratio = sqrt(required_squared)
    required_odd = simplify(required_ratio / b)
    required_even = simplify(boundary_even.subs(transport_ratio, required_ratio))
    required_flux = simplify(required_ratio / a)
    rates = _three_rates(required_odd, required_even, required_flux)
    rate_checks = _three_rate_checks(rates)
    physical_conditions, physical_violations, unresolved_conditions = (
        _condition_partition(
            {
                "a_positive": _positive_check(a),
                "b_positive": _positive_check(b),
                "compatibility_radicand_positive": _positive_check(
                    required_squared
                ),
                "sigma_odd_positive": _positive_check(required_odd),
                "sigma_even_positive": _positive_check(required_even),
                "sigma_flux_positive": _positive_check(required_flux),
                **rate_checks,
            }
        )
    )
    substitutions.update(
        {
            "required_transport_ratio": required_ratio,
            "required_kappa": simplify(dt * required_ratio),
            "sigma_odd": required_odd,
            "sigma_even": required_even,
            "sigma_flux": required_flux,
        }
    )
    if physical_violations:
        return ParameterReport(
            status="no_feasible_solution",
            case=case,
            exact_substitutions=substitutions,
            collision_rates=rates,
            open_interval_checks=rate_checks,
            violated_constraints=physical_violations,
            assumptions=(
                "task5_frozen_coefficient_complete_quartic_cancellation",
                "positive_physical_main_branch",
            )
            + branch_assumptions
            + boundary.assumptions,
            minimal_extension="choose_positive_physical_scalar_parameters",
            consumed_evidence=evidence,
            feasibility_conditions=physical_conditions,
        )
    return ParameterReport(
        status="feasible_restricted",
        case=case,
        exact_substitutions=substitutions,
        collision_rates=rates,
        open_interval_checks=rate_checks,
        violated_constraints=(),
        assumptions=(
            "task5_frozen_coefficient_complete_quartic_cancellation",
            "positive_physical_main_branch",
            "a_nonzero",
            "b_positive",
            "transport_ratio_positive",
        )
        + branch_assumptions
        + boundary.assumptions,
        minimal_extension=None,
        consumed_evidence=evidence,
        is_conditional=bool(unresolved_conditions),
        feasibility_conditions=unresolved_conditions or None,
    )


def solve_scalar_parameters(
    *,
    case: str,
    pressure_ratio: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    kappa: Expr,
    dt: Expr,
    require_bulk_quartic: bool,
    boundary_policy: str,
) -> ParameterReport:
    """Check one target against bulk quartic and the requested wall policy."""

    if not require_bulk_quartic:
        raise ValueError("this solver entry point requires bulk quartic cancellation")
    if boundary_policy not in {"restricted_abb", "explicit_correction"}:
        raise ValueError("unknown scalar boundary policy")
    pressure_ratio, chi_kappa, cs2, kappa, dt = map(
        sympify, (pressure_ratio, chi_kappa, cs2, kappa, dt)
    )
    a = simplify(cs2 + pressure_ratio)
    b = simplify((1 - chi_kappa) * cs2)
    target_ratio = simplify(kappa / dt)
    input_conditions, input_violations, unresolved_input = _condition_partition(
        {
            "a_positive": _positive_check(a),
            "b_positive": _positive_check(b),
            "transport_ratio_positive": _positive_check(target_ratio),
        }
    )
    feedback_singular = (
        case == "feedback" and simplify(a + 2 * target_ratio) == 0
    )
    if input_violations:
        violations = list(input_violations)
        if feedback_singular:
            violations.append("feedback_closure_singular")
        return ParameterReport(
            status="no_feasible_solution",
            case=case,
            exact_substitutions={
                "a": a,
                "b": b,
                "target_kappa": kappa,
                "target_transport_ratio": target_ratio,
            },
            collision_rates={},
            open_interval_checks={},
            violated_constraints=tuple(violations),
            assumptions=("positive_physical_scalar_transport_required",),
            minimal_extension="change_nonpositive_scalar_transport_inputs",
            consumed_evidence=ConsumedEvidence(None, None),
            feasibility_conditions=input_conditions,
        )
    if feedback_singular:
        return _degenerate_scalar_report(
            case=case,
            a=a,
            b=b,
            violation="feedback_closure_singular",
            extension="change_feedback_closure_or_target_transport",
            extra_substitutions={
                "target_kappa": kappa,
                "target_transport_ratio": target_ratio,
            },
        )
    derived = derive_scalar_compatibility(
        case=case,
        pressure_ratio=pressure_ratio,
        chi_kappa=chi_kappa,
        cs2=cs2,
        dt=dt,
    )
    if chi_kappa.free_symbols:
        substitutions = dict(derived.exact_substitutions)
        substitutions.update(
            {"target_kappa": kappa, "target_transport_ratio": target_ratio}
        )
        return ParameterReport(
            status="no_feasible_solution",
            case=case,
            exact_substitutions=substitutions,
            collision_rates={},
            open_interval_checks={},
            violated_constraints=("chi_kappa_must_be_specified",),
            assumptions=derived.assumptions,
            minimal_extension="specify_chi_kappa_then_recompute_nominal_odd_rate",
            consumed_evidence=derived.consumed_evidence,
            is_conditional=True,
            feasibility_conditions={
                **(derived.feasibility_conditions or {}),
                **unresolved_input,
            },
        )
    if derived.status == "degenerate_branch":
        return derived
    if boundary_policy == "restricted_abb" and derived.status == "no_feasible_solution":
        return derived

    sigma_odd = simplify(target_ratio / b)
    sigma_flux = simplify(target_ratio / a)
    if boundary_policy == "explicit_correction":
        bulk_family = derived.exact_substitutions.get("bulk_sigma_even_family")
        if bulk_family is None:
            return ParameterReport(
                status="mrt_extension_required",
                case=case,
                exact_substitutions=dict(derived.exact_substitutions),
                collision_rates={},
                open_interval_checks={},
                violated_constraints=("a_one_even_shift_not_fixed",),
                assumptions=derived.assumptions,
                minimal_extension=(
                    "split_even_mrt_candidate_requiring_mode_jacobian_derivation"
                ),
                consumed_evidence=derived.consumed_evidence,
            )
        sigma_even = simplify(bulk_family.subs(symbols("sigma_o"), sigma_odd))
    else:
        boundary = derived.consumed_evidence.boundary_residual
        assert boundary is not None
        boundary_product = boundary.calibration_parameters[0]
        boundary_value = boundary.solve_zero_conditions()[boundary_product]
        sigma_even_roots = solve(
            boundary_product.subs(symbols("sigma_o"), sigma_odd) - boundary_value,
            symbols("sigma_e"),
            dict=False,
        )
        if len(sigma_even_roots) != 1:
            raise ValueError("restricted ABB did not determine one target even shift")
        sigma_even = simplify(sigma_even_roots[0])
    rates = _three_rates(sigma_odd, sigma_even, sigma_flux)
    rate_checks = _three_rate_checks(rates)
    compatible = (
        boundary_policy == "explicit_correction"
        or simplify(
            target_ratio**2
            - derived.exact_substitutions["required_transport_ratio_squared"]
        )
        == 0
    )
    physical_checks, physical_violations, unresolved_physical = _condition_partition(
        {
            "a_positive": _positive_check(a),
            "b_positive": _positive_check(b),
            "transport_ratio_positive": _positive_check(target_ratio),
            "sigma_odd_positive": _positive_check(sigma_odd),
            "sigma_even_positive": _positive_check(sigma_even),
            "sigma_flux_positive": _positive_check(sigma_flux),
        }
    )
    _, rate_violations, unresolved_rates = _condition_partition(rate_checks)
    all_physical = not physical_violations and not unresolved_physical
    all_rates = not rate_violations and not unresolved_rates
    if compatible and all_physical and all_rates:
        if boundary_policy == "explicit_correction":
            status = "boundary_correction_required"
            violated = ("restricted_abb_product",)
            extension = "explicit_restricted_abb_boundary_correction"
        else:
            status = "feasible_restricted"
            violated = ()
            extension = None
    else:
        status = "no_feasible_solution"
        violated_items = []
        if not compatible:
            violated_items.append("bulk_quartic_and_restricted_abb")
        violated_items.extend(physical_violations)
        violated_items.extend(rate_violations)
        violated = tuple(violated_items)
        extension = "explicit_restricted_abb_boundary_correction"
    substitutions = dict(derived.exact_substitutions)
    substitutions.update(
        {
            "target_kappa": kappa,
            "target_transport_ratio": target_ratio,
            "sigma_odd": sigma_odd,
            "sigma_even": sigma_even,
            "sigma_flux": sigma_flux,
        }
    )
    return ParameterReport(
        status=status,
        case=case,
        exact_substitutions=substitutions,
        collision_rates=rates,
        open_interval_checks=rate_checks,
        violated_constraints=violated,
        assumptions=derived.assumptions,
        minimal_extension=extension,
        consumed_evidence=derived.consumed_evidence,
        is_conditional=bool(
            unresolved_physical or unresolved_input or unresolved_rates
        ),
        feasibility_conditions={
            **unresolved_input,
            **unresolved_physical,
            **unresolved_rates,
        }
        or None,
    )


def solve_flow_parameters(
    *,
    nu: Expr,
    cs2: Expr,
    dt: Expr,
    chi_s: Expr,
    chi_b: Expr,
    retain_trace_jets: bool,
) -> ParameterReport:
    """Solve the Task 6 restricted shear candidate or preserve its general rejection."""

    nu, cs2, dt, chi_s, chi_b = map(sympify, (nu, cs2, dt, chi_s, chi_b))
    sigma_plus_symbol, sigma_minus_symbol = symbols(
        "sigma_f_plus sigma_f_minus"
    )
    shear_denominator = simplify(cs2 * dt * (1 - chi_s))
    if shear_denominator == 0:
        return ParameterReport(
            status="no_feasible_solution",
            case="flow",
            exact_substitutions={
                "nu": nu,
                "cs2": cs2,
                "dt": dt,
                "chi_s": chi_s,
                "chi_b": chi_b,
            },
            collision_rates={},
            open_interval_checks={},
            violated_constraints=("zero_shear_transport_denominator",),
            assumptions=("positive_physical_flow_transport_required",),
            minimal_extension="restore_positive_shear_transport_denominator",
            consumed_evidence=ConsumedEvidence(None, None),
        )
    sigma_plus = simplify(nu / shear_denominator)
    physical_shear_shift = simplify((1 - chi_s) * sigma_plus)
    physical_bulk_shift = simplify((1 - chi_b) * sigma_plus)
    nu_bulk_2d = simplify(cs2 * dt * physical_bulk_shift)
    partial_rates = {
        "nominal_even": simplify(1 / (sigma_plus + Rational(1, 2))),
        "physical_shear": simplify(
            1 / (physical_shear_shift + Rational(1, 2))
        ),
        "physical_bulk": simplify(
            1 / (physical_bulk_shift + Rational(1, 2))
        ),
    }
    partial_checks = {
        "nominal_even_rate": _open_rate_check(partial_rates["nominal_even"]),
        "physical_shear_rate": _open_rate_check(partial_rates["physical_shear"]),
        "physical_bulk_rate": _open_rate_check(partial_rates["physical_bulk"]),
    }
    flow_conditions, flow_violations, unresolved_flow = _condition_partition(
        {
            "nu_positive": _positive_check(nu),
            "cs2_positive": _positive_check(cs2),
            "dt_positive": _positive_check(dt),
            "one_minus_chi_s_positive": _positive_check(1 - chi_s),
            "one_minus_chi_b_positive": _positive_check(1 - chi_b),
            "sigma_f_plus_positive": _positive_check(sigma_plus),
            "sigma_shear_physical_positive": _positive_check(
                physical_shear_shift
            ),
            "sigma_bulk_physical_positive": _positive_check(
                physical_bulk_shift
            ),
            "nu_bulk_2d_positive": _positive_check(nu_bulk_2d),
            **partial_checks,
        }
    )
    if flow_violations:
        return ParameterReport(
            status="no_feasible_solution",
            case="flow",
            exact_substitutions={
                "nu": nu,
                "sigma_f_plus": sigma_plus,
                "sigma_shear_physical": physical_shear_shift,
                "sigma_bulk_physical": physical_bulk_shift,
                "nu_bulk_2d": nu_bulk_2d,
            },
            collision_rates=partial_rates,
            open_interval_checks=partial_checks,
            violated_constraints=flow_violations,
            assumptions=("physical_feasibility_precedes_wall_classification",),
            minimal_extension="change_nonpositive_or_out_of_range_flow_parameters",
            consumed_evidence=ConsumedEvidence(None, None),
            feasibility_conditions=flow_conditions,
        )

    if retain_trace_jets:
        residual = velocity_bb_residual(
            sigma_f_plus=sigma_plus_symbol,
            sigma_f_minus=sigma_minus_symbol,
            chi_s=chi_s,
            chi_b=chi_b,
            force_mode="general_source_aware",
            normalization="source_feedback_henon_product",
            momentum_gauge="half_force",
            geometry="flat_grid_aligned_halfway",
        )
        classification = classify_magic(residual)
        evidence = ConsumedEvidence(
            bulk_route=None,
            bulk_conditions=None,
            boundary_residual=residual,
            boundary_classification=classification,
        )
        substitutions = {
            "nu": nu,
            "sigma_f_plus": sigma_plus,
            "sigma_shear_physical": physical_shear_shift,
            "sigma_bulk_physical": physical_bulk_shift,
            "nu_bulk_2d": nu_bulk_2d,
        }
        if classification.rate_compatibility_status == "no_single_magic":
            return ParameterReport(
                status="mrt_extension_required",
                case="flow",
                exact_substitutions=substitutions,
                collision_rates=partial_rates,
                open_interval_checks=partial_checks,
                violated_constraints=("no_single_shear_bulk_magic",),
                assumptions=classification.assumptions,
                minimal_extension=(
                    "separate_shear_bulk_relaxation_and_general_velocity_"
                    "boundary_correction"
                ),
                consumed_evidence=evidence,
                is_conditional=bool(unresolved_flow),
                feasibility_conditions=unresolved_flow or None,
            )
        return ParameterReport(
            status="boundary_correction_required",
            case="flow",
            exact_substitutions=substitutions,
            collision_rates=partial_rates,
            open_interval_checks=partial_checks,
            violated_constraints=("general_velocity_wall_jets",),
            assumptions=classification.assumptions,
            minimal_extension="explicit_general_velocity_boundary_correction",
            consumed_evidence=evidence,
            is_conditional=bool(unresolved_flow),
            feasibility_conditions=unresolved_flow or None,
        )

    residual = velocity_bb_residual(
        sigma_f_plus=sigma_plus_symbol,
        sigma_f_minus=sigma_minus_symbol,
        chi_s=chi_s,
        chi_b=chi_b,
        force_mode="uniform_body_force",
        normalization="source_feedback_henon_product",
        momentum_gauge="half_force",
        geometry="flat_grid_aligned_halfway",
    )
    classification = classify_magic(residual)
    conditions = residual.solve_zero_conditions()
    if classification.status != "restricted_calibration":
        raise ValueError("Task 6 did not classify the restricted shear calibration")
    if len(residual.calibration_parameters) != 1:
        raise ValueError("Task 6 did not expose one restricted shear product")
    product = residual.calibration_parameters[0]
    if product not in conditions:
        raise ValueError("Task 6 did not solve the restricted shear product")
    odd_roots = solve(
        product.subs(sigma_plus_symbol, sigma_plus) - conditions[product],
        sigma_minus_symbol,
        dict=False,
    )
    if len(odd_roots) != 1:
        raise ValueError("restricted shear product did not determine one odd shift")
    sigma_minus = simplify(odd_roots[0])
    rates = {
        "nominal_even": simplify(1 / (sigma_plus + Rational(1, 2))),
        "nominal_odd": simplify(1 / (sigma_minus + Rational(1, 2))),
        "physical_shear": simplify(
            1 / (physical_shear_shift + Rational(1, 2))
        ),
        "physical_bulk": simplify(
            1 / (physical_bulk_shift + Rational(1, 2))
        ),
    }
    checks = {
        "nominal_even_rate": _open_rate_check(rates["nominal_even"]),
        "nominal_odd_rate": _open_rate_check(rates["nominal_odd"]),
        "physical_shear_rate": _open_rate_check(rates["physical_shear"]),
        "physical_bulk_rate": _open_rate_check(rates["physical_bulk"]),
    }
    nu_positive = symbols("nu_positive", positive=True)
    odd_rate_as_nu = simplify(
        1
        / (
            3 * cs2 * dt / (16 * nu_positive)
            + Rational(1, 2)
        )
    )
    substitutions = {
        "nu": nu,
        "sigma_f_plus": sigma_plus,
        "sigma_f_minus": sigma_minus,
        "sigma_shear_physical": physical_shear_shift,
        "sigma_bulk_physical": physical_bulk_shift,
        "nu_bulk_2d": nu_bulk_2d,
        "s_f_minus_as_nu": odd_rate_as_nu,
        "s_f_minus_limit_nu_to_zero_positive": limit(
            odd_rate_as_nu, nu_positive, 0, dir="+"
        ),
    }
    physical, physical_violations, unresolved_physical = _condition_partition(
        {
            "nu_positive": _positive_check(nu),
            "cs2_positive": _positive_check(cs2),
            "dt_positive": _positive_check(dt),
            "one_minus_chi_s_positive": _positive_check(1 - chi_s),
            "one_minus_chi_b_positive": _positive_check(1 - chi_b),
            "sigma_f_plus_positive": _positive_check(sigma_plus),
            "sigma_f_minus_positive": _positive_check(sigma_minus),
            "sigma_shear_physical_positive": _positive_check(
                physical_shear_shift
            ),
            "sigma_bulk_physical_positive": _positive_check(
                physical_bulk_shift
            ),
            "nu_bulk_2d_positive": _positive_check(nu_bulk_2d),
            **checks,
        }
    )
    valid = not physical_violations and not unresolved_physical
    return ParameterReport(
        status="feasible_restricted" if valid else "no_feasible_solution",
        case="flow",
        exact_substitutions=substitutions,
        collision_rates=rates,
        open_interval_checks=checks,
        violated_constraints=physical_violations,
        assumptions=classification.assumptions
        + ("open_interval_does_not_imply_good_conditioning",),
        minimal_extension=None if valid else "change_positive_flow_transport_target",
        consumed_evidence=ConsumedEvidence(
            bulk_route=None,
            bulk_conditions=None,
            boundary_residual=residual,
            boundary_classification=classification,
        ),
        is_conditional=bool(unresolved_physical),
        feasibility_conditions=unresolved_physical or None,
    )


def _open_rate_check(rate: Expr) -> object:
    if rate.is_real is False:
        return False
    if rate.free_symbols and rate.is_real is not True:
        return And(Q.positive(rate), Q.positive(2 - rate))
    try:
        check = And(rate > 0, rate < 2)
    except TypeError:
        return And(Q.positive(rate), Q.positive(2 - rate))
    return bool(check) if check in (True, False) else check


def _three_rates(
    sigma_odd: Expr, sigma_even: Expr, sigma_flux: Expr
) -> Mapping[str, Expr]:
    return {
        "nominal_odd": simplify(1 / (sigma_odd + Rational(1, 2))),
        "nominal_even": simplify(1 / (sigma_even + Rational(1, 2))),
        "physical_flux": simplify(1 / (sigma_flux + Rational(1, 2))),
    }


def _three_rate_checks(rates: Mapping[str, Expr]) -> Mapping[str, object]:
    return {
        "nominal_odd_rate": _open_rate_check(rates["nominal_odd"]),
        "nominal_even_rate": _open_rate_check(rates["nominal_even"]),
        "physical_flux_rate": _open_rate_check(rates["physical_flux"]),
    }


def _truth(expression: object) -> object:
    return bool(expression) if expression in (True, False) else expression


def _positive_check(expression: Expr) -> object:
    expression = sympify(expression)
    if expression.is_positive is True:
        return True
    if expression.is_positive is False:
        return False
    return Q.positive(expression)


def _condition_partition(
    conditions: Mapping[str, object],
) -> tuple[Mapping[str, object], tuple[str, ...], Mapping[str, object]]:
    normalized = {name: _truth(value) for name, value in conditions.items()}
    violations = tuple(
        _failed_condition_name(name)
        for name, value in normalized.items()
        if value is False
    )
    unresolved = {
        name: value
        for name, value in normalized.items()
        if value is not True and value is not False
    }
    return normalized, violations, unresolved


def _failed_condition_name(name: str) -> str:
    if name.endswith("_positive"):
        return f"{name[:-len('_positive')]}_not_positive"
    return name


def _degenerate_scalar_report(
    *,
    case: str,
    a: Expr,
    b: Expr,
    violation: str | tuple[str, ...],
    extension: str,
    extra_substitutions: Mapping[str, Expr] | None = None,
) -> ParameterReport:
    violations = (violation,) if isinstance(violation, str) else violation
    substitutions: dict[str, Expr] = {"a": a, "b": b}
    if extra_substitutions:
        substitutions.update(extra_substitutions)
    return ParameterReport(
        status="degenerate_branch",
        case=case,
        exact_substitutions=substitutions,
        collision_rates={},
        open_interval_checks={},
        violated_constraints=violations,
        assumptions=("classified_before_divided_main_branch",),
        minimal_extension=extension,
        consumed_evidence=ConsumedEvidence(
            bulk_route=None,
            bulk_conditions=None,
        ),
    )


def _degenerate_flow_report(
    *,
    nu: Expr,
    cs2: Expr,
    dt: Expr,
    chi_s: Expr,
    chi_b: Expr,
    violation: str,
) -> ParameterReport:
    return ParameterReport(
        status="degenerate_branch",
        case="flow",
        exact_substitutions={
            "nu": nu,
            "cs2": cs2,
            "dt": dt,
            "chi_s": chi_s,
            "chi_b": chi_b,
        },
        collision_rates={},
        open_interval_checks={},
        violated_constraints=(violation,),
        assumptions=("classified_before_flow_rate_division",),
        minimal_extension="restore_positive_shear_transport_denominator",
        consumed_evidence=ConsumedEvidence(
            bulk_route=None,
            bulk_conditions=None,
        ),
    )
