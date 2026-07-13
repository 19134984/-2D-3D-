"""Exact source-aware wall residuals for the transformed LBM-CDE TRT model.

The module keeps wall Taylor jets as named independent coefficients.  A
reported calibration always carries its normalization, forcing gauge, and
geometry; no bare relaxation product is treated as universal.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from sympy import (
    Expr,
    Matrix,
    Poly,
    Rational,
    Symbol,
    Tuple,
    expand,
    simplify,
    solve,
    sqrt,
    symbols,
    sympify,
)


MAGIC_STATUSES = {
    "universal_magic",
    "restricted_calibration",
    "boundary_correction_required",
    "no_single_magic",
    "corner_closure_conflict",
}


@dataclass(frozen=True)
class BoundaryResidual:
    """Named independent wall-jet coefficients and their derivation context."""

    boundary: str
    coefficients: Mapping[str, Expr]
    assumptions: tuple[str, ...]
    parameter_mapping: Mapping[str, Expr]
    calibration_parameters: tuple[Expr, ...]
    recommended_status: str
    derivation: str

    def solve_zero_conditions(self) -> Mapping[Expr, Expr]:
        """Solve every nonzero retained coefficient for each calibration product."""

        equations = [simplify(value) for value in self.coefficients.values()]
        equations = [value for value in equations if value != 0]
        if not equations:
            return {}
        conditions: dict[Expr, Expr] = {}
        for index, parameter in enumerate(self.calibration_parameters):
            placeholder = Symbol(f"_boundary_product_{index}")
            reduced = [equation.subs(parameter, placeholder) for equation in equations]
            roots = solve(reduced, placeholder, dict=False)
            if isinstance(roots, dict):
                roots = [roots[placeholder]] if placeholder in roots else []
            elif roots and isinstance(roots[0], dict):
                roots = [root[placeholder] for root in roots if placeholder in root]
            if len(roots) == 1:
                conditions[parameter] = simplify(roots[0])
                continue
            for factor_symbol in sorted(
                parameter.free_symbols, key=lambda item: item.name, reverse=True
            ):
                factor_roots = solve(equations, factor_symbol, dict=False)
                if isinstance(factor_roots, dict):
                    factor_roots = (
                        [factor_roots[factor_symbol]]
                        if factor_symbol in factor_roots
                        else []
                    )
                elif factor_roots and isinstance(factor_roots[0], dict):
                    factor_roots = [
                        root[factor_symbol]
                        for root in factor_roots
                        if factor_symbol in root
                    ]
                if len(factor_roots) == 1:
                    conditions[parameter] = simplify(
                        parameter.subs(factor_symbol, factor_roots[0])
                    )
                    break
        return conditions

    def substitute(self, replacements: Mapping[Expr, Expr]) -> "BoundaryResidual":
        """Return a residual with exact substitutions applied to every coefficient."""

        return BoundaryResidual(
            boundary=self.boundary,
            coefficients={
                name: simplify(value.subs(replacements))
                for name, value in self.coefficients.items()
            },
            assumptions=self.assumptions,
            parameter_mapping={
                name: simplify(value.subs(replacements))
                for name, value in self.parameter_mapping.items()
            },
            calibration_parameters=self.calibration_parameters,
            recommended_status=self.recommended_status,
            derivation=self.derivation,
        )

    @property
    def unsatisfied_jets(self) -> tuple[str, ...]:
        return tuple(
            name
            for name, value in self.coefficients.items()
            if simplify(value) != 0
        )

    def has_independent(self, jet: str) -> bool:
        return jet in self.coefficients and simplify(self.coefficients[jet]) != 0

    @property
    def cancelled_jets(self) -> tuple[str, ...]:
        """Rows retained in the audit table that cancel in the stated chain."""

        return tuple(
            name
            for name, value in self.coefficients.items()
            if simplify(value) == 0
        )


@dataclass(frozen=True)
class MagicClassification:
    """Classification plus the evidence needed to audit it."""

    status: str
    coefficients: Mapping[str, Expr]
    assumptions: tuple[str, ...]
    parameter_mapping: Mapping[str, Expr]
    unsatisfied_jets: tuple[str, ...]
    cancelled_jets: tuple[str, ...]
    rate_compatibility_status: str | None = None
    rate_compatibility_conditions: Mapping[Expr, Expr] | None = None


@dataclass(frozen=True)
class ManufacturedStencilCheck:
    """Auditable stencil result or explicitly labeled source-only evidence."""

    profile: str
    first_nonzero_coefficient: Expr
    probes_quadratic_slip: bool
    uses_half_force_reconstruction: bool
    derivation: Expr
    includes_pressure_wall_equilibrium: bool = False
    uses_half_source_reconstruction: bool = False
    verification_status: str = "verified"
    consumed_inputs: tuple[str, ...] = ()
    unresolved_jets: tuple[str, ...] = ()


@dataclass(frozen=True)
class CornerStencilCheck:
    """Rank, overwrite, source-count, and geometry audit for one corner."""

    profile: str
    closure_rank: int
    augmented_rank: int
    dirichlet_then_adiabatic: Expr
    adiabatic_then_dirichlet: Expr
    first_nonzero_coefficient: Expr
    simultaneously_satisfied: bool
    naive_overwrite_heat_source_count: int
    shared_population_heat_source_count: int
    diagonal_wall_distance: Expr
    derivation: Expr
    source_assignment_population_ids: tuple[str, ...] = ()
    shared_source_population_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class PopulationChainResult:
    """Executable population equations and the coefficient table they generate."""

    coefficients: Mapping[str, Expr]
    equations: tuple[Expr, ...]
    assignments: tuple[tuple[str, str, Expr], ...] = ()
    metadata: Mapping[str, Expr] | None = None
    unresolved_jets: tuple[str, ...] = ()


def velocity_bb_residual(
    *,
    sigma_f_plus: Expr,
    sigma_f_minus: Expr,
    chi_s: Expr = 0,
    chi_b: Expr = 0,
    force_mode: str,
    normalization: str,
    momentum_gauge: str,
    geometry: str,
) -> BoundaryResidual:
    """Return the classical halfway-BB quadratic slip residual.

    The uniform-body-force branch is Multireflection Eqs. (41)--(42),
    expressed before any LBM-CDE feedback is introduced.  The pressure branch
    records Luo's distinct pressure-boundary calibration.
    """

    sigma_f_plus = sympify(sigma_f_plus)
    sigma_f_minus = sympify(sigma_f_minus)
    chi_s = sympify(chi_s)
    chi_b = sympify(chi_b)
    product = sigma_f_plus * sigma_f_minus

    if geometry != "flat_grid_aligned_halfway":
        raise ValueError("the classical calibration requires a flat halfway wall")

    common = (
        "d2q9_two_rate_collision",
        "steady_linear_stokes",
        "flat_grid_aligned_halfway_wall",
    )
    if (
        force_mode == "uniform_body_force"
        and normalization == "source_feedback_henon_product"
    ):
        if momentum_gauge != "half_force":
            raise ValueError("source-feedback shear uses half-force momentum")
        shear_shift = simplify((1 - chi_s) * sigma_f_plus)
        candidate = simplify(shear_shift * sigma_f_minus)
        coefficients = {
            "shear_curvature": simplify(Rational(16, 3) * candidate - 1)
        }
        assumptions = common + (
            "uniform_body_force",
            "half_force_momentum_gauge",
            "one_dimensional_shear_only",
            "frozen_source_feedback",
            "no_bulk_time_or_tangential_jets",
        )
        parameter_mapping = _flow_feedback_mapping(
            sigma_f_plus,
            sigma_f_minus,
            chi_s,
            chi_b,
        )
        derivation = (
            "Mode-projected Eq. (42) calibration with the Task 3 physical "
            "shear shift; nominal odd ghost shift is retained"
        )
        calibration_parameters = (candidate,)
        recommended_status = "restricted_calibration"
    elif force_mode == "general_source_aware":
        if normalization != "source_feedback_henon_product":
            raise ValueError("general feedback requires explicit split mapping")
        if momentum_gauge != "half_force":
            raise ValueError("general feedback requires half-force momentum")
        shear_candidate = simplify((1 - chi_s) * product)
        bulk_candidate = simplify((1 - chi_b) * product)
        coefficients = {
            "shear_curvature": simplify(
                Rational(16, 3) * shear_candidate - 1
            ),
            "bulk_curvature": simplify(
                Rational(16, 3) * bulk_candidate - 1
            ),
        }
        unresolved_jets = (
            "normal_velocity_gradient",
            "tangential_velocity_gradient",
            "normal_velocity_curvature",
            "mixed_velocity_curvature",
            "tangential_velocity_curvature",
            "wall_time_velocity",
            "wall_time_normal_velocity",
            "wall_time_tangential_velocity",
            "pressure",
            "normal_pressure_gradient",
            "tangential_pressure_gradient",
            "normal_force",
            "tangential_force",
            "normal_force_gradient",
            "tangential_force_gradient",
            "force_time",
            "normal_velocity_source",
            "tangential_velocity_source",
            "source_time",
        )
        coefficients.update(
            {
                jet: Symbol(f"C_velocity_{jet}_unresolved", nonzero=True)
                for jet in unresolved_jets
            }
        )
        assumptions = common + (
            "source_aware_transformed_collision",
            "half_force_momentum_gauge",
            "independent_shear_and_bulk_jets",
            "independent_pressure_force_time_velocity_and_source_jets",
            "unresolved_full_velocity_wall_chain",
            "unresolved_rows_require_explicit_boundary_correction",
        )
        parameter_mapping = _flow_feedback_mapping(
            sigma_f_plus,
            sigma_f_minus,
            chi_s,
            chi_b,
        )
        derivation = (
            "Independent shear and trace projections of the halfway chain are "
            "resolved; all other retained pressure, force, time, velocity, and "
            "source jets are explicit unresolved coefficients pending a full "
            "population-chain closure"
        )
        calibration_parameters = (product,)
        recommended_status = "boundary_correction_required"
    elif force_mode == "uniform_body_force":
        if normalization != "glt2003":
            raise ValueError("uniform-force quarter mapping requires glt2003")
        if momentum_gauge != "half_force":
            raise ValueError("uniform-force mapping requires half-force momentum")
        if chi_s != 0 or chi_b != 0:
            raise ValueError("the paper normalization audit is the no-feedback limit")
        lambda_squared = simplify(Rational(4, 3) * product)
        coefficients = {"quadratic_slip": simplify(4 * lambda_squared - 1)}
        assumptions = common + (
            "uniform_body_force",
            "half_force_momentum_gauge",
            "no_lbm_cde_feedback",
        )
        parameter_mapping = {
            "lambda_nu": simplify(1 / (sigma_f_plus + Rational(1, 2))),
            "lambda_2": simplify(1 / (sigma_f_minus + Rational(1, 2))),
            "henon_even": sigma_f_plus,
            "henon_odd": sigma_f_minus,
            "Lambda_squared": lambda_squared,
            "henon_product": product,
        }
        derivation = (
            "Multireflection Eq. (41) inserted into Eq. (42): "
            "D_eff^2-D_half^2=4*Lambda^2-1"
        )
        calibration_parameters = (product,)
        recommended_status = "restricted_calibration"
    elif force_mode == "pressure_boundary":
        if normalization != "wang_luo_henon_product":
            raise ValueError("pressure drive uses the Wang/Luo Hénon product")
        if momentum_gauge != "pressure_boundary":
            raise ValueError("pressure drive must not use the body-force gauge")
        if chi_s != 0 or chi_b != 0:
            raise ValueError("Luo's pressure calibration is the no-feedback limit")
        coefficients = {"quadratic_slip": simplify(product - Rational(3, 8))}
        assumptions = common + (
            "pressure_boundary_drive",
            "pressure_boundary_momentum_gauge",
            "no_lbm_cde_feedback",
        )
        parameter_mapping = {
            "henon_product": product,
            "source_specific_value": Rational(3, 8),
        }
        derivation = "Luo 2014 discussion following Eq. (50), pressure drive"
        calibration_parameters = (product,)
        recommended_status = "restricted_calibration"
    else:
        raise ValueError(f"unsupported force mode: {force_mode}")

    return BoundaryResidual(
        boundary="velocity_halfway_bounce_back",
        coefficients=coefficients,
        assumptions=assumptions,
        parameter_mapping=parameter_mapping,
        calibration_parameters=calibration_parameters,
        recommended_status=recommended_status,
        derivation=derivation,
    )


def _flow_feedback_mapping(
    sigma_f_plus: Expr,
    sigma_f_minus: Expr,
    chi_s: Expr,
    chi_b: Expr,
) -> Mapping[str, Expr]:
    """Return physical blocks and untouched nominal TRT ghost shifts."""

    s_f_plus = simplify(1 / (sigma_f_plus + Rational(1, 2)))
    s_f_minus = simplify(1 / (sigma_f_minus + Rational(1, 2)))
    return {
        "flow_shear_physical": simplify((1 - chi_s) * sigma_f_plus),
        "flow_bulk_physical": simplify((1 - chi_b) * sigma_f_plus),
        "flow_even_ghost_nominal": sigma_f_plus,
        "flow_odd_ghost_nominal": sigma_f_minus,
        "s_f_plus_nominal": s_f_plus,
        "s_f_minus_nominal": s_f_minus,
        "flow_even_source_factor": simplify(1 - s_f_plus / 2),
        "flow_odd_source_factor": simplify(1 - s_f_minus / 2),
        "momentum_half_source": Rational(1, 2),
    }


def temperature_abb_residual(
    *,
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
    regime: str,
    normalization: str,
    scalar_gauge: str,
    geometry: str,
    normal_temperature_gradient: Expr = 0,
) -> BoundaryResidual:
    """Return D2Q9 Dirichlet ABB coefficients from the transformed link chain.

    The restricted row is obtained by solving the steady D1Q3 aggregate of
    the three normal D2Q9 links with a quadratic temperature polynomial.  The
    general table keeps the remaining Taylor jets instead of applying the CDE
    to identify them prematurely.
    """

    sigma_g_plus, sigma_g_minus, chi_kappa = map(
        sympify, (sigma_g_plus, sigma_g_minus, chi_kappa)
    )
    cs2, pressure_ratio = map(sympify, (cs2, pressure_ratio))
    normal_temperature_gradient = sympify(normal_temperature_gradient)
    _validate_scalar_boundary_context(
        cs2=cs2,
        normalization=normalization,
        scalar_gauge=scalar_gauge,
        geometry=geometry,
    )
    mapping = _scalar_feedback_mapping(
        sigma_g_plus,
        sigma_g_minus,
        chi_kappa,
        cs2,
        pressure_ratio,
    )
    physical_product = mapping["scalar_flux_physical"] * sigma_g_plus
    chain = _temperature_abb_population_chain(
        sigma_g_plus,
        sigma_g_minus,
        chi_kappa,
        cs2,
        pressure_ratio,
        regime,
        normal_temperature_gradient,
    )
    common = (
        "d2q9_two_rate_scalar_collision",
        "flat_grid_aligned_halfway_wall",
        "full_pressure_equilibrium_wall_term",
        "half_source_temperature_reconstruction",
        "parity_specific_scalar_source_factors",
    )

    if regime == "steady_1d_quadratic":
        coefficients = chain.coefficients
        assumptions = common + (
            "steady_one_dimensional_quadratic_temperature",
            "constant_pressure_ratio",
            "zero_velocity_and_force",
            "cde_consistent_uniform_heat_source",
            "no_tangential_or_wall_time_jets",
            "external_exact_gradient_source_realization",
        )
        calibration_parameters = (physical_product,)
        status = "restricted_calibration"
        derivation = (
            "Exact D1Q3 aggregate of the three normal D2Q9 population "
            "recurrences; ABB uses 2*(w_n+lambda_n*pi/cs2)*T_w"
        )
    elif regime == "steady_1d_quadratic_local_feedback":
        coefficients = chain.coefficients
        assumptions = common + (
            "steady_one_dimensional_quadratic_temperature",
            "constant_pressure_ratio",
            "zero_velocity_and_force",
            "cde_consistent_uniform_heat_source",
            "no_tangential_or_wall_time_jets",
            "local_nonequilibrium_feedback",
            "homogeneous_physical_scalar_flux_block",
            "nominal_even_ghost_block",
        )
        calibration_parameters = (physical_product,)
        status = "restricted_calibration"
        derivation = (
            "Exact homogeneous D1Q3 recurrence after local feedback "
            "elimination: physical scalar-flux shift and nominal even ghost"
        )
    elif regime == "general_source_aware":
        coefficients = chain.coefficients
        assumptions = common + (
            "independent_wall_taylor_jets",
            "normal_pressure_row_comes_from_affine_equilibrium_product_chain",
            "bulk_cde_not_substituted_into_boundary_table",
            "external_exact_gradient_source_rows_not_extrapolated_to_local_feedback",
            "unresolved_rows_require_explicit_boundary_correction",
        )
        calibration_parameters = (physical_product,)
        status = "boundary_correction_required"
        derivation = (
            "Quadratic and time-dependent normal-link recurrences plus an "
            "independent affine T*pi hydrostatic recurrence; rate-independent "
            "pressure/source/time rows are retained"
        )
    else:
        raise ValueError(f"unsupported temperature ABB regime: {regime}")

    return BoundaryResidual(
        boundary="temperature_dirichlet_anti_bounce_back",
        coefficients=coefficients,
        assumptions=assumptions,
        parameter_mapping=mapping,
        calibration_parameters=calibration_parameters,
        recommended_status=status,
        derivation=derivation,
    )


def adiabatic_bb_residual(
    *,
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
    regime: str,
    normalization: str,
    scalar_gauge: str,
    geometry: str,
) -> BoundaryResidual:
    """Return kinetic odd-flux rows for D2Q9 adiabatic bounce-back."""

    sigma_g_plus, sigma_g_minus, chi_kappa = map(
        sympify, (sigma_g_plus, sigma_g_minus, chi_kappa)
    )
    cs2, pressure_ratio = map(sympify, (cs2, pressure_ratio))
    _validate_scalar_boundary_context(
        cs2=cs2,
        normalization=normalization,
        scalar_gauge=scalar_gauge,
        geometry=geometry,
    )
    mapping = _scalar_feedback_mapping(
        sigma_g_plus,
        sigma_g_minus,
        chi_kappa,
        cs2,
        pressure_ratio,
    )
    chain = _adiabatic_population_chain(
        sigma_g_plus,
        sigma_g_minus,
        chi_kappa,
        cs2,
        pressure_ratio,
        regime,
    )
    common = (
        "d2q9_two_rate_scalar_collision",
        "flat_grid_aligned_halfway_wall",
        "kinetic_reflected_odd_flux_derived_before_neumann_data",
        "half_source_temperature_reconstruction",
        "parity_specific_scalar_source_factors",
    )
    if regime == "kinetic_reflected_flux":
        coefficients = chain.coefficients
        assumptions = common + (
            "steady_one_dimensional_affine_temperature",
            "constant_pressure_ratio",
            "zero_tangential_jets",
        )
        status = "restricted_calibration"
        derivation = (
            "Exact normal D1Q3 reflected-link recurrence; the row is the "
            "kinetic odd flux and is prescribed zero, not tuned by a product"
        )
    elif regime == "general_source_aware":
        coefficients = chain.coefficients
        assumptions = common + (
            "smooth_homogeneous_neumann_data_imposed_after_chain",
            "axis_and_two_diagonal_links_summed",
            "force_only_and_hydrostatic_pressure_force_rows_kept_separate",
            "unresolved_time_normal_curvature_and_source_rows_are_explicit",
            "unresolved_rows_require_explicit_boundary_correction",
        )
        status = "boundary_correction_required"
        derivation = (
            "Normal link plus both diagonal reflected recurrences; diagonal "
            "tangential gradients cancel but their curvature sum remains"
        )
    else:
        raise ValueError(f"unsupported adiabatic BB regime: {regime}")

    return BoundaryResidual(
        boundary="temperature_adiabatic_bounce_back",
        coefficients=coefficients,
        assumptions=assumptions,
        parameter_mapping=mapping,
        calibration_parameters=(),
        recommended_status=status,
        derivation=derivation,
    )


def mixed_corner_residual(
    *,
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
    grid_spacing: Expr,
    corner_type: str,
    normalization: str,
    scalar_gauge: str,
    geometry: str,
) -> BoundaryResidual:
    """Return the compatibility defect of one shared diagonal corner value."""

    sigma_g_plus, sigma_g_minus, chi_kappa = map(
        sympify, (sigma_g_plus, sigma_g_minus, chi_kappa)
    )
    cs2, pressure_ratio, grid_spacing = map(
        sympify, (cs2, pressure_ratio, grid_spacing)
    )
    if corner_type != "dirichlet_adiabatic":
        raise ValueError("only the mixed Dirichlet/adiabatic corner is audited")
    if geometry != "grid_aligned_right_angle_corner":
        raise ValueError("corner audit requires a grid-aligned right angle")
    _validate_scalar_boundary_context(
        cs2=cs2,
        normalization=normalization,
        scalar_gauge=scalar_gauge,
        geometry="flat_grid_aligned_halfway",
    )
    chain = _corner_population_chain(
        sigma_g_plus,
        sigma_g_minus,
        chi_kappa,
        cs2,
        pressure_ratio,
        grid_spacing,
    )
    coefficients = chain.coefficients
    mapping = dict(
        _scalar_feedback_mapping(
            sigma_g_plus,
            sigma_g_minus,
            chi_kappa,
            cs2,
            pressure_ratio,
        )
    )
    mapping.update(chain.metadata or {})
    return BoundaryResidual(
        boundary="mixed_dirichlet_adiabatic_corner",
        coefficients=coefficients,
        assumptions=(
            "one_diagonal_population_shared_by_two_walls",
            "dirichlet_abb_and_adiabatic_bb_both_constrain_shared_unknown",
            "grid_aligned_right_angle_corner",
            "corner_taylor_polynomial_evaluated_at_half_diagonal_link",
            "sequential_overwrite_is_not_simultaneous_satisfaction",
            "shared_population_source_counted_once",
            "unresolved_corner_source_and_time_rows_are_explicit",
        ),
        parameter_mapping=mapping,
        calibration_parameters=(),
        recommended_status="corner_closure_conflict",
        derivation=(
            "A=[1;1] acts on one shared diagonal population; generic distinct "
            "Dirichlet-ABB and adiabatic-BB right-hand sides give rank([A|b])=2"
        ),
    )


def _validate_scalar_boundary_context(
    *, cs2: Expr, normalization: str, scalar_gauge: str, geometry: str
) -> None:
    if cs2 != Rational(1, 3):
        raise ValueError("the audited D2Q9 chain requires cs2=1/3")
    if normalization != "transformed_cde_chain":
        raise ValueError("scalar boundary audit requires transformed_cde_chain")
    if scalar_gauge != "half_source":
        raise ValueError("scalar boundary audit requires half-source reconstruction")
    if geometry != "flat_grid_aligned_halfway":
        raise ValueError("scalar boundary audit requires a flat halfway wall")


def _scalar_feedback_mapping(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
) -> Mapping[str, Expr]:
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    s_minus = simplify(1 / (sigma_g_minus + Rational(1, 2)))
    return {
        "scalar_flux_physical": (
            (1 - chi_kappa)
            * cs2
            / (cs2 + pressure_ratio)
            * sigma_g_minus
        ),
        "scalar_even_ghost_nominal": sigma_g_plus,
        "scalar_odd_ghost_nominal": sigma_g_minus,
        "s_g_plus_nominal": s_plus,
        "s_g_minus_nominal": s_minus,
        "scalar_even_source_factor": simplify(1 - s_plus / 2),
        "scalar_odd_source_factor": simplify(1 - s_minus / 2),
        "scalar_half_source": Rational(1, 2),
        "pressure_equilibrium_normal_weight": simplify(
            Rational(1, 6) + pressure_ratio / 2
        ),
    }


def _temperature_abb_population_chain(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
    regime: str,
    normal_temperature_gradient: Expr,
) -> PopulationChainResult:
    """Generate ABB rows from finite population recurrences, not a final table."""

    local_feedback = regime == "steady_1d_quadratic_local_feedback"
    if regime not in {
        "steady_1d_quadratic",
        "steady_1d_quadratic_local_feedback",
        "general_source_aware",
    }:
        raise ValueError(f"unsupported temperature ABB regime: {regime}")
    homogeneous_flux_shift = None
    if local_feedback:
        homogeneous_flux_shift = simplify(
            (1 - chi_kappa)
            * cs2
            / (cs2 + pressure_ratio)
            * sigma_g_minus
        )
    curvature, equations = _primary_quadratic_abb_population_equations(
        sigma_g_plus,
        sigma_g_minus,
        chi_kappa,
        pressure_ratio,
        homogeneous_flux_shift,
    )
    if regime != "general_source_aware":
        return PopulationChainResult(
            coefficients={"temperature_curvature": curvature},
            equations=equations,
        )

    pressure_coefficient, pressure_equation = _primary_affine_pressure_chain(
        sympify(normal_temperature_gradient)
    )
    unresolved = (
        "temperature",
        "normal_temperature_gradient",
        "tangential_temperature_gradient",
        "mixed_temperature_curvature",
        "tangential_temperature_curvature",
        "wall_time",
        "wall_time_time",
        "wall_time_normal",
        "wall_time_tangential",
        "pressure",
        "tangential_pressure_gradient",
        "normal_force",
        "tangential_force",
        "heat_source",
        "normal_heat_source_gradient",
        "tangential_heat_source_gradient",
        "normal_velocity",
        "tangential_velocity",
        "pressure_temperature_gradient_source",
        "temperature_force_source",
        "heat_velocity_source",
    )
    coefficients = {
        name: Symbol(f"C_abb_{name}_unresolved", nonzero=True)
        for name in unresolved
    }
    coefficients.update(
        {
            "temperature_curvature": curvature,
            "normal_pressure_gradient": pressure_coefficient,
        }
    )
    return PopulationChainResult(
        coefficients=coefficients,
        equations=equations + (pressure_equation,),
        unresolved_jets=unresolved,
    )


def _primary_quadratic_abb_population_equations(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    pressure_ratio: Expr,
    homogeneous_flux_shift: Expr | None,
) -> tuple[Expr, tuple[Expr, ...]]:
    """Solve a finite D1Q3 aggregate used only by the primary ABB route."""

    x = Symbol("x_primary_abb")
    temperature_0, temperature_1 = symbols("T0_primary T1_primary")
    temperature_2 = Symbol("T2_primary", nonzero=True)
    heat = Symbol("Q_primary")
    temperature = temperature_0 + temperature_1 * x + temperature_2 * x**2
    gradient = temperature.diff(x)
    populations = symbols("primary_q0:9")
    g_zero = sum(populations[index] * x**index for index in range(3))
    g_plus = sum(populations[3 + index] * x**index for index in range(3))
    g_minus = sum(populations[6 + index] * x**index for index in range(3))
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    odd_shift = (
        sigma_g_minus
        if homogeneous_flux_shift is None
        else sympify(homogeneous_flux_shift)
    )
    s_minus = simplify(1 / (odd_shift + Rational(1, 2)))
    even = (g_plus + g_minus) / 2
    odd = (g_plus - g_minus) / 2
    equilibrium_zero = (Rational(2, 3) - pressure_ratio) * temperature
    equilibrium_normal = (
        Rational(1, 6) + pressure_ratio / 2
    ) * temperature
    odd_source = (
        Rational(1, 6)
        * (pressure_ratio + chi_kappa * Rational(1, 3))
        / Rational(1, 3)
        * gradient
        if homogeneous_flux_shift is None
        else sympify(0)
    )
    post_zero = (
        g_zero
        - s_plus * (g_zero - equilibrium_zero)
        + (1 - s_plus / 2) * Rational(2, 3) * heat
    )
    post_plus = (
        g_plus
        - s_plus * (even - equilibrium_normal)
        - s_minus * odd
        + (1 - s_plus / 2) * Rational(1, 6) * heat
        + (1 - s_minus / 2) * odd_source
    )
    post_minus = (
        g_minus
        - s_plus * (even - equilibrium_normal)
        + s_minus * odd
        + (1 - s_plus / 2) * Rational(1, 6) * heat
        - (1 - s_minus / 2) * odd_source
    )
    equations = tuple(
        _polynomial_coefficients(
            (
                g_zero - post_zero,
                g_plus.subs(x, x + 1) - post_plus,
                g_minus.subs(x, x - 1) - post_minus,
                g_zero + g_plus + g_minus - (temperature - heat / 2),
            ),
            (x,),
        )
    )
    solution = solve(
        equations, (*populations, heat), dict=True, simplify=False
    )[0]
    wall_equilibrium = (
        Rational(1, 6) + pressure_ratio / 2
    ) * temperature_0
    link_residual = simplify(
        (
            g_plus.subs(x, Rational(1, 2))
            + post_minus.subs(x, Rational(1, 2))
            - 2 * wall_equilibrium
        ).subs(solution)
    )
    return simplify(link_residual / temperature_2), equations


def _primary_affine_pressure_chain(
    normal_temperature_gradient: Expr,
) -> tuple[Expr, Expr]:
    """Extract the affine T*pi cross jet from the finite wall equilibrium."""

    x = Symbol("x_primary_pressure")
    h = Symbol("h_primary_pressure", nonzero=True)
    temperature_0, pressure_0 = symbols("T0_primary_pressure pi0_primary")
    pressure_gradient = Symbol("pi_n_primary", nonzero=True)
    temperature = temperature_0 + normal_temperature_gradient * x
    pressure = pressure_0 + pressure_gradient * x
    equilibrium = (Rational(1, 6) + pressure / 2) * temperature
    wall_pair_defect = simplify(
        -(
            equilibrium.subs(x, h / 2)
            + equilibrium.subs(x, -h / 2)
            - 2 * equilibrium.subs(x, 0)
        )
    )
    coefficient = simplify(
        expand(wall_pair_defect).coeff(pressure_gradient) / h**2
    )
    return coefficient, wall_pair_defect


def _adiabatic_population_chain(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
    regime: str,
) -> PopulationChainResult:
    """Solve collision-streaming-reflection equations for normal/diagonal links."""

    if regime not in {"kinetic_reflected_flux", "general_source_aware"}:
        raise ValueError(f"unsupported adiabatic BB regime: {regime}")
    normal_flux, normal_equations, normal_residual = (
        _primary_adiabatic_normal_population_equations(
            sigma_g_plus,
            sigma_g_minus,
            chi_kappa,
            cs2,
            pressure_ratio,
        )
    )
    if regime == "kinetic_reflected_flux":
        return PopulationChainResult(
            coefficients={"normal_temperature_gradient": normal_flux},
            equations=normal_equations,
            assignments=(("normal_pair", "bounce_back", normal_residual),),
        )

    tangential_curvature, diagonal_equations, diagonal_residuals = (
        _primary_adiabatic_diagonal_population_equations(
            sigma_g_plus,
            sigma_g_minus,
            chi_kappa,
            cs2,
            pressure_ratio,
        )
    )
    unresolved = (
        "temperature_curvature",
        "wall_time",
        "wall_time_normal",
        "wall_time_tangential",
        "pressure",
        "normal_pressure_gradient",
        "tangential_pressure_gradient",
        "normal_force",
        "tangential_force",
        "heat_source",
        "normal_heat_source_gradient",
        "tangential_heat_source_gradient",
        "temperature_force_source",
    )
    coefficients = {
        name: Symbol(f"C_adiabatic_{name}_unresolved", nonzero=True)
        for name in unresolved
    }
    coefficients.update(
        {
            "normal_temperature_gradient": normal_flux,
            "tangential_temperature_curvature": tangential_curvature,
        }
    )
    return PopulationChainResult(
        coefficients=coefficients,
        equations=normal_equations + diagonal_equations,
        assignments=(
            ("normal_pair", "bounce_back", normal_residual),
            *(
                (f"diagonal_{direction}", "bounce_back", residual)
                for direction, residual in zip((-1, 1), diagonal_residuals)
            ),
        ),
        unresolved_jets=unresolved,
    )


def _primary_adiabatic_normal_population_equations(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
) -> tuple[Expr, tuple[Expr, ...], Expr]:
    """Solve a steady affine normal pair and extract its reflected odd flux."""

    x = Symbol("normal_primary_adiabatic")
    temperature_0 = Symbol("T0_normal_primary")
    temperature_1 = Symbol("T1_normal_primary", nonzero=True)
    temperature = temperature_0 + temperature_1 * x
    populations = symbols("normal_primary_population_0:4")
    g_plus = populations[0] + populations[1] * x
    g_minus = populations[2] + populations[3] * x
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    s_minus = simplify(1 / (sigma_g_minus + Rational(1, 2)))
    even = (g_plus + g_minus) / 2
    odd = (g_plus - g_minus) / 2
    equilibrium = (Rational(1, 6) + pressure_ratio / 2) * temperature
    odd_source = simplify(
        Rational(1, 6)
        * (pressure_ratio / cs2 + chi_kappa)
        * temperature.diff(x)
    )
    post_plus = (
        g_plus
        - s_plus * (even - equilibrium)
        - s_minus * odd
        + (1 - s_minus / 2) * odd_source
    )
    post_minus = (
        g_minus
        - s_plus * (even - equilibrium)
        + s_minus * odd
        - (1 - s_minus / 2) * odd_source
    )
    equations = tuple(
        _polynomial_coefficients(
            (
                g_plus.subs(x, x + 1) - post_plus,
                g_minus.subs(x, x - 1) - post_minus,
            ),
            (x,),
        )
    )
    solution = solve(equations, populations, dict=True, simplify=False)[0]
    reflected_residual = simplify((g_plus - post_minus).subs(solution))
    coefficient = simplify(reflected_residual / temperature_1)
    return coefficient, equations, reflected_residual


def _primary_adiabatic_diagonal_population_equations(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
) -> tuple[Expr, tuple[Expr, ...], tuple[Expr, ...]]:
    """Solve both diagonal polynomial recurrences and sum reflected defects."""

    x = Symbol("tangential_primary_adiabatic")
    temperature_0, temperature_1 = symbols("T0_diagonal_primary T1_diagonal_primary")
    temperature_2 = Symbol("T2_diagonal_primary", nonzero=True)
    temperature = temperature_0 + temperature_1 * x + temperature_2 * x**2
    gradient = temperature.diff(x)
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    s_minus = simplify(1 / (sigma_g_minus + Rational(1, 2)))
    equations: list[Expr] = []
    reflected_residuals: list[Expr] = []
    for direction, label in ((-1, "minus"), (1, "plus")):
        populations = symbols(f"diagonal_primary_{label}_0:6")
        g_plus = sum(populations[index] * x**index for index in range(3))
        g_minus = sum(populations[3 + index] * x**index for index in range(3))
        even = (g_plus + g_minus) / 2
        odd = (g_plus - g_minus) / 2
        equilibrium = (
            Rational(1, 36) + pressure_ratio / 12
        ) * temperature
        odd_source = simplify(
            Rational(1, 36)
            * direction
            * (pressure_ratio / cs2 + chi_kappa)
            * gradient
        )
        post_plus = (
            g_plus
            - s_plus * (even - equilibrium)
            - s_minus * odd
            + (1 - s_minus / 2) * odd_source
        )
        post_minus = (
            g_minus
            - s_plus * (even - equilibrium)
            + s_minus * odd
            - (1 - s_minus / 2) * odd_source
        )
        link_equations = tuple(
            _polynomial_coefficients(
                (
                    g_plus.subs(x, x + direction) - post_plus,
                    g_minus.subs(x, x - direction) - post_minus,
                ),
                (x,),
            )
        )
        solution = solve(
            link_equations, populations, dict=True, simplify=False
        )[0]
        equations.extend(link_equations)
        reflected_residuals.append(
            simplify((g_plus - post_minus).subs(solution))
        )
    combined = simplify(sum(reflected_residuals))
    return (
        simplify(combined / temperature_2),
        tuple(equations),
        tuple(reflected_residuals),
    )


def _corner_population_chain(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
    grid_spacing: Expr,
) -> PopulationChainResult:
    """Generate the competing finite assignments for one shared diagonal link."""

    _ = (sigma_g_minus, chi_kappa, cs2)
    temperature, tx, ty, txx, txy, tyy = symbols(
        "T_corner Tx_corner Ty_corner Txx_corner Txy_corner Tyy_corner"
    )
    diagonal_weight = simplify(Rational(1, 36) + pressure_ratio / 12)
    fluid_temperature = (
        temperature
        + grid_spacing * (tx + ty) / 2
        + grid_spacing**2 * (txx + 2 * txy + tyy) / 8
    )
    outgoing = simplify(diagonal_weight * fluid_temperature)
    wall_equilibrium = simplify(diagonal_weight * temperature)
    dirichlet_rhs = simplify(-outgoing + 2 * wall_equilibrium)
    adiabatic_rhs = outgoing
    incoming = Symbol("g_corner_shared_incoming")
    equations = (
        simplify(incoming - dirichlet_rhs),
        simplify(incoming - adiabatic_rhs),
    )
    compatibility_defect = simplify(dirichlet_rhs - adiabatic_rhs)
    resolved_jets = {
        "temperature": temperature,
        "x_temperature_gradient": tx,
        "y_temperature_gradient": ty,
        "xx_temperature_curvature": txx,
        "mixed_temperature_curvature": txy,
        "yy_temperature_curvature": tyy,
    }
    coefficients = {
        name: simplify(expand(compatibility_defect).coeff(jet))
        for name, jet in resolved_jets.items()
    }
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    heat = Symbol("Q_corner_primary")
    source_increment = simplify(
        Rational(1, 36) * (1 - s_plus / 2) * heat
    )
    coefficients["heat_source"] = simplify(source_increment.coeff(heat))
    unresolved = (
        "wall_time",
        "wall_time_x",
        "wall_time_y",
        "pressure",
        "x_pressure_gradient",
        "y_pressure_gradient",
        "x_force",
        "y_force",
        "x_velocity",
        "y_velocity",
    )
    coefficients.update(
        {
            name: Symbol(f"C_corner_{name}_unresolved", nonzero=True)
            for name in unresolved
        }
    )
    assignments = (
        ("shared_diagonal", "dirichlet_abb", source_increment),
        ("shared_diagonal", "adiabatic_bb", source_increment),
    )
    shared_source_populations = {
        population for population, _, source in assignments if source.has(heat)
    }
    source_applications = [source for _, _, source in assignments if source.has(heat)]
    shared_unknowns = sorted({population for population, _, _ in assignments})
    compatibility_matrix = Matrix([[1] for _ in assignments])
    rhs_dirichlet, rhs_adiabatic = symbols(
        "corner_rhs_dirichlet corner_rhs_adiabatic"
    )
    augmented = compatibility_matrix.row_join(
        Matrix([rhs_dirichlet, rhs_adiabatic])
    )
    metadata = {
        "shared_diagonal_unknown_count": sympify(len(shared_unknowns)),
        "wall_equation_count": sympify(len(equations)),
        "closure_rank": sympify(compatibility_matrix.rank()),
        "generic_augmented_rank": sympify(augmented.rank()),
        "diagonal_equilibrium_wall_weight": diagonal_weight,
        "diagonal_wall_distance": simplify(sqrt(2) * grid_spacing / 2),
        "naive_overwrite_heat_source_count": sympify(len(source_applications)),
        "shared_population_heat_source_count": sympify(
            len(shared_source_populations)
        ),
    }
    return PopulationChainResult(
        coefficients=coefficients,
        equations=equations,
        assignments=assignments,
        metadata=metadata,
        unresolved_jets=unresolved,
    )


def manufactured_velocity_stencil(
    *,
    profile: str,
    geometry: str,
    sigma_f_plus: Expr = 0,
    sigma_f_minus: Expr = 0,
    force_amplitude: Expr = 0,
    pressure_gradient: Expr = 0,
    normalization: str = "",
    momentum_gauge: str = "",
) -> ManufacturedStencilCheck:
    """Return direct polynomial checks or labeled source-only calibrations.

    The Poiseuille entries are paper-evidence records, not independently
    regenerated quadratic-slip probes; their status fields make that explicit.
    """

    if geometry != "flat_grid_aligned_halfway":
        raise ValueError("manufactured checks currently use a halfway wall")
    h = Symbol("h", nonzero=True)
    if profile == "linear_couette":
        wall_value, slope = symbols("u_wall shear")
        x = Symbol("x")
        velocity = wall_value + slope * x
        reflection = simplify(
            velocity.subs(x, -h / 2)
            + velocity.subs(x, h / 2)
            - 2 * wall_value
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=reflection,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=reflection,
        )
    if profile == "hydrostatic_rest":
        force_normal = Symbol("F_n")
        transformed_momentum = -force_normal / 2
        reconstructed_velocity = simplify(
            transformed_momentum + force_normal / 2
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=reconstructed_velocity,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=True,
            derivation=reconstructed_velocity,
        )
    if profile == "uniform_force_quadratic_poiseuille":
        if normalization != "glt2003" or momentum_gauge != "half_force":
            raise ValueError("uniform-force Poiseuille requires the GLT half-force gauge")
        sigma_f_plus, sigma_f_minus = map(
            sympify, (sigma_f_plus, sigma_f_minus)
        )
        lambda_squared = simplify(
            Rational(4, 3) * sigma_f_plus * sigma_f_minus
        )
        coefficient = simplify(4 * lambda_squared - 1)
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=coefficient,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=True,
            derivation=coefficient,
            verification_status="source_evidence_only",
        )
    if profile == "pressure_driven_quadratic_poiseuille":
        if (
            normalization != "wang_luo_henon_product"
            or momentum_gauge != "pressure_boundary"
        ):
            raise ValueError("pressure Poiseuille requires the Wang/Luo pressure gauge")
        sigma_f_plus, sigma_f_minus = map(
            sympify, (sigma_f_plus, sigma_f_minus)
        )
        coefficient = simplify(
            sigma_f_plus * sigma_f_minus - Rational(3, 8)
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=coefficient,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=coefficient,
            includes_pressure_wall_equilibrium=True,
            verification_status="source_evidence_only",
        )
    raise ValueError(f"unsupported velocity profile: {profile}")


def manufactured_temperature_stencil(
    *,
    profile: str,
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
    geometry: str,
    normal_temperature_gradient: Expr = 0,
    normal_pressure_gradient: Expr = 0,
    normal_force: Expr = 0,
    heat_source: Expr = 0,
    wall_time_derivative: Expr = 0,
    normal_temperature_curvature: Expr = 0,
) -> ManufacturedStencilCheck:
    """Direct polynomial link checks independent of the residual APIs."""

    sigma_g_plus, sigma_g_minus, chi_kappa = map(
        sympify, (sigma_g_plus, sigma_g_minus, chi_kappa)
    )
    cs2, pressure_ratio = map(sympify, (cs2, pressure_ratio))
    if geometry != "flat_grid_aligned_halfway":
        raise ValueError("manufactured scalar checks require a halfway wall")
    if cs2 != Rational(1, 3):
        raise ValueError("manufactured D2Q9 checks require cs2=1/3")

    if profile == "dirichlet_nonstationary_abb":
        wall_time_derivative = sympify(wall_time_derivative)
        unresolved = Symbol("C_direct_abb_wall_time_unresolved", nonzero=True)
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=unresolved,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=Tuple(wall_time_derivative, unresolved),
            uses_half_source_reconstruction=True,
            verification_status="unsupported_unresolved",
            consumed_inputs=("wall_time_derivative",),
            unresolved_jets=("wall_time",),
        )
    if profile == "adiabatic_wall_time_normal_curvature":
        wall_time_derivative = sympify(wall_time_derivative)
        normal_temperature_curvature = sympify(normal_temperature_curvature)
        unresolved = Symbol(
            "C_direct_adiabatic_time_normal_curvature_unresolved",
            nonzero=True,
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=unresolved,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=Tuple(
                wall_time_derivative,
                normal_temperature_curvature,
                unresolved,
            ),
            uses_half_source_reconstruction=True,
            verification_status="unsupported_unresolved",
            consumed_inputs=(
                "wall_time_derivative",
                "normal_temperature_curvature",
            ),
            unresolved_jets=("wall_time", "temperature_curvature"),
        )

    if profile == "dirichlet_quadratic":
        coefficient, chain = _direct_quadratic_abb_stencil(
            sigma_g_plus,
            sigma_g_minus,
            chi_kappa,
            pressure_ratio,
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=coefficient,
            probes_quadratic_slip=True,
            uses_half_force_reconstruction=False,
            derivation=chain,
            includes_pressure_wall_equilibrium=True,
            uses_half_source_reconstruction=True,
        )
    if profile == "dirichlet_quadratic_local_feedback":
        scalar_flux = (
            (1 - chi_kappa)
            * cs2
            / (cs2 + pressure_ratio)
            * sigma_g_minus
        )
        coefficient, chain = _direct_quadratic_abb_stencil(
            sigma_g_plus,
            sigma_g_minus,
            chi_kappa,
            pressure_ratio,
            homogeneous_flux_shift=scalar_flux,
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=coefficient,
            probes_quadratic_slip=True,
            uses_half_force_reconstruction=False,
            derivation=chain,
            includes_pressure_wall_equilibrium=True,
            uses_half_source_reconstruction=True,
        )
    if profile == "dirichlet_affine":
        h = Symbol("h", nonzero=True)
        wall_temperature, normal_slope = symbols("T_wall T_n")
        x = Symbol("x")
        temperature = wall_temperature + normal_slope * x
        midpoint_pair = simplify(
            temperature.subs(x, -h / 2)
            + temperature.subs(x, h / 2)
            - 2 * wall_temperature
        )
        wall_weight = simplify(Rational(1, 6) + pressure_ratio / 2)
        chain = simplify(wall_weight * midpoint_pair)
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=chain,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=chain,
            includes_pressure_wall_equilibrium=True,
            uses_half_source_reconstruction=True,
        )
    if profile == "dirichlet_variable_pressure":
        residual = _direct_variable_pressure_abb_stencil(
            sigma_g_plus,
            sigma_g_minus,
            chi_kappa,
            sympify(normal_temperature_gradient),
            sympify(normal_pressure_gradient),
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=residual,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=residual,
            includes_pressure_wall_equilibrium=True,
            uses_half_source_reconstruction=True,
        )
    if profile == "adiabatic_tangential_quadratic":
        coefficient, chain = _direct_adiabatic_diagonal_stencil(
            sigma_g_plus,
            sigma_g_minus,
            chi_kappa,
            pressure_ratio,
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=coefficient,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=chain,
            uses_half_source_reconstruction=True,
        )
    if profile == "adiabatic_force_only":
        residual, temperature, force = _direct_adiabatic_force_stencil(
            sigma_g_plus,
            sigma_g_minus,
            sympify(normal_force),
            sympify(0),
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=simplify(residual / (temperature * force)),
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=residual,
            uses_half_source_reconstruction=True,
        )
    if profile == "adiabatic_hydrostatic_pair":
        residual, _, _ = _direct_adiabatic_force_stencil(
            sigma_g_plus,
            sigma_g_minus,
            sympify(normal_force),
            sympify(normal_pressure_gradient),
        )
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=residual,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=residual,
            uses_half_source_reconstruction=True,
        )
    if profile == "adiabatic_uniform_heat":
        # Q is consumed by both populations with the same even-source factor;
        # their difference proves that it cancels from the reflected odd flux.
        heat_source = sympify(heat_source)
        s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
        even_before = Symbol("g_even_uniform_heat")
        source_increment = simplify(
            (1 - s_plus / 2) * Rational(1, 6) * heat_source
        )
        post_plus = simplify(even_before + source_increment)
        post_minus = simplify(even_before + source_increment)
        odd_after = simplify(post_plus - post_minus)
        return ManufacturedStencilCheck(
            profile=profile,
            first_nonzero_coefficient=odd_after,
            probes_quadratic_slip=False,
            uses_half_force_reconstruction=False,
            derivation=Tuple(post_plus, post_minus, odd_after),
            uses_half_source_reconstruction=True,
            consumed_inputs=("heat_source",),
        )
    raise ValueError(f"unsupported temperature profile: {profile}")


def manufactured_corner_stencil(
    *,
    profile: str,
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    cs2: Expr,
    pressure_ratio: Expr,
    grid_spacing: Expr,
    mixed_temperature_curvature: Expr,
    geometry: str,
) -> CornerStencilCheck:
    """Directly apply both corner closures to the same diagonal unknown."""

    sigma_g_plus, sigma_g_minus, chi_kappa = map(
        sympify, (sigma_g_plus, sigma_g_minus, chi_kappa)
    )
    cs2, pressure_ratio, grid_spacing, mixed_temperature_curvature = map(
        sympify,
        (cs2, pressure_ratio, grid_spacing, mixed_temperature_curvature),
    )
    if profile != "mixed_txy":
        raise ValueError("manufactured corner requires the mixed_txy profile")
    if cs2 != Rational(1, 3):
        raise ValueError("manufactured corner requires D2Q9 cs2=1/3")
    if geometry != "grid_aligned_right_angle_corner":
        raise ValueError("manufactured corner requires a right-angle corner")

    wall_temperature = Symbol("T_corner_wall")
    diagonal_weight = simplify(Rational(1, 36) + pressure_ratio / 12)
    fluid_temperature = simplify(
        wall_temperature
        + mixed_temperature_curvature * grid_spacing**2 / 4
    )
    outgoing_post = simplify(diagonal_weight * fluid_temperature)
    wall_equilibrium = simplify(diagonal_weight * wall_temperature)
    dirichlet_rhs = simplify(-outgoing_post + 2 * wall_equilibrium)
    adiabatic_rhs = outgoing_post

    closure_matrix = Matrix([[1], [1]])
    augmented = closure_matrix.row_join(Matrix([dirichlet_rhs, adiabatic_rhs]))
    # Applying Dirichlet then adiabatic leaves the latter value, and vice versa.
    dirichlet_then_adiabatic = adiabatic_rhs
    adiabatic_then_dirichlet = dirichlet_rhs
    order_difference = simplify(
        adiabatic_then_dirichlet - dirichlet_then_adiabatic
    )
    coefficient = simplify(
        order_difference
        / (mixed_temperature_curvature * grid_spacing**2)
    )
    heat = Symbol("Q_corner_direct")
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    source_increment = simplify(
        Rational(1, 36) * (1 - s_plus / 2) * heat
    )
    source_assignments = (
        ("shared_diagonal", "dirichlet_abb", source_increment),
        ("shared_diagonal", "adiabatic_bb", source_increment),
    )
    assignment_population_ids = tuple(
        population
        for population, _, source in source_assignments
        if source.has(heat)
    )
    shared_population_ids = tuple(dict.fromkeys(assignment_population_ids))
    return CornerStencilCheck(
        profile=profile,
        closure_rank=closure_matrix.rank(),
        augmented_rank=augmented.rank(),
        dirichlet_then_adiabatic=dirichlet_then_adiabatic,
        adiabatic_then_dirichlet=adiabatic_then_dirichlet,
        first_nonzero_coefficient=coefficient,
        simultaneously_satisfied=simplify(order_difference) == 0,
        naive_overwrite_heat_source_count=len(assignment_population_ids),
        shared_population_heat_source_count=len(shared_population_ids),
        diagonal_wall_distance=simplify(sqrt(2) * grid_spacing / 2),
        derivation=Tuple(
            order_difference,
            *(source for _, _, source in source_assignments),
        ),
        source_assignment_population_ids=assignment_population_ids,
        shared_source_population_ids=shared_population_ids,
    )


def _direct_quadratic_abb_stencil(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    pressure_ratio: Expr,
    homogeneous_flux_shift: Expr | None = None,
) -> tuple[Expr, Expr]:
    """Substitute a quadratic polynomial into the finite D1Q3 link chain."""

    x = Symbol("x")
    temperature_0, temperature_1 = symbols("T0_direct T1_direct")
    temperature_2 = Symbol("T2_direct", nonzero=True)
    heat = Symbol("Q_direct")
    temperature = temperature_0 + temperature_1 * x + temperature_2 * x**2
    gradient = temperature.diff(x)
    coefficients = symbols("direct_q0:9")
    g_zero = sum(coefficients[index] * x**index for index in range(3))
    g_plus = sum(coefficients[3 + index] * x**index for index in range(3))
    g_minus = sum(coefficients[6 + index] * x**index for index in range(3))
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    odd_shift = (
        sigma_g_minus
        if homogeneous_flux_shift is None
        else sympify(homogeneous_flux_shift)
    )
    s_minus = simplify(1 / (odd_shift + Rational(1, 2)))
    even = (g_plus + g_minus) / 2
    odd = (g_plus - g_minus) / 2
    equilibrium_zero = (Rational(2, 3) - pressure_ratio) * temperature
    equilibrium_normal = (
        Rational(1, 6) + pressure_ratio / 2
    ) * temperature
    odd_source = (
        Rational(1, 6)
        * (pressure_ratio + chi_kappa * Rational(1, 3))
        / Rational(1, 3)
        * gradient
        if homogeneous_flux_shift is None
        else sympify(0)
    )
    post_zero = (
        g_zero
        - s_plus * (g_zero - equilibrium_zero)
        + (1 - s_plus / 2) * Rational(2, 3) * heat
    )
    post_plus = (
        g_plus
        - s_plus * (even - equilibrium_normal)
        - s_minus * odd
        + (1 - s_plus / 2) * Rational(1, 6) * heat
        + (1 - s_minus / 2) * odd_source
    )
    post_minus = (
        g_minus
        - s_plus * (even - equilibrium_normal)
        + s_minus * odd
        + (1 - s_plus / 2) * Rational(1, 6) * heat
        - (1 - s_minus / 2) * odd_source
    )
    equations = _polynomial_coefficients(
        (
            g_zero - post_zero,
            g_plus.subs(x, x + 1) - post_plus,
            g_minus.subs(x, x - 1) - post_minus,
            g_zero + g_plus + g_minus - (temperature - heat / 2),
        ),
        (x,),
    )
    solution = solve(equations, (*coefficients, heat), dict=True, simplify=False)[0]
    fluid_midpoint = Rational(1, 2)
    wall_equilibrium = (
        Rational(1, 6) + pressure_ratio / 2
    ) * temperature_0
    link_residual = simplify(
        (
            g_plus.subs(x, fluid_midpoint)
            + post_minus.subs(x, fluid_midpoint)
            - 2 * wall_equilibrium
        ).subs(solution)
    )
    return simplify(link_residual / temperature_2), link_residual


def _direct_variable_pressure_abb_stencil(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    normal_temperature_gradient: Expr,
    normal_pressure_gradient: Expr,
) -> Expr:
    """Solve the affine T*pi link product with F_n/rho0=d_n pi."""

    x = Symbol("x")
    temperature_0, pressure_0 = symbols("T0_pressure pi0_pressure")
    temperature = temperature_0 + normal_temperature_gradient * x
    pressure = pressure_0 + normal_pressure_gradient * x
    coefficients = symbols("direct_p0:9")
    g_zero = sum(coefficients[index] * x**index for index in range(3))
    g_plus = sum(coefficients[3 + index] * x**index for index in range(3))
    g_minus = sum(coefficients[6 + index] * x**index for index in range(3))
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    s_minus = simplify(1 / (sigma_g_minus + Rational(1, 2)))
    even = (g_plus + g_minus) / 2
    odd = (g_plus - g_minus) / 2
    equilibrium_zero = (Rational(2, 3) - pressure) * temperature
    equilibrium_normal = (Rational(1, 6) + pressure / 2) * temperature
    product_gradient = simplify(
        pressure * normal_temperature_gradient
        + temperature * normal_pressure_gradient
    )
    odd_source = Rational(1, 6) * (
        product_gradient / Rational(1, 3)
        + chi_kappa * normal_temperature_gradient
    )
    post_zero = g_zero - s_plus * (g_zero - equilibrium_zero)
    post_plus = (
        g_plus
        - s_plus * (even - equilibrium_normal)
        - s_minus * odd
        + (1 - s_minus / 2) * odd_source
    )
    post_minus = (
        g_minus
        - s_plus * (even - equilibrium_normal)
        + s_minus * odd
        - (1 - s_minus / 2) * odd_source
    )
    equations = _polynomial_coefficients(
        (
            g_zero - post_zero,
            g_plus.subs(x, x + 1) - post_plus,
            g_minus.subs(x, x - 1) - post_minus,
            g_zero + g_plus + g_minus - temperature,
        ),
        (x,),
    )
    solution = solve(equations, coefficients, dict=True, simplify=False)[0]
    fluid_midpoint = Rational(1, 2)
    wall_equilibrium = (
        Rational(1, 6) + pressure_0 / 2
    ) * temperature_0
    return simplify(
        (
            g_plus.subs(x, fluid_midpoint)
            + post_minus.subs(x, fluid_midpoint)
            - 2 * wall_equilibrium
        ).subs(solution)
    )


def _direct_adiabatic_diagonal_stencil(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    chi_kappa: Expr,
    pressure_ratio: Expr,
) -> tuple[Expr, Expr]:
    """Sum the two diagonal reflected links for a tangential quadratic."""

    x = Symbol("tau")
    temperature_0, temperature_1 = symbols("T0_tau T1_tau")
    temperature_2 = Symbol("T2_tau", nonzero=True)
    temperature = temperature_0 + temperature_1 * x + temperature_2 * x**2
    gradient = temperature.diff(x)
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    s_minus = simplify(1 / (sigma_g_minus + Rational(1, 2)))
    residuals = []
    for direction in (-1, 1):
        coefficients = symbols(f"diag_{direction}_0:6")
        g_plus = sum(coefficients[index] * x**index for index in range(3))
        g_minus = sum(coefficients[3 + index] * x**index for index in range(3))
        even = (g_plus + g_minus) / 2
        odd = (g_plus - g_minus) / 2
        equilibrium = (
            Rational(1, 36) + pressure_ratio / 12
        ) * temperature
        odd_source = (
            Rational(1, 36)
            * direction
            * (pressure_ratio / Rational(1, 3) + chi_kappa)
            * gradient
        )
        post_plus = (
            g_plus
            - s_plus * (even - equilibrium)
            - s_minus * odd
            + (1 - s_minus / 2) * odd_source
        )
        post_minus = (
            g_minus
            - s_plus * (even - equilibrium)
            + s_minus * odd
            - (1 - s_minus / 2) * odd_source
        )
        equations = _polynomial_coefficients(
            (
                g_plus.subs(x, x + direction) - post_plus,
                g_minus.subs(x, x - direction) - post_minus,
            ),
            (x,),
        )
        solution = solve(equations, coefficients, dict=True, simplify=False)[0]
        residuals.append(simplify((g_plus - post_minus).subs(solution)))
    combined = simplify(sum(residuals))
    return simplify(combined / temperature_2), combined


def _direct_adiabatic_force_stencil(
    sigma_g_plus: Expr,
    sigma_g_minus: Expr,
    normal_force: Expr,
    normal_pressure_gradient: Expr,
) -> tuple[Expr, Expr, Expr]:
    """Solve a uniform-T normal pair with independent F_n and d_n pi."""

    x = Symbol("normal_coordinate")
    temperature = Symbol("T_uniform", nonzero=True)
    pressure_0 = Symbol("pi_wall")
    pressure = pressure_0 + normal_pressure_gradient * x
    coefficients = symbols("force_pair_0:6")
    g_zero = coefficients[0] + coefficients[1] * x
    g_plus = coefficients[2] + coefficients[3] * x
    g_minus = coefficients[4] + coefficients[5] * x
    s_plus = simplify(1 / (sigma_g_plus + Rational(1, 2)))
    s_minus = simplify(1 / (sigma_g_minus + Rational(1, 2)))
    even = (g_plus + g_minus) / 2
    odd = (g_plus - g_minus) / 2
    equilibrium_zero = (Rational(2, 3) - pressure) * temperature
    equilibrium_normal = (Rational(1, 6) + pressure / 2) * temperature
    odd_source = Rational(1, 6) * temperature * normal_force / Rational(1, 3)
    post_zero = g_zero - s_plus * (g_zero - equilibrium_zero)
    post_plus = (
        g_plus
        - s_plus * (even - equilibrium_normal)
        - s_minus * odd
        + (1 - s_minus / 2) * odd_source
    )
    post_minus = (
        g_minus
        - s_plus * (even - equilibrium_normal)
        + s_minus * odd
        - (1 - s_minus / 2) * odd_source
    )
    equations = _polynomial_coefficients(
        (
            g_zero - post_zero,
            g_plus.subs(x, x + 1) - post_plus,
            g_minus.subs(x, x - 1) - post_minus,
            g_zero + g_plus + g_minus - temperature,
        ),
        (x,),
    )
    solution = solve(equations, coefficients, dict=True, simplify=False)[0]
    midpoint = Rational(1, 2)
    reflected_flux = simplify(
        (
            g_plus.subs(x, midpoint) - post_minus.subs(x, midpoint)
        ).subs(solution)
    )
    return reflected_flux, temperature, normal_force


def _polynomial_coefficients(
    expressions: Sequence[Expr], variables: Sequence[Symbol]
) -> list[Expr]:
    coefficients: list[Expr] = []
    for expression in expressions:
        coefficients.extend(Poly(expand(expression), *variables).coeffs())
    return coefficients


def classify_magic(residual: BoundaryResidual) -> MagicClassification:
    """Return the declared class together with every retained coefficient row."""

    status = residual.recommended_status
    if status not in MAGIC_STATUSES:
        raise ValueError(f"unknown boundary classification: {status}")
    unsatisfied_jets = residual.unsatisfied_jets
    (
        rate_compatibility_status,
        rate_compatibility_conditions,
    ) = _classify_shear_bulk_rate_subsystem(residual)
    if status == "universal_magic":
        conditions = residual.solve_zero_conditions()
        if len(conditions) != 1:
            status = "no_single_magic"
        else:
            calibrated = residual.substitute(conditions)
            if calibrated.unsatisfied_jets:
                status = "no_single_magic"
            else:
                unsatisfied_jets = ()
    return MagicClassification(
        status=status,
        coefficients=residual.coefficients,
        assumptions=residual.assumptions,
        parameter_mapping=residual.parameter_mapping,
        unsatisfied_jets=unsatisfied_jets,
        cancelled_jets=residual.cancelled_jets,
        rate_compatibility_status=rate_compatibility_status,
        rate_compatibility_conditions=rate_compatibility_conditions,
    )


def _classify_shear_bulk_rate_subsystem(
    residual: BoundaryResidual,
) -> tuple[str | None, Mapping[Expr, Expr] | None]:
    """Classify only the resolved shear/bulk rows, never the full wall closure."""

    rate_rows = (
        residual.coefficients.get("shear_curvature"),
        residual.coefficients.get("bulk_curvature"),
    )
    if any(row is None for row in rate_rows) or not residual.calibration_parameters:
        return None, None
    parameter = residual.calibration_parameters[0]
    placeholder = Symbol("_rate_compatibility_product")
    equations = [simplify(row.subs(parameter, placeholder)) for row in rate_rows]
    roots = solve(equations, placeholder, dict=False)
    if isinstance(roots, dict):
        roots = [roots[placeholder]] if placeholder in roots else []
    elif roots and isinstance(roots[0], dict):
        roots = [root[placeholder] for root in roots if placeholder in root]
    if len(roots) == 1:
        return "restricted_calibration", {parameter: simplify(roots[0])}
    return "no_single_magic", {}
