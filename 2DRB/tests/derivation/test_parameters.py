"""Exact parameter-compatibility solver tests."""

from __future__ import annotations

import inspect
from pathlib import Path
import unittest
from unittest.mock import patch

from sympy import Rational, simplify, sqrt, symbols

from tools.derivation.boundary import BoundaryResidual, MagicClassification

from tools.derivation.d2q9_temperature import (
    CanonicalQuarticCondition,
    D2Q9EquivalentCoefficients,
    QuarticConditionSystem,
    amplification_route,
    canonical_quartic_condition,
    quartic_condition_system,
)
import tools.derivation.parameters as parameters_module


try:
    from tools.derivation.parameters import (
        ParameterReport,
        derive_scalar_compatibility,
        recover_baseline_quartic_family,
        solve_flow_parameters,
        solve_scalar_parameters,
    )
except (ImportError, ModuleNotFoundError) as exc:
    PARAMETERS_IMPORT_ERROR = exc
else:
    PARAMETERS_IMPORT_ERROR = None


class BaselineQuarticFamilyTests(unittest.TestCase):
    """Recover the full Task 5 family instead of selecting a known point."""

    def require_parameters(self) -> None:
        self.assertIsNone(
            PARAMETERS_IMPORT_ERROR,
            f"parameter solver is unavailable: {PARAMETERS_IMPORT_ERROR}",
        )

    def test_baseline_family_is_generated_for_symbolic_nominal_odd_shift(self):
        self.require_parameters()
        sigma_odd = symbols("sigma_o", positive=True)
        try:
            report = recover_baseline_quartic_family(sigma_odd=sigma_odd)
        except Exception as exc:
            self.fail(f"baseline family generation failed: {exc}")
        self.assertIsInstance(report, ParameterReport)
        self.assertEqual(report.status, "feasible_exact")
        self.assertEqual(
            simplify(
                report.exact_substitutions["sigma_even"]
                - (4 * sigma_odd**2 + 1) / (8 * sigma_odd)
            ),
            0,
        )
        self.assertNotEqual(
            simplify(report.exact_substitutions["sigma_even"] - 1 / 3),
            0,
            "the solver must return the family, not one memorized TRT point",
        )

    def test_baseline_calls_and_retains_actual_task5_route_and_condition(self):
        self.require_parameters()
        sigma_odd = symbols("sigma_o", positive=True)
        with patch(
            "tools.derivation.parameters.amplification_route",
            wraps=amplification_route,
        ) as route_call, patch(
            "tools.derivation.parameters.quartic_condition_system",
            wraps=quartic_condition_system,
        ) as condition_call:
            report = recover_baseline_quartic_family(sigma_odd=sigma_odd)
        route_call.assert_called_once()
        condition_call.assert_called_once()
        evidence = report.consumed_evidence
        self.assertIsInstance(evidence.bulk_route, D2Q9EquivalentCoefficients)
        self.assertIsInstance(evidence.bulk_conditions, QuarticConditionSystem)
        sigma_even = report.exact_substitutions["sigma_even"]
        self.assertEqual(simplify(evidence.bulk_route.pde_c40.subs(
            symbols("sigma_e"), sigma_even
        )), 0)
        self.assertEqual(simplify(evidence.bulk_route.pde_c22.subs(
            symbols("sigma_e"), sigma_even
        )), 0)

    def test_rate_checks_are_separate_for_nominal_and_physical_blocks(self):
        self.require_parameters()
        report = recover_baseline_quartic_family(sigma_odd=Rational(1, 4))
        self.assertEqual(
            set(report.open_interval_checks),
            {
                "nominal_odd_rate",
                "nominal_even_rate",
                "physical_flux_rate",
            },
        )
        self.assertTrue(all(report.open_interval_checks.values()))
        self.assertEqual(
            report.collision_rates["nominal_odd"], Rational(4, 3)
        )
        self.assertEqual(
            report.collision_rates["physical_flux"], Rational(4, 3)
        )


class ScalarCompatibilityTests(unittest.TestCase):
    """Combine Task 5 quartic families with Task 6 restricted ABB only."""

    def require_parameters(self) -> None:
        self.assertIsNone(
            PARAMETERS_IMPORT_ERROR,
            f"parameter solver is unavailable: {PARAMETERS_IMPORT_ERROR}",
        )

    def test_external_and_feedback_identities_are_generated_by_exact_elimination(self):
        self.require_parameters()
        pi, chi = symbols("pi chi_kappa")
        a = Rational(1, 3) + pi
        b = (1 - chi) * Rational(1, 3)
        try:
            external = derive_scalar_compatibility(
                case="external",
                pressure_ratio=pi,
                chi_kappa=chi,
                cs2=Rational(1, 3),
            )
            feedback = derive_scalar_compatibility(
                case="feedback",
                pressure_ratio=pi,
                chi_kappa=chi,
                cs2=Rational(1, 3),
            )
        except Exception as exc:
            self.fail(f"scalar exact elimination is unavailable: {exc}")
        self.assertEqual(
            simplify(
                feedback.exact_substitutions["required_transport_ratio_squared"]
                - a * (3 * a + 1) / 48
            ),
            0,
        )
        self.assertEqual(
            simplify(
                external.exact_substitutions["required_transport_ratio_squared"]
                - (12 * a * b + 5 * a - 4 * b - 9 * a**2) / 48
            ),
            0,
        )
        for report in (external, feedback):
            self.assertTrue(
                getattr(report, "is_conditional", False),
                "a symbolic compatibility family must not be reported unconditional",
            )
            self.assertTrue(
                getattr(report, "feasibility_conditions", None),
                "a symbolic compatibility family must expose its physical conditions",
            )
            self.assertIn("a_positive", report.feasibility_conditions)
            self.assertIn("b_positive", report.feasibility_conditions)
            self.assertFalse(
                any(
                    name.endswith("_not_positive")
                    for name in report.feasibility_conditions
                ),
                "feasibility conditions must name the predicate being required",
            )
            self.assertIsInstance(
                report.consumed_evidence.boundary_residual, BoundaryResidual
            )
            self.assertIsInstance(
                report.consumed_evidence.boundary_classification,
                MagicClassification,
            )
            self.assertEqual(
                report.consumed_evidence.boundary_classification.status,
                "restricted_calibration",
            )

    def test_pi_zero_positive_compatible_point_fixes_physical_and_even_shifts(self):
        self.require_parameters()
        chi = Rational(1, 4)
        report = solve_scalar_parameters(
            case="feedback",
            pressure_ratio=0,
            chi_kappa=chi,
            cs2=Rational(1, 3),
            kappa=1 / sqrt(72),
            dt=1,
            require_bulk_quartic=True,
            boundary_policy="restricted_abb",
        )
        self.assertEqual(report.status, "feasible_restricted")
        self.assertEqual(
            simplify(report.exact_substitutions["sigma_flux"] - sqrt(2) / 4),
            0,
        )
        self.assertEqual(
            simplify(report.exact_substitutions["sigma_even"] - 3 * sqrt(2) / 8),
            0,
        )
        self.assertEqual(
            simplify(
                report.exact_substitutions["sigma_odd"]
                - 1 / (2 * sqrt(2) * (1 - chi))
            ),
            0,
        )
        self.assertTrue(all(report.open_interval_checks.values()))
        self.assertIn(
            "flat_grid_aligned_halfway_wall", report.assumptions
        )

    def test_general_dt_scales_the_compatible_diffusivity(self):
        self.require_parameters()
        report = solve_scalar_parameters(
            case="external",
            pressure_ratio=0,
            chi_kappa=Rational(1, 4),
            cs2=Rational(1, 3),
            kappa=3 / sqrt(72),
            dt=3,
            require_bulk_quartic=True,
            boundary_policy="restricted_abb",
        )
        self.assertEqual(report.status, "feasible_restricted")
        self.assertEqual(
            report.exact_substitutions["required_kappa"], 3 / sqrt(72)
        )

    def test_low_positive_kappa_is_rejected_when_both_constraints_are_mandatory(self):
        self.require_parameters()
        report = solve_scalar_parameters(
            case="feedback",
            pressure_ratio=0,
            chi_kappa=Rational(1, 4),
            cs2=Rational(1, 3),
            kappa=Rational(1, 1000),
            dt=1,
            require_bulk_quartic=True,
            boundary_policy="restricted_abb",
        )
        self.assertEqual(report.status, "no_feasible_solution")
        self.assertIn("bulk_quartic_and_restricted_abb", report.violated_constraints)
        self.assertEqual(
            report.minimal_extension,
            "explicit_restricted_abb_boundary_correction",
        )
        self.assertNotIn("split_even", report.minimal_extension)


class CanonicalEvidenceConsumptionTests(unittest.TestCase):
    """The parameter layer consumes Task 5 canonical conditions without copies."""

    def require_parameters(self) -> None:
        self.assertIsNone(
            PARAMETERS_IMPORT_ERROR,
            f"parameter solver is unavailable: {PARAMETERS_IMPORT_ERROR}",
        )

    def test_parameter_source_contains_no_private_copy_of_task5_bulk_family(self):
        self.require_parameters()
        source = inspect.getsource(parameters_module)
        self.assertNotIn("_reviewed_task5_bulk_family", source)

    def test_derive_retains_the_exact_task5_canonical_objects_it_consumes(self):
        self.require_parameters()
        parameter_api = getattr(
            parameters_module,
            "canonical_quartic_condition",
            None,
        )
        self.assertTrue(
            callable(parameter_api),
            "parameters must import and consume Task 5 canonical conditions",
        )
        captured: list[CanonicalQuarticCondition] = []

        def capture(**kwargs):
            result = canonical_quartic_condition(**kwargs)
            captured.append(result)
            return result

        with patch.object(
            parameters_module,
            "canonical_quartic_condition",
            side_effect=capture,
        ) as canonical_call:
            report = derive_scalar_compatibility(
                case="feedback",
                pressure_ratio=Rational(2, 3),
                chi_kappa=Rational(1, 4),
            )
        canonical_call.assert_called_once()
        self.assertEqual(len(captured), 1)
        self.assertIs(report.consumed_evidence.bulk_route, captured[0].coefficients)
        self.assertIs(report.consumed_evidence.bulk_conditions, captured[0].conditions)
        self.assertIsNone(report.consumed_evidence.bulk_specialization)
        self.assertEqual(
            report.exact_substitutions["required_transport_ratio_squared"],
            Rational(1, 12),
        )


class ScalarCorrectionAndDegeneracyTests(unittest.TestCase):
    """Keep low-diffusivity correction and singular branches explicit."""

    def require_parameters(self) -> None:
        self.assertIsNone(
            PARAMETERS_IMPORT_ERROR,
            f"parameter solver is unavailable: {PARAMETERS_IMPORT_ERROR}",
        )

    def test_low_kappa_bulk_quartic_with_explicit_wall_correction_is_feasible(self):
        self.require_parameters()
        try:
            report = solve_scalar_parameters(
                case="feedback",
                pressure_ratio=0,
                chi_kappa=Rational(1, 4),
                cs2=Rational(1, 3),
                kappa=Rational(1, 1000),
                dt=1,
                require_bulk_quartic=True,
                boundary_policy="explicit_correction",
            )
        except Exception as exc:
            self.fail(f"explicit-correction branch is unavailable: {exc}")
        self.assertEqual(report.status, "boundary_correction_required")
        self.assertEqual(report.exact_substitutions["sigma_odd"], Rational(1, 250))
        self.assertEqual(
            report.exact_substitutions["sigma_even"], Rational(250009, 6000)
        )
        self.assertEqual(
            report.exact_substitutions["sigma_flux"], Rational(3, 1000)
        )
        self.assertEqual(report.collision_rates["nominal_odd"], Rational(125, 63))
        self.assertEqual(
            report.collision_rates["nominal_even"], Rational(6000, 253009)
        )
        self.assertEqual(
            report.collision_rates["physical_flux"], Rational(1000, 503)
        )
        self.assertTrue(all(report.open_interval_checks.values()))
        self.assertIn("restricted_abb_product", report.violated_constraints)
        self.assertEqual(
            report.minimal_extension,
            "explicit_restricted_abb_boundary_correction",
        )

    def test_a_one_special_branches_are_derived_without_dividing_by_a_minus_one(self):
        self.require_parameters()
        try:
            feedback = derive_scalar_compatibility(
                case="feedback",
                pressure_ratio=Rational(2, 3),
                chi_kappa=Rational(1, 4),
            )
            external = derive_scalar_compatibility(
                case="external",
                pressure_ratio=Rational(2, 3),
                chi_kappa=-1,
            )
        except Exception as exc:
            self.fail(f"a=1 special branches are unavailable: {exc}")
        self.assertEqual(
            feedback.exact_substitutions["required_transport_ratio_squared"],
            Rational(1, 12),
        )
        self.assertEqual(
            external.exact_substitutions["required_transport_ratio_squared"],
            Rational(1, 36),
        )
        self.assertIn("a_equals_one_direct_condition", feedback.assumptions)
        self.assertIn("a_equals_one_direct_condition", external.assumptions)

    def test_a_zero_and_b_zero_are_infeasible_not_abb_limits(self):
        self.require_parameters()
        try:
            a_zero = derive_scalar_compatibility(
                case="external",
                pressure_ratio=Rational(-1, 3),
                chi_kappa=Rational(1, 4),
            )
            b_zero = derive_scalar_compatibility(
                case="feedback",
                pressure_ratio=0,
                chi_kappa=1,
            )
        except Exception as exc:
            self.fail(f"degenerate scalar branches raised instead of reporting: {exc}")
        self.assertEqual(a_zero.status, "no_feasible_solution")
        self.assertIn("a_not_positive", a_zero.violated_constraints)
        self.assertEqual(b_zero.status, "no_feasible_solution")
        self.assertIn("b_not_positive", b_zero.violated_constraints)

    def test_negative_b_is_rejected_by_derive_and_target_solver(self):
        self.require_parameters()
        try:
            family = derive_scalar_compatibility(
                case="feedback",
                pressure_ratio=0,
                chi_kappa=2,
            )
            target = solve_scalar_parameters(
                case="feedback",
                pressure_ratio=0,
                chi_kappa=2,
                cs2=Rational(1, 3),
                kappa=1 / sqrt(72),
                dt=1,
                require_bulk_quartic=True,
                boundary_policy="restricted_abb",
            )
        except Exception as exc:
            self.fail(f"negative-b inputs raised instead of reporting: {exc}")
        for report in (family, target):
            self.assertEqual(report.status, "no_feasible_solution")
            self.assertIn("b_not_positive", report.violated_constraints)

    def test_zero_target_transport_is_rejected_before_rate_generation(self):
        self.require_parameters()
        try:
            report = solve_scalar_parameters(
                case="external",
                pressure_ratio=0,
                chi_kappa=Rational(1, 4),
                cs2=Rational(1, 3),
                kappa=0,
                dt=1,
                require_bulk_quartic=True,
                boundary_policy="restricted_abb",
            )
        except Exception as exc:
            self.fail(f"zero transport raised instead of reporting: {exc}")
        self.assertEqual(report.status, "no_feasible_solution")
        self.assertIn("transport_ratio_not_positive", report.violated_constraints)

    def test_negative_external_radicand_returns_no_feasible_solution(self):
        self.require_parameters()
        try:
            report = derive_scalar_compatibility(
                case="external",
                pressure_ratio=Rational(7, 15),
                chi_kappa=Rational(1, 4),
            )
        except Exception as exc:
            self.fail(f"negative-radicand branch raised instead of reporting: {exc}")
        self.assertEqual(report.status, "no_feasible_solution")
        self.assertIn("negative_compatibility_radicand", report.violated_constraints)
        self.assertTrue(report.minimal_extension)

    def test_feedback_closure_singularity_is_rejected_before_bulk_route(self):
        self.require_parameters()
        try:
            report = solve_scalar_parameters(
                case="feedback",
                pressure_ratio=Rational(-2, 9),
                chi_kappa=Rational(1, 4),
                cs2=Rational(1, 3),
                kappa=Rational(-1, 18),
                dt=1,
                require_bulk_quartic=True,
                boundary_policy="restricted_abb",
            )
        except Exception as exc:
            self.fail(f"feedback closure singularity raised instead of reporting: {exc}")
        self.assertEqual(report.status, "no_feasible_solution")
        self.assertIn("feedback_closure_singular", report.violated_constraints)
        self.assertIn("transport_ratio_not_positive", report.violated_constraints)

    def test_a_one_split_even_extension_is_explicitly_only_a_candidate(self):
        self.require_parameters()
        report = solve_scalar_parameters(
            case="feedback",
            pressure_ratio=Rational(2, 3),
            chi_kappa=Rational(1, 4),
            cs2=Rational(1, 3),
            kappa=1 / sqrt(12),
            dt=1,
            require_bulk_quartic=True,
            boundary_policy="explicit_correction",
        )
        self.assertEqual(report.status, "mrt_extension_required")
        self.assertEqual(
            report.minimal_extension,
            "split_even_mrt_candidate_requiring_mode_jacobian_derivation",
        )

    def test_target_solver_requires_a_specified_chi_kappa(self):
        self.require_parameters()
        chi = symbols("chi_kappa")
        try:
            report = solve_scalar_parameters(
                case="feedback",
                pressure_ratio=0,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                kappa=1 / sqrt(72),
                dt=1,
                require_bulk_quartic=True,
                boundary_policy="restricted_abb",
            )
        except Exception as exc:
            self.fail(f"unspecified chi_kappa raised instead of reporting: {exc}")
        self.assertEqual(report.status, "no_feasible_solution")
        self.assertIn("chi_kappa_must_be_specified", report.violated_constraints)


class FlowCompatibilityTests(unittest.TestCase):
    """Consume the Task 6 shear product without promoting it to a universal wall."""

    def require_parameters(self) -> None:
        self.assertIsNone(
            PARAMETERS_IMPORT_ERROR,
            f"parameter solver is unavailable: {PARAMETERS_IMPORT_ERROR}",
        )

    def test_restricted_shear_product_generates_nominal_odd_rate_and_limit(self):
        self.require_parameters()
        try:
            report = solve_flow_parameters(
                nu=Rational(1, 1000),
                cs2=Rational(1, 3),
                dt=1,
                chi_s=Rational(1, 4),
                chi_b=Rational(1, 5),
                retain_trace_jets=False,
            )
        except Exception as exc:
            self.fail(f"restricted flow compatibility is unavailable: {exc}")
        self.assertEqual(report.status, "feasible_restricted")
        self.assertEqual(
            report.exact_substitutions["sigma_f_plus"], Rational(1, 250)
        )
        self.assertEqual(
            report.exact_substitutions["sigma_f_minus"], Rational(125, 2)
        )
        self.assertEqual(report.collision_rates["nominal_odd"], Rational(1, 63))
        self.assertEqual(
            report.exact_substitutions["s_f_minus_limit_nu_to_zero_positive"],
            0,
        )
        self.assertTrue(all(report.open_interval_checks.values()))
        evidence = report.consumed_evidence
        self.assertIsInstance(evidence.boundary_residual, BoundaryResidual)
        self.assertIsInstance(evidence.boundary_classification, MagicClassification)
        self.assertEqual(evidence.boundary_classification.status, "restricted_calibration")

    def test_unequal_shear_and_bulk_feedback_rejects_one_magic_with_trace_jets(self):
        self.require_parameters()
        try:
            report = solve_flow_parameters(
                nu=Rational(1, 1000),
                cs2=Rational(1, 3),
                dt=1,
                chi_s=Rational(1, 4),
                chi_b=Rational(1, 5),
                retain_trace_jets=True,
            )
        except Exception as exc:
            self.fail(f"general flow classification is unavailable: {exc}")
        self.assertEqual(report.status, "mrt_extension_required")
        classification = report.consumed_evidence.boundary_classification
        self.assertIsInstance(classification, MagicClassification)
        self.assertEqual(classification.rate_compatibility_status, "no_single_magic")
        self.assertIn("no_single_shear_bulk_magic", report.violated_constraints)
        self.assertEqual(
            report.minimal_extension,
            "separate_shear_bulk_relaxation_and_general_velocity_boundary_correction",
        )
        self.assertNotIn("_or_", report.minimal_extension)

    def test_equal_feedback_still_preserves_general_unresolved_wall_classification(self):
        self.require_parameters()
        try:
            report = solve_flow_parameters(
                nu=Rational(1, 1000),
                cs2=Rational(1, 3),
                dt=1,
                chi_s=Rational(1, 4),
                chi_b=Rational(1, 4),
                retain_trace_jets=True,
            )
        except Exception as exc:
            self.fail(f"general equal-feedback classification is unavailable: {exc}")
        classification = report.consumed_evidence.boundary_classification
        self.assertEqual(report.status, "boundary_correction_required")
        self.assertEqual(classification.status, "boundary_correction_required")
        self.assertEqual(
            classification.rate_compatibility_status, "restricted_calibration"
        )
        self.assertIn("general_velocity_wall_jets", report.violated_constraints)

    def test_negative_nu_precedes_general_trace_classification(self):
        self.require_parameters()
        try:
            report = solve_flow_parameters(
                nu=Rational(-1, 1000),
                cs2=Rational(1, 3),
                dt=1,
                chi_s=Rational(1, 4),
                chi_b=Rational(1, 5),
                retain_trace_jets=True,
            )
        except Exception as exc:
            self.fail(f"negative nu raised instead of reporting: {exc}")
        self.assertEqual(report.status, "no_feasible_solution")
        self.assertIn("nu_not_positive", report.violated_constraints)
        self.assertIn("nominal_even_rate", report.violated_constraints)

    def test_out_of_interval_actual_flow_rate_precedes_wall_status(self):
        self.require_parameters()
        try:
            report = solve_flow_parameters(
                nu=Rational(1, 1000),
                cs2=Rational(1, 3),
                dt=1,
                chi_s=2,
                chi_b=Rational(1, 5),
                retain_trace_jets=True,
            )
        except Exception as exc:
            self.fail(f"out-of-range flow rate raised instead of reporting: {exc}")
        self.assertEqual(report.status, "no_feasible_solution")
        self.assertIn("one_minus_chi_s_not_positive", report.violated_constraints)
        self.assertIn("nominal_even_rate", report.violated_constraints)

    def test_negative_bulk_transport_precedes_restricted_or_general_wall_status(self):
        self.require_parameters()
        for retain_trace_jets in (False, True):
            with self.subTest(retain_trace_jets=retain_trace_jets):
                report = solve_flow_parameters(
                    nu=Rational(1, 1000),
                    cs2=Rational(1, 3),
                    dt=1,
                    chi_s=Rational(1, 4),
                    chi_b=2,
                    retain_trace_jets=retain_trace_jets,
                )
                self.assertEqual(report.status, "no_feasible_solution")
                self.assertEqual(
                    report.exact_substitutions["sigma_bulk_physical"],
                    Rational(-1, 250),
                )
                self.assertEqual(
                    report.exact_substitutions["nu_bulk_2d"],
                    Rational(-1, 750),
                )
                self.assertEqual(
                    report.collision_rates["physical_bulk"], Rational(125, 62)
                )
                self.assertFalse(
                    report.open_interval_checks["physical_bulk_rate"]
                )
                self.assertIn(
                    "one_minus_chi_b_not_positive", report.violated_constraints
                )
                self.assertIn(
                    "sigma_bulk_physical_not_positive",
                    report.violated_constraints,
                )
                self.assertIn("nu_bulk_2d_not_positive", report.violated_constraints)
                self.assertIn("physical_bulk_rate", report.violated_constraints)


class AlgorithmDocumentationTests(unittest.TestCase):
    def test_flow_feedback_affine_term_uses_latex_thin_space_not_comma(self):
        root = Path(__file__).resolve().parents[2]
        algorithm = (root / "docs/derivation/chapters/08-algorithm.md").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\Delta t\,a_{uF}", algorithm)
        self.assertNotIn(r"\Delta t,a_{uF}", algorithm)


if __name__ == "__main__":
    unittest.main()
