"""Exact source-aware boundary residual tests."""

from __future__ import annotations

import unittest
from unittest.mock import patch

from sympy import Rational, simplify, sqrt, symbols


try:
    import tools.derivation.boundary as boundary_module
    from tools.derivation.boundary import (
        BoundaryResidual,
        adiabatic_bb_residual,
        classify_magic,
        manufactured_corner_stencil,
        manufactured_temperature_stencil,
        manufactured_velocity_stencil,
        mixed_corner_residual,
        temperature_abb_residual,
        velocity_bb_residual,
    )
except (ImportError, ModuleNotFoundError) as exc:
    BOUNDARY_IMPORT_ERROR = exc
else:
    BOUNDARY_IMPORT_ERROR = None


class ClassicalVelocityBoundaryTests(unittest.TestCase):
    """Lock the paper normalization before source-feedback extensions."""

    def require_boundary_module(self) -> None:
        self.assertIsNone(
            BOUNDARY_IMPORT_ERROR,
            f"boundary derivation is unavailable: {BOUNDARY_IMPORT_ERROR}",
        )

    def test_multireflection_quarter_maps_exactly_to_three_sixteenths(self):
        self.require_boundary_module()
        sigma_even, sigma_odd = symbols("sigma_f_plus sigma_f_minus")
        residual = velocity_bb_residual(
            sigma_f_plus=sigma_even,
            sigma_f_minus=sigma_odd,
            chi_s=0,
            force_mode="uniform_body_force",
            normalization="glt2003",
            momentum_gauge="half_force",
            geometry="flat_grid_aligned_halfway",
        )
        self.assertEqual(
            simplify(
                residual.parameter_mapping["Lambda_squared"]
                - Rational(4, 3) * sigma_even * sigma_odd
            ),
            0,
        )
        self.assertEqual(
            simplify(
                residual.parameter_mapping["lambda_nu"]
                - 1 / (sigma_even + Rational(1, 2))
            ),
            0,
        )
        self.assertEqual(
            simplify(
                residual.parameter_mapping["lambda_2"]
                - 1 / (sigma_odd + Rational(1, 2))
            ),
            0,
        )
        self.assertEqual(residual.parameter_mapping["henon_even"], sigma_even)
        self.assertEqual(residual.parameter_mapping["henon_odd"], sigma_odd)
        condition = residual.solve_zero_conditions()
        self.assertEqual(condition[sigma_even * sigma_odd], Rational(3, 16))
        self.assertEqual(classify_magic(residual).status, "restricted_calibration")

    def test_quarter_and_three_sixteenths_are_not_naked_constants(self):
        self.require_boundary_module()
        sigma_even, sigma_odd = symbols("sigma_f_plus sigma_f_minus")
        residual = velocity_bb_residual(
            sigma_f_plus=sigma_even,
            sigma_f_minus=sigma_odd,
            chi_s=0,
            force_mode="uniform_body_force",
            normalization="glt2003",
            momentum_gauge="half_force",
            geometry="flat_grid_aligned_halfway",
        )
        assumptions = set(residual.assumptions)
        self.assertIn("steady_linear_stokes", assumptions)
        self.assertIn("uniform_body_force", assumptions)
        self.assertIn("half_force_momentum_gauge", assumptions)
        self.assertIn("flat_grid_aligned_halfway_wall", assumptions)
        self.assertIn("no_lbm_cde_feedback", assumptions)

    def test_pressure_boundary_drive_is_a_distinct_three_eighths_case(self):
        self.require_boundary_module()
        sigma_even, sigma_odd = symbols("sigma_f_plus sigma_f_minus")
        residual = velocity_bb_residual(
            sigma_f_plus=sigma_even,
            sigma_f_minus=sigma_odd,
            chi_s=0,
            force_mode="pressure_boundary",
            normalization="wang_luo_henon_product",
            momentum_gauge="pressure_boundary",
            geometry="flat_grid_aligned_halfway",
        )
        condition = residual.solve_zero_conditions()
        self.assertEqual(condition[sigma_even * sigma_odd], Rational(3, 8))
        self.assertIn("pressure_boundary_drive", residual.assumptions)
        self.assertNotIn("uniform_body_force", residual.assumptions)
        self.assertEqual(classify_magic(residual).status, "restricted_calibration")

    def test_linear_couette_checks_halfway_sign_but_not_quadratic_magic(self):
        self.require_boundary_module()
        check = manufactured_velocity_stencil(
            profile="linear_couette",
            geometry="flat_grid_aligned_halfway",
        )
        self.assertEqual(simplify(check.first_nonzero_coefficient), 0)
        self.assertFalse(check.probes_quadratic_slip)

    def test_hydrostatic_half_force_reconstruction_has_no_spurious_slip(self):
        self.require_boundary_module()
        check = manufactured_velocity_stencil(
            profile="hydrostatic_rest",
            geometry="flat_grid_aligned_halfway",
        )
        self.assertEqual(simplify(check.first_nonzero_coefficient), 0)
        self.assertTrue(check.uses_half_force_reconstruction)

    def test_poiseuille_calibrations_are_source_evidence_not_independent_probes(self):
        self.require_boundary_module()
        sigma_even, sigma_odd, force, pressure_gradient = symbols(
            "sigma_f_plus sigma_f_minus F_drive dpdx"
        )
        with patch(
            "tools.derivation.boundary.velocity_bb_residual",
            side_effect=AssertionError("manufactured route called residual solver"),
        ):
            try:
                forced = manufactured_velocity_stencil(
                    profile="uniform_force_quadratic_poiseuille",
                    geometry="flat_grid_aligned_halfway",
                    sigma_f_plus=sigma_even,
                    sigma_f_minus=sigma_odd,
                    force_amplitude=force,
                    normalization="glt2003",
                    momentum_gauge="half_force",
                )
                pressured = manufactured_velocity_stencil(
                    profile="pressure_driven_quadratic_poiseuille",
                    geometry="flat_grid_aligned_halfway",
                    sigma_f_plus=sigma_even,
                    sigma_f_minus=sigma_odd,
                    pressure_gradient=pressure_gradient,
                    normalization="wang_luo_henon_product",
                    momentum_gauge="pressure_boundary",
                )
            except Exception as exc:
                self.fail(
                    f"Poiseuille source-evidence records are unavailable: {exc}"
                )
        self.assertEqual(
            simplify(
                forced.first_nonzero_coefficient
                - (Rational(16, 3) * sigma_even * sigma_odd - 1)
            ),
            0,
        )
        self.assertEqual(
            simplify(
                pressured.first_nonzero_coefficient
                - (sigma_even * sigma_odd - Rational(3, 8))
            ),
            0,
        )
        self.assertEqual(forced.verification_status, "source_evidence_only")
        self.assertEqual(pressured.verification_status, "source_evidence_only")
        self.assertFalse(forced.probes_quadratic_slip)
        self.assertFalse(pressured.probes_quadratic_slip)


class SourceFeedbackVelocityBoundaryTests(unittest.TestCase):
    """Keep physical shear/bulk shifts distinct from nominal TRT ghosts."""

    def call_residual(self, **kwargs):
        self.assertIsNone(
            BOUNDARY_IMPORT_ERROR,
            f"boundary derivation is unavailable: {BOUNDARY_IMPORT_ERROR}",
        )
        try:
            return velocity_bb_residual(**kwargs)
        except Exception as exc:  # RED converts a missing branch to assertion FAIL.
            self.fail(f"source-feedback velocity residual is unavailable: {exc}")

    def test_one_dimensional_feedback_candidate_uses_physical_shear_only(self):
        sigma_even, sigma_odd, chi_s = symbols(
            "sigma_f_plus sigma_f_minus chi_s"
        )
        residual = self.call_residual(
            sigma_f_plus=sigma_even,
            sigma_f_minus=sigma_odd,
            chi_s=chi_s,
            chi_b=symbols("chi_b"),
            force_mode="uniform_body_force",
            normalization="source_feedback_henon_product",
            momentum_gauge="half_force",
            geometry="flat_grid_aligned_halfway",
        )
        candidate = (1 - chi_s) * sigma_even * sigma_odd
        condition = residual.solve_zero_conditions()
        self.assertEqual(condition[candidate], Rational(3, 16))
        self.assertEqual(classify_magic(residual).status, "restricted_calibration")
        self.assertIn("one_dimensional_shear_only", residual.assumptions)

    def test_feedback_mapping_does_not_replace_nominal_ghost_rates(self):
        sigma_even, sigma_odd, chi_s, chi_b = symbols(
            "sigma_f_plus sigma_f_minus chi_s chi_b"
        )
        residual = self.call_residual(
            sigma_f_plus=sigma_even,
            sigma_f_minus=sigma_odd,
            chi_s=chi_s,
            chi_b=chi_b,
            force_mode="general_source_aware",
            normalization="source_feedback_henon_product",
            momentum_gauge="half_force",
            geometry="flat_grid_aligned_halfway",
        )
        mapping = residual.parameter_mapping
        self.assertEqual(mapping["flow_shear_physical"], (1 - chi_s) * sigma_even)
        self.assertEqual(mapping["flow_bulk_physical"], (1 - chi_b) * sigma_even)
        self.assertEqual(mapping["flow_even_ghost_nominal"], sigma_even)
        self.assertEqual(mapping["flow_odd_ghost_nominal"], sigma_odd)

    def test_distinct_shear_and_bulk_feedback_has_no_single_magic(self):
        sigma_even, sigma_odd = symbols("sigma_f_plus sigma_f_minus")
        residual = self.call_residual(
            sigma_f_plus=sigma_even,
            sigma_f_minus=sigma_odd,
            chi_s=Rational(1, 4),
            chi_b=Rational(1, 2),
            force_mode="general_source_aware",
            normalization="source_feedback_henon_product",
            momentum_gauge="half_force",
            geometry="flat_grid_aligned_halfway",
        )
        self.assertEqual(residual.solve_zero_conditions(), {})
        classification = classify_magic(residual)
        self.assertEqual(classification.status, "boundary_correction_required")
        self.assertEqual(
            getattr(classification, "rate_compatibility_status", None),
            "no_single_magic",
        )
        self.assertEqual(
            getattr(classification, "rate_compatibility_conditions", None),
            {},
        )
        self.assertTrue(
            {"shear_curvature", "bulk_curvature"}.issubset(
                classification.unsatisfied_jets
            )
        )

    def test_equal_shear_and_bulk_feedback_has_one_restricted_rate_calibration(self):
        sigma_even, sigma_odd = symbols("sigma_f_plus sigma_f_minus")
        residual = self.call_residual(
            sigma_f_plus=sigma_even,
            sigma_f_minus=sigma_odd,
            chi_s=Rational(1, 4),
            chi_b=Rational(1, 4),
            force_mode="general_source_aware",
            normalization="source_feedback_henon_product",
            momentum_gauge="half_force",
            geometry="flat_grid_aligned_halfway",
        )
        classification = classify_magic(residual)
        self.assertEqual(classification.status, "boundary_correction_required")
        self.assertEqual(
            getattr(classification, "rate_compatibility_status", None),
            "restricted_calibration",
        )
        self.assertEqual(
            getattr(classification, "rate_compatibility_conditions", None),
            {sigma_even * sigma_odd: Rational(1, 4)},
        )
        self.assertTrue(classification.unsatisfied_jets)

    def test_general_velocity_table_retains_independent_unresolved_wall_jets(self):
        sigma_even, sigma_odd, chi_s, chi_b = symbols(
            "sigma_f_plus sigma_f_minus chi_s chi_b"
        )
        residual = self.call_residual(
            sigma_f_plus=sigma_even,
            sigma_f_minus=sigma_odd,
            chi_s=chi_s,
            chi_b=chi_b,
            force_mode="general_source_aware",
            normalization="source_feedback_henon_product",
            momentum_gauge="half_force",
            geometry="flat_grid_aligned_halfway",
        )
        required = {
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
        }
        self.assertTrue(required.issubset(residual.coefficients))
        self.assertTrue(
            all(residual.coefficients[name] != 0 for name in required),
            "unknown wall-jet coefficients must be explicit unresolved symbols",
        )
        self.assertIn("unresolved_full_velocity_wall_chain", residual.assumptions)
        self.assertEqual(
            classify_magic(residual).status, "boundary_correction_required"
        )


class TemperatureAndAdiabaticBoundaryTests(unittest.TestCase):
    """Audit the scalar population chain without hiding independent jets."""

    def require_boundary_module(self) -> None:
        self.assertIsNone(
            BOUNDARY_IMPORT_ERROR,
            f"scalar boundary derivation is unavailable: {BOUNDARY_IMPORT_ERROR}",
        )

    def call_temperature(self, **kwargs):
        self.require_boundary_module()
        try:
            return temperature_abb_residual(**kwargs)
        except Exception as exc:
            self.fail(f"temperature ABB residual is unavailable: {exc}")

    def call_adiabatic(self, **kwargs):
        self.require_boundary_module()
        try:
            return adiabatic_bb_residual(**kwargs)
        except Exception as exc:
            self.fail(f"adiabatic BB residual is unavailable: {exc}")

    def test_quadratic_abb_chain_maps_physical_flux_times_nominal_even_ghost(self):
        sigma_even, sigma_odd, chi, pi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi"
        )
        residual = self.call_temperature(
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=chi,
            cs2=Rational(1, 3),
            pressure_ratio=pi,
            regime="steady_1d_quadratic",
            normalization="transformed_cde_chain",
            scalar_gauge="half_source",
            geometry="flat_grid_aligned_halfway",
        )
        physical_product = (
            (1 - chi)
            * Rational(1, 3)
            / (Rational(1, 3) + pi)
            * sigma_odd
            * sigma_even
        )
        self.assertEqual(
            residual.solve_zero_conditions()[physical_product], Rational(3, 16)
        )
        self.assertEqual(
            simplify(
                residual.coefficients["temperature_curvature"]
                - (
                    16 * (1 - chi) * sigma_even * sigma_odd
                    - 3 * (1 + 3 * pi)
                )
                / 36
            ),
            0,
        )
        mapping = residual.parameter_mapping
        self.assertEqual(mapping["scalar_even_ghost_nominal"], sigma_even)
        self.assertEqual(mapping["scalar_odd_ghost_nominal"], sigma_odd)
        self.assertEqual(mapping["scalar_flux_physical"], physical_product / sigma_even)
        self.assertEqual(classify_magic(residual).status, "restricted_calibration")

    def test_temperature_primary_api_consumes_executable_population_chain(self):
        generator = getattr(
            boundary_module, "_temperature_abb_population_chain", None
        )
        self.assertTrue(callable(generator))
        sigma_even, sigma_odd, chi, pi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi"
        )
        generated_results = []

        def record_chain(*args, **kwargs):
            result = generator(*args, **kwargs)
            generated_results.append(result)
            return result

        with patch.object(
            boundary_module, generator.__name__, side_effect=record_chain
        ) as traced:
            residual = self.call_temperature(
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                pressure_ratio=pi,
                regime="steady_1d_quadratic",
                normalization="transformed_cde_chain",
                scalar_gauge="half_source",
                geometry="flat_grid_aligned_halfway",
            )
        self.assertEqual(traced.call_count, 1)
        generated = generated_results[0]
        self.assertTrue(generated.equations)
        self.assertEqual(residual.coefficients, generated.coefficients)

    def test_local_feedback_homogeneous_chain_is_not_extrapolated_from_external_source(self):
        sigma_even, sigma_odd, chi, pi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi"
        )
        local = self.call_temperature(
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=chi,
            cs2=Rational(1, 3),
            pressure_ratio=pi,
            regime="steady_1d_quadratic_local_feedback",
            normalization="transformed_cde_chain",
            scalar_gauge="half_source",
            geometry="flat_grid_aligned_halfway",
        )
        physical_product = (
            (1 - chi)
            * Rational(1, 3)
            / (Rational(1, 3) + pi)
            * sigma_odd
            * sigma_even
        )
        self.assertEqual(
            local.solve_zero_conditions()[physical_product],
            Rational(3, 16),
        )
        self.assertEqual(
            simplify(
                local.solve_zero_conditions()[physical_product].subs(pi, 0)
                - Rational(3, 16)
            ),
            0,
        )
        self.assertIn("local_nonequilibrium_feedback", local.assumptions)
        try:
            direct = manufactured_temperature_stencil(
                profile="dirichlet_quadratic_local_feedback",
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                pressure_ratio=pi,
                geometry="flat_grid_aligned_halfway",
            )
        except Exception as exc:
            self.fail(f"local-feedback homogeneous ABB chain is unavailable: {exc}")
        self.assertEqual(
            simplify(
                direct.first_nonzero_coefficient
                - local.coefficients["temperature_curvature"]
            ),
            0,
        )

    def test_general_abb_retains_full_jet_names_and_rate_independent_pressure_term(self):
        sigma_even, sigma_odd, chi, pi, temperature_gradient = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi T_n"
        )
        residual = self.call_temperature(
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=chi,
            cs2=Rational(1, 3),
            pressure_ratio=pi,
            normal_temperature_gradient=temperature_gradient,
            regime="general_source_aware",
            normalization="transformed_cde_chain",
            scalar_gauge="half_source",
            geometry="flat_grid_aligned_halfway",
        )
        required = {
            "temperature",
            "normal_temperature_gradient",
            "tangential_temperature_gradient",
            "temperature_curvature",
            "mixed_temperature_curvature",
            "tangential_temperature_curvature",
            "wall_time",
            "wall_time_time",
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
            "normal_velocity",
            "tangential_velocity",
            "pressure_temperature_gradient_source",
            "temperature_force_source",
            "heat_velocity_source",
        }
        self.assertTrue(required.issubset(residual.coefficients))
        self.assertEqual(
            residual.coefficients["normal_pressure_gradient"],
            -temperature_gradient / 4,
        )
        self.assertEqual(
            classify_magic(residual).status, "boundary_correction_required"
        )

    def test_direct_quadratic_temperature_stencil_is_independent_and_agrees(self):
        self.require_boundary_module()
        sigma_even, sigma_odd, chi, pi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi"
        )
        expected = (
            16 * (1 - chi) * sigma_even * sigma_odd
            - 3 * (1 + 3 * pi)
        ) / 36
        with patch(
            "tools.derivation.boundary.temperature_abb_residual",
            side_effect=AssertionError("manufactured route called symbolic solver"),
        ):
            check = manufactured_temperature_stencil(
                profile="dirichlet_quadratic",
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                pressure_ratio=pi,
                geometry="flat_grid_aligned_halfway",
            )
        self.assertEqual(simplify(check.first_nonzero_coefficient - expected), 0)
        self.assertTrue(check.includes_pressure_wall_equilibrium)
        self.assertTrue(check.uses_half_source_reconstruction)

    def test_affine_abb_and_variable_pressure_manufactured_probes_are_distinct(self):
        self.require_boundary_module()
        sigma_even, sigma_odd, normal_gradient, pressure_gradient = symbols(
            "sigma_g_plus sigma_g_minus T_n pi_n"
        )
        affine = manufactured_temperature_stencil(
            profile="dirichlet_affine",
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=0,
            cs2=Rational(1, 3),
            pressure_ratio=0,
            geometry="flat_grid_aligned_halfway",
        )
        variable_pressure = manufactured_temperature_stencil(
            profile="dirichlet_variable_pressure",
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=0,
            cs2=Rational(1, 3),
            pressure_ratio=0,
            normal_temperature_gradient=normal_gradient,
            normal_pressure_gradient=pressure_gradient,
            geometry="flat_grid_aligned_halfway",
        )
        self.assertEqual(affine.first_nonzero_coefficient, 0)
        self.assertEqual(
            variable_pressure.first_nonzero_coefficient,
            -normal_gradient * pressure_gradient / 4,
        )

    def test_adiabatic_chain_starts_from_kinetic_odd_flux_not_a_magic_product(self):
        sigma_even, sigma_odd, chi, pi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi"
        )
        residual = self.call_adiabatic(
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=chi,
            cs2=Rational(1, 3),
            pressure_ratio=pi,
            regime="kinetic_reflected_flux",
            normalization="transformed_cde_chain",
            scalar_gauge="half_source",
            geometry="flat_grid_aligned_halfway",
        )
        self.assertEqual(
            simplify(
                residual.coefficients["normal_temperature_gradient"]
                + (1 - chi) * sigma_odd / 3
            ),
            0,
        )
        self.assertEqual(residual.calibration_parameters, ())
        self.assertEqual(residual.solve_zero_conditions(), {})
        self.assertEqual(residual.parameter_mapping["scalar_even_ghost_nominal"], sigma_even)
        self.assertEqual(residual.parameter_mapping["scalar_odd_ghost_nominal"], sigma_odd)

    def test_adiabatic_primary_api_consumes_executable_population_chain(self):
        generator = getattr(boundary_module, "_adiabatic_population_chain", None)
        self.assertTrue(callable(generator))
        sigma_even, sigma_odd, chi, pi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi"
        )
        generated_results = []

        def record_chain(*args, **kwargs):
            result = generator(*args, **kwargs)
            generated_results.append(result)
            return result

        with patch.object(
            boundary_module, generator.__name__, side_effect=record_chain
        ) as traced:
            residual = self.call_adiabatic(
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                pressure_ratio=pi,
                regime="general_source_aware",
                normalization="transformed_cde_chain",
                scalar_gauge="half_source",
                geometry="flat_grid_aligned_halfway",
            )
        self.assertEqual(traced.call_count, 1)
        generated = generated_results[0]
        self.assertTrue(generated.equations)
        self.assertEqual(residual.coefficients, generated.coefficients)
        self.assertTrue(any(equation.has(sigma_even) for equation in generated.equations))
        self.assertTrue(any(equation.has(pi) for equation in generated.equations))

        pressure_perturbed = generator(
            sigma_even,
            sigma_odd,
            chi,
            Rational(1, 3),
            pi + Rational(1, 7),
            "general_source_aware",
        )
        even_rate_perturbed = generator(
            sigma_even + Rational(1, 5),
            sigma_odd,
            chi,
            Rational(1, 3),
            pi,
            "general_source_aware",
        )
        self.assertNotEqual(generated.equations, pressure_perturbed.equations)
        self.assertNotEqual(generated.equations, even_rate_perturbed.equations)
        self.assertEqual(
            simplify(
                generated.coefficients["tangential_temperature_curvature"]
                - pressure_perturbed.coefficients[
                    "tangential_temperature_curvature"
                ]
            ),
            0,
        )
        self.assertEqual(
            simplify(
                generated.coefficients["normal_temperature_gradient"]
                - even_rate_perturbed.coefficients[
                    "normal_temperature_gradient"
                ]
            ),
            0,
        )

    def test_general_adiabatic_retains_post_neumann_jets_and_needs_correction(self):
        sigma_even, sigma_odd, chi, pi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi"
        )
        residual = self.call_adiabatic(
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=chi,
            cs2=Rational(1, 3),
            pressure_ratio=pi,
            regime="general_source_aware",
            normalization="transformed_cde_chain",
            scalar_gauge="half_source",
            geometry="flat_grid_aligned_halfway",
        )
        retained_after_neumann = {
            "temperature_curvature",
            "tangential_temperature_curvature",
            "wall_time",
            "normal_pressure_gradient",
            "normal_force",
            "heat_source",
            "normal_heat_source_gradient",
            "tangential_heat_source_gradient",
        }
        self.assertTrue(retained_after_neumann.issubset(residual.coefficients))
        classification = classify_magic(residual)
        self.assertEqual(classification.status, "boundary_correction_required")
        self.assertIn(
            "tangential_temperature_curvature", classification.unsatisfied_jets
        )
        self.assertNotEqual(residual.coefficients["temperature_force_source"], 0)

    def test_direct_adiabatic_diagonal_pair_exposes_tangential_curvature(self):
        self.require_boundary_module()
        sigma_even, sigma_odd, chi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa"
        )
        with patch(
            "tools.derivation.boundary.adiabatic_bb_residual",
            side_effect=AssertionError("manufactured route called symbolic solver"),
        ):
            tangential = manufactured_temperature_stencil(
                profile="adiabatic_tangential_quadratic",
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                pressure_ratio=0,
                geometry="flat_grid_aligned_halfway",
            )
            try:
                force_only = manufactured_temperature_stencil(
                    profile="adiabatic_force_only",
                    sigma_g_plus=sigma_even,
                    sigma_g_minus=sigma_odd,
                    chi_kappa=chi,
                    cs2=Rational(1, 3),
                    pressure_ratio=0,
                    normal_force=symbols("F_n"),
                    geometry="flat_grid_aligned_halfway",
                )
                hydrostatic = manufactured_temperature_stencil(
                    profile="adiabatic_hydrostatic_pair",
                    sigma_g_plus=sigma_even,
                    sigma_g_minus=sigma_odd,
                    chi_kappa=chi,
                    cs2=Rational(1, 3),
                    pressure_ratio=0,
                    normal_force=symbols("F_n"),
                    normal_pressure_gradient=symbols("F_n"),
                    geometry="flat_grid_aligned_halfway",
                )
                uniform_heat = manufactured_temperature_stencil(
                    profile="adiabatic_uniform_heat",
                    sigma_g_plus=sigma_even,
                    sigma_g_minus=sigma_odd,
                    chi_kappa=chi,
                    cs2=Rational(1, 3),
                    pressure_ratio=0,
                    heat_source=symbols("Q"),
                    geometry="flat_grid_aligned_halfway",
                )
            except Exception as exc:
                self.fail(f"separate adiabatic source probes are unavailable: {exc}")
        self.assertEqual(
            simplify(
                tangential.first_nonzero_coefficient
                - (1 - chi) * sigma_odd / 9
            ),
            0,
        )
        self.assertEqual(force_only.first_nonzero_coefficient, sigma_odd)
        self.assertEqual(hydrostatic.first_nonzero_coefficient, 0)
        self.assertEqual(uniform_heat.first_nonzero_coefficient, 0)
        self.assertTrue(hydrostatic.uses_half_source_reconstruction)
        self.assertTrue(uniform_heat.uses_half_source_reconstruction)
        self.assertTrue(uniform_heat.derivation.has(symbols("Q")))
        self.assertIn(
            "heat_source", getattr(uniform_heat, "consumed_inputs", ())
        )

    def test_nonstationary_abb_and_adiabatic_time_curvature_are_explicitly_unresolved(self):
        self.require_boundary_module()
        sigma_even, sigma_odd, chi, pi = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi"
        )
        wall_time, normal_curvature = symbols("T_t T_nn")
        try:
            abb = manufactured_temperature_stencil(
                profile="dirichlet_nonstationary_abb",
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                pressure_ratio=pi,
                wall_time_derivative=wall_time,
                geometry="flat_grid_aligned_halfway",
            )
            adiabatic = manufactured_temperature_stencil(
                profile="adiabatic_wall_time_normal_curvature",
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                pressure_ratio=pi,
                wall_time_derivative=wall_time,
                normal_temperature_curvature=normal_curvature,
                geometry="flat_grid_aligned_halfway",
            )
        except Exception as exc:
            self.fail(f"unsteady manufactured coverage is unavailable: {exc}")
        self.assertEqual(
            getattr(abb, "verification_status", None),
            "unsupported_unresolved",
        )
        self.assertEqual(
            getattr(adiabatic, "verification_status", None),
            "unsupported_unresolved",
        )
        self.assertTrue(abb.derivation.has(wall_time))
        self.assertTrue(adiabatic.derivation.has(wall_time, normal_curvature))
        self.assertTrue(abb.first_nonzero_coefficient != 0)
        self.assertTrue(adiabatic.first_nonzero_coefficient != 0)


class MixedCornerBoundaryTests(unittest.TestCase):
    """A shared diagonal population cannot be independently overwritten twice."""

    def require_boundary_module(self) -> None:
        self.assertIsNone(
            BOUNDARY_IMPORT_ERROR,
            f"corner boundary derivation is unavailable: {BOUNDARY_IMPORT_ERROR}",
        )

    def call_corner(self, **kwargs):
        self.require_boundary_module()
        try:
            return mixed_corner_residual(**kwargs)
        except Exception as exc:
            self.fail(f"mixed-corner residual is unavailable: {exc}")

    def test_one_diagonal_unknown_two_wall_equations_are_rank_incompatible(self):
        sigma_even, sigma_odd, chi, pi, h = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi h", nonzero=True
        )
        residual = self.call_corner(
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=chi,
            cs2=Rational(1, 3),
            pressure_ratio=pi,
            grid_spacing=h,
            corner_type="dirichlet_adiabatic",
            normalization="transformed_cde_chain",
            scalar_gauge="half_source",
            geometry="grid_aligned_right_angle_corner",
        )
        mapping = residual.parameter_mapping
        self.assertEqual(mapping["shared_diagonal_unknown_count"], 1)
        self.assertEqual(mapping["wall_equation_count"], 2)
        self.assertEqual(mapping["closure_rank"], 1)
        self.assertEqual(mapping["generic_augmented_rank"], 2)
        self.assertEqual(classify_magic(residual).status, "corner_closure_conflict")
        self.assertEqual(residual.calibration_parameters, ())

    def test_corner_primary_api_consumes_finite_assignment_chain(self):
        generator = getattr(boundary_module, "_corner_population_chain", None)
        self.assertTrue(callable(generator))
        sigma_even, sigma_odd, chi, pi, h = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi h", nonzero=True
        )
        generated_results = []

        def record_chain(*args, **kwargs):
            result = generator(*args, **kwargs)
            generated_results.append(result)
            return result

        with patch.object(
            boundary_module, generator.__name__, side_effect=record_chain
        ) as traced:
            residual = self.call_corner(
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=chi,
                cs2=Rational(1, 3),
                pressure_ratio=pi,
                grid_spacing=h,
                corner_type="dirichlet_adiabatic",
                normalization="transformed_cde_chain",
                scalar_gauge="half_source",
                geometry="grid_aligned_right_angle_corner",
            )
        self.assertEqual(traced.call_count, 1)
        generated = generated_results[0]
        self.assertTrue(generated.equations)
        self.assertTrue(generated.assignments)
        self.assertEqual(residual.coefficients, generated.coefficients)

    def test_corner_table_retains_2d_temperature_and_source_jets(self):
        sigma_even, sigma_odd, chi, pi, h = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi h", nonzero=True
        )
        residual = self.call_corner(
            sigma_g_plus=sigma_even,
            sigma_g_minus=sigma_odd,
            chi_kappa=chi,
            cs2=Rational(1, 3),
            pressure_ratio=pi,
            grid_spacing=h,
            corner_type="dirichlet_adiabatic",
            normalization="transformed_cde_chain",
            scalar_gauge="half_source",
            geometry="grid_aligned_right_angle_corner",
        )
        required = {
            "temperature",
            "x_temperature_gradient",
            "y_temperature_gradient",
            "xx_temperature_curvature",
            "mixed_temperature_curvature",
            "yy_temperature_curvature",
            "wall_time",
            "wall_time_x",
            "wall_time_y",
            "pressure",
            "x_pressure_gradient",
            "y_pressure_gradient",
            "x_force",
            "y_force",
            "heat_source",
            "x_velocity",
            "y_velocity",
        }
        self.assertTrue(required.issubset(residual.coefficients))
        mixed = residual.coefficients["mixed_temperature_curvature"]
        self.assertNotEqual(mixed, 0)
        self.assertFalse(
            bool(mixed.free_symbols & {sigma_even, sigma_odd, chi}),
            "the equilibrium-level T_xy compatibility defect is not tunable",
        )

    def test_direct_txy_corner_is_order_dependent_and_not_simultaneous(self):
        self.require_boundary_module()
        sigma_even, sigma_odd, chi, pi, h, mixed = symbols(
            "sigma_g_plus sigma_g_minus chi_kappa pi h T_xy", nonzero=True
        )
        with patch(
            "tools.derivation.boundary.mixed_corner_residual",
            side_effect=AssertionError("manufactured corner called residual solver"),
        ):
            try:
                check = manufactured_corner_stencil(
                    profile="mixed_txy",
                    sigma_g_plus=sigma_even,
                    sigma_g_minus=sigma_odd,
                    chi_kappa=chi,
                    cs2=Rational(1, 3),
                    pressure_ratio=pi,
                    grid_spacing=h,
                    mixed_temperature_curvature=mixed,
                    geometry="grid_aligned_right_angle_corner",
                )
            except Exception as exc:
                self.fail(f"manufactured mixed-corner stencil is unavailable: {exc}")
        self.assertEqual(check.closure_rank, 1)
        self.assertEqual(check.augmented_rank, 2)
        self.assertNotEqual(
            simplify(check.dirichlet_then_adiabatic - check.adiabatic_then_dirichlet),
            0,
        )
        self.assertEqual(
            simplify(
                check.first_nonzero_coefficient
                + (1 + 3 * pi) / 72
            ),
            0,
        )
        self.assertFalse(check.simultaneously_satisfied)

    def test_corner_counts_shared_source_once_and_uses_diagonal_distance(self):
        self.require_boundary_module()
        sigma_even, sigma_odd, h = symbols(
            "sigma_g_plus sigma_g_minus h", nonzero=True
        )
        try:
            check = manufactured_corner_stencil(
                profile="mixed_txy",
                sigma_g_plus=sigma_even,
                sigma_g_minus=sigma_odd,
                chi_kappa=0,
                cs2=Rational(1, 3),
                pressure_ratio=0,
                grid_spacing=h,
                mixed_temperature_curvature=symbols("T_xy"),
                geometry="grid_aligned_right_angle_corner",
            )
        except Exception as exc:
            self.fail(f"corner source/distance audit is unavailable: {exc}")
        self.assertEqual(check.naive_overwrite_heat_source_count, 2)
        self.assertEqual(check.shared_population_heat_source_count, 1)
        assignment_ids = getattr(check, "source_assignment_population_ids", ())
        shared_ids = getattr(check, "shared_source_population_ids", ())
        self.assertEqual(
            check.naive_overwrite_heat_source_count, len(assignment_ids)
        )
        self.assertEqual(
            check.shared_population_heat_source_count, len(set(shared_ids))
        )
        self.assertEqual(check.diagonal_wall_distance, sqrt(2) * h / 2)


class ClassificationSafetyTests(unittest.TestCase):
    """A universal label must be proved by every retained coefficient row."""

    def test_incompatible_rows_cannot_be_declared_universal(self):
        self.assertIsNone(BOUNDARY_IMPORT_ERROR)
        product = symbols("candidate_product")
        residual = BoundaryResidual(
            boundary="adversarial_universal_claim",
            coefficients={"jet_a": product - 1, "jet_b": product - 2},
            assumptions=("adversarial_test",),
            parameter_mapping={"candidate_product": product},
            calibration_parameters=(product,),
            recommended_status="universal_magic",
            derivation="two incompatible retained jets",
        )
        classification = classify_magic(residual)
        self.assertEqual(classification.status, "no_single_magic")
        self.assertEqual(set(classification.unsatisfied_jets), {"jet_a", "jet_b"})

    def test_one_condition_zeroing_every_row_can_remain_universal(self):
        self.assertIsNone(BOUNDARY_IMPORT_ERROR)
        product = symbols("candidate_product")
        residual = BoundaryResidual(
            boundary="consistent_universal_claim",
            coefficients={"jet_a": product - 1, "jet_b": 2 * (product - 1)},
            assumptions=("adversarial_test",),
            parameter_mapping={"candidate_product": product},
            calibration_parameters=(product,),
            recommended_status="universal_magic",
            derivation="one condition removes all retained jets",
        )
        self.assertEqual(classify_magic(residual).status, "universal_magic")


if __name__ == "__main__":
    unittest.main()
