"""Effective-rate elimination and second-order CE residual tests."""

from __future__ import annotations

from dataclasses import replace
import unittest

from sympy import Integer, Rational, diff, simplify, symbols


_EFFECTIVE_RATES_IMPORT_ERROR: str | None = None

try:
    from tools.derivation.effective_rates import (
        EquilibriumMomentConstraints,
        SourceMomentConstraints,
        TrapezoidalFactors,
        bulk_effective_rate,
        effective_operator_blocks,
        scalar_flux_effective_rate,
        scalar_variable_pressure_residual,
        second_order_residual_table,
        shear_effective_rate,
    )
except (ImportError, ModuleNotFoundError) as exc:
    _EFFECTIVE_RATES_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"


class EffectiveRateTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.assertIsNone(
            _EFFECTIVE_RATES_IMPORT_ERROR,
            "effective-rate API is unavailable: "
            f"{_EFFECTIVE_RATES_IMPORT_ERROR}",
        )

        self.sigma, self.chi, self.cs2, self.dt = symbols(
            "sigma chi c_s2 dt", nonzero=True
        )
        self.pi = symbols("pi")

    def assert_expr_equal(self, actual, expected) -> None:
        self.assertEqual(simplify(actual - expected), 0)


class ScalarFluxEffectiveRateTests(EffectiveRateTestCase):
    def test_scalar_flux_feedback_is_reduced_from_direct_collision(self) -> None:
        nonequilibrium, affine = symbols("j_neq a_flux")
        result = scalar_flux_effective_rate(
            self.sigma,
            self.chi,
            pressure_ratio=self.pi,
            cs2=self.cs2,
            dt=self.dt,
        )
        nominal_rate = 1 / (self.sigma + Rational(1, 2))
        diffusivity = (1 - self.chi) * self.cs2 * self.dt * self.sigma
        gradient_closure = -(
            2 * nonequilibrium + self.dt * affine
        ) / (2 * diffusivity + self.dt * (self.cs2 + self.pi))
        direct_collision = (
            -nominal_rate * nonequilibrium
            + self.dt
            * (1 - nominal_rate / 2)
            * ((self.pi + self.chi * self.cs2) * gradient_closure + affine)
        )
        reduced_collision = (
            -result.rate_effective * nonequilibrium
            + self.dt * (1 - result.rate_effective / 2) * affine
        )

        self.assert_expr_equal(direct_collision, reduced_collision)
        self.assert_expr_equal(
            result.sigma_effective,
            (1 - self.chi)
            * self.cs2
            * self.sigma
            / (self.cs2 + self.pi),
        )

    def test_frozen_pressure_cancels_from_physical_diffusivity(self) -> None:
        result = scalar_flux_effective_rate(
            self.sigma,
            self.chi,
            pressure_ratio=self.pi,
            cs2=self.cs2,
            dt=self.dt,
        )

        self.assert_expr_equal(
            result.transport_coefficient,
            (1 - self.chi) * self.cs2 * self.dt * self.sigma,
        )
        self.assert_expr_equal(
            (self.cs2 + self.pi) * self.dt * result.sigma_effective,
            result.transport_coefficient,
        )

    def test_variable_pressure_residual_is_not_claimed_zero(self) -> None:
        delta_pi, grad_pi_dot_grad_t, laplacian_t = symbols(
            "delta_pi grad_pi_dot_grad_T laplacian_T"
        )
        result = scalar_flux_effective_rate(
            self.sigma,
            self.chi,
            pressure_ratio=self.pi,
            cs2=self.cs2,
            dt=self.dt,
        )
        residual = scalar_variable_pressure_residual(
            result,
            pressure_perturbation=delta_pi,
            grad_pressure_ratio_dot_grad_temperature=grad_pi_dot_grad_t,
            laplacian_temperature=laplacian_t,
        )

        self.assert_expr_equal(
            residual,
            self.dt
            * result.sigma_effective
            * (grad_pi_dot_grad_t + delta_pi * laplacian_t),
        )
        self.assertNotEqual(residual, 0)


class FlowEffectiveRateTests(EffectiveRateTestCase):
    def _assert_direct_mode_elimination(self, *, chi, mode_result) -> None:
        nonequilibrium, affine, rho0 = symbols("Pi_neq a_mode rho_0")
        nominal_rate = 1 / (self.sigma + Rational(1, 2))
        physical_modulus = (
            rho0 * self.cs2 * self.dt * mode_result.sigma_effective
        )
        strain_closure = -(
            2 * nonequilibrium + self.dt * affine
        ) / (4 * physical_modulus + 2 * self.dt * rho0 * self.cs2)
        direct_collision = (
            -nominal_rate * nonequilibrium
            + self.dt
            * (1 - nominal_rate / 2)
            * (affine + 2 * rho0 * self.cs2 * chi * strain_closure)
        )
        reduced_collision = (
            -mode_result.rate_effective * nonequilibrium
            + self.dt * (1 - mode_result.rate_effective / 2) * affine
        )
        self.assert_expr_equal(direct_collision, reduced_collision)

    def test_off_diagonal_and_deviatoric_diagonal_use_shear_block(self) -> None:
        for mode in ("off_diagonal", "deviatoric_diagonal"):
            with self.subTest(mode=mode):
                result = shear_effective_rate(
                    self.sigma,
                    self.chi,
                    cs2=self.cs2,
                    dt=self.dt,
                    mode=mode,
                )
                self._assert_direct_mode_elimination(
                    chi=self.chi,
                    mode_result=result,
                )
                self.assertEqual(result.mode, mode)
                self.assert_expr_equal(
                    result.sigma_effective,
                    (1 - self.chi) * self.sigma,
                )

    def test_trace_mode_uses_bulk_block_and_dimension_normalization(self) -> None:
        dimension = symbols("D", positive=True)
        result = bulk_effective_rate(
            self.sigma,
            self.chi,
            dimension=dimension,
            cs2=self.cs2,
            dt=self.dt,
        )

        self._assert_direct_mode_elimination(
            chi=self.chi,
            mode_result=result,
        )
        self.assertEqual(result.mode, "trace_bulk")
        self.assert_expr_equal(
            result.transport_coefficient,
            2
            * (1 - self.chi)
            * self.cs2
            * self.dt
            * self.sigma
            / dimension,
        )

    def test_default_bgk_feedback_limit_and_high_nominal_rate_limit(self) -> None:
        nominal = shear_effective_rate(self.sigma, Integer(0))
        self.assert_expr_equal(nominal.sigma_effective, self.sigma)
        self.assert_expr_equal(
            nominal.rate_effective,
            1 / (self.sigma + Rational(1, 2)),
        )

        small = symbols("eta", positive=True)
        low_transport = shear_effective_rate(
            Rational(1, 2),
            1 - small,
        )
        self.assert_expr_equal(low_transport.sigma_effective, small / 2)
        self.assert_expr_equal(low_transport.rate_effective, 2 / (1 + small))

    def test_trt_skeleton_splits_physical_and_ghost_blocks(self) -> None:
        sigma_f_minus, sigma_g_plus, sigma_g_minus = symbols(
            "sigma_f_minus sigma_g_plus sigma_g_minus"
        )
        blocks = effective_operator_blocks(
            sigma_f_plus=self.sigma,
            sigma_f_minus=sigma_f_minus,
            sigma_g_plus=sigma_g_plus,
            sigma_g_minus=sigma_g_minus,
            chi_s=self.chi,
            chi_b=self.chi,
            chi_kappa=self.chi,
            pressure_ratio=self.pi,
            cs2=self.cs2,
            dt=self.dt,
        )

        self.assert_expr_equal(
            blocks["flow_shear"].sigma_effective,
            (1 - self.chi) * self.sigma,
        )
        self.assert_expr_equal(
            blocks["scalar_flux"].sigma_effective,
            (1 - self.chi)
            * self.cs2
            * sigma_g_minus
            / (self.cs2 + self.pi),
        )
        self.assert_expr_equal(
            blocks["flow_odd_ghost"].sigma_nominal,
            sigma_f_minus,
        )
        self.assert_expr_equal(
            blocks["scalar_odd_ghost"].sigma_nominal,
            sigma_g_minus,
        )
        self.assertNotEqual(
            simplify(
                blocks["scalar_flux"].sigma_effective
                - blocks["scalar_odd_ghost"].sigma_nominal
            ),
            0,
        )

    def test_equal_nominal_rates_and_zero_feedback_recover_bgk_blocks(self) -> None:
        blocks = effective_operator_blocks(
            sigma_f_plus=self.sigma,
            sigma_f_minus=self.sigma,
            sigma_g_plus=self.sigma,
            sigma_g_minus=self.sigma,
            chi_s=Integer(0),
            chi_b=Integer(0),
            chi_kappa=Integer(0),
            pressure_ratio=Integer(0),
            cs2=self.cs2,
            dt=self.dt,
        )

        for name, block in blocks.items():
            with self.subTest(block=name):
                shift = getattr(
                    block,
                    "sigma_effective",
                    block.sigma_nominal,
                )
                self.assert_expr_equal(shift, self.sigma)


class SecondOrderResidualTests(EffectiveRateTestCase):
    def test_default_named_moments_generate_all_second_order_cancellations(self) -> None:
        s_f_plus, s_f_minus, s_g_plus, s_g_minus = symbols(
            "s_f_plus s_f_minus s_g_plus s_g_minus", nonzero=True
        )
        table = second_order_residual_table(
            s_f_plus=s_f_plus,
            s_f_minus=s_f_minus,
            s_g_plus=s_g_plus,
            s_g_minus=s_g_minus,
        )

        for name in (
            "p_grad_T",
            "T_F",
            "Q_u",
            "u_F",
            "d_t_F",
            "d_t_Q",
            "s_f_minus_transport",
            "s_g_plus_transport",
        ):
            with self.subTest(name=name):
                self.assertEqual(simplify(table[name]), 0)
        self.assertEqual(table["first_omitted_epsilon_order"], 3)
        self.assertEqual(table["first_omitted_mach_order"], 3)

    def test_removing_each_inverse_design_source_moment_exposes_residual(self) -> None:
        defaults = SourceMomentConstraints()
        cases = {
            "p_grad_T": "scalar_first_p_grad_t",
            "T_F": "scalar_first_t_force",
            "Q_u": "scalar_first_q_velocity",
            "u_F": "flow_second_u_force",
        }

        for residual_name, field_name in cases.items():
            with self.subTest(residual=residual_name):
                perturbed = replace(defaults, **{field_name: Integer(0)})
                residual = simplify(
                    second_order_residual_table(sources=perturbed)[residual_name]
                )
                self.assertNotEqual(residual, 0)

    def test_wrong_parity_specific_trapezoidal_factors_expose_residuals(self) -> None:
        s_f_plus, s_f_minus, s_g_plus, s_g_minus = symbols(
            "s_f_plus s_f_minus s_g_plus s_g_minus", nonzero=True
        )
        correct = TrapezoidalFactors.from_rates(
            s_f_plus=s_f_plus,
            s_f_minus=s_f_minus,
            s_g_plus=s_g_plus,
            s_g_minus=s_g_minus,
        )
        cases = (
            (replace(correct, scalar_odd=Integer(1)), "p_grad_T"),
            (replace(correct, flow_even=Integer(1)), "u_F"),
            (replace(correct, flow_odd=Integer(1)), "d_t_F"),
            (replace(correct, scalar_even=Integer(1)), "d_t_Q"),
        )

        for factors, residual_name in cases:
            with self.subTest(residual=residual_name):
                table = second_order_residual_table(
                    s_f_plus=s_f_plus,
                    s_f_minus=s_f_minus,
                    s_g_plus=s_g_plus,
                    s_g_minus=s_g_minus,
                    trapezoidal=factors,
                )
                self.assertNotEqual(simplify(table[residual_name]), 0)

    def test_wrong_conserved_source_factor_has_exact_half_source_residual(self) -> None:
        s_f_plus, s_f_minus, s_g_plus, s_g_minus = symbols(
            "s_f_plus s_f_minus s_g_plus s_g_minus", nonzero=True
        )
        correct = TrapezoidalFactors.from_rates(
            s_f_plus=s_f_plus,
            s_f_minus=s_f_minus,
            s_g_plus=s_g_plus,
            s_g_minus=s_g_minus,
        )
        wrong_force = second_order_residual_table(
            s_f_plus=s_f_plus,
            s_f_minus=s_f_minus,
            s_g_plus=s_g_plus,
            s_g_minus=s_g_minus,
            trapezoidal=replace(correct, flow_odd=Integer(1)),
        )
        wrong_heat = second_order_residual_table(
            s_f_plus=s_f_plus,
            s_f_minus=s_f_minus,
            s_g_plus=s_g_plus,
            s_g_minus=s_g_minus,
            trapezoidal=replace(correct, scalar_even=Integer(1)),
        )
        wrong_scalar_flux = second_order_residual_table(
            s_f_plus=s_f_plus,
            s_f_minus=s_f_minus,
            s_g_plus=s_g_plus,
            s_g_minus=s_g_minus,
            trapezoidal=replace(correct, scalar_odd=Integer(1)),
        )
        wrong_flow_stress = second_order_residual_table(
            s_f_plus=s_f_plus,
            s_f_minus=s_f_minus,
            s_g_plus=s_g_plus,
            s_g_minus=s_g_minus,
            trapezoidal=replace(correct, flow_even=Integer(1)),
        )

        self.assert_expr_equal(wrong_force["T_F"], s_f_minus / 2)
        self.assert_expr_equal(wrong_force["u_F"], s_f_minus / 2)
        self.assert_expr_equal(wrong_force["d_t_F"], Rational(1, 2))
        self.assert_expr_equal(wrong_heat["Q_u"], s_g_plus / 2)
        self.assert_expr_equal(wrong_heat["d_t_Q"], Rational(1, 2))
        for name in ("p_grad_T", "T_F", "Q_u"):
            with self.subTest(block="scalar_odd", residual=name):
                self.assert_expr_equal(
                    wrong_scalar_flux[name],
                    s_g_minus / (s_g_minus - 2),
                )
        self.assert_expr_equal(
            wrong_flow_stress["u_F"],
            s_f_plus / (s_f_plus - 2),
        )

    def test_transport_blocks_exclude_opposite_parity_rates(self) -> None:
        s_f_plus, s_f_minus, s_g_plus, s_g_minus = symbols(
            "s_f_plus s_f_minus s_g_plus s_g_minus", nonzero=True
        )
        table = second_order_residual_table(
            s_f_plus=s_f_plus,
            s_f_minus=s_f_minus,
            s_g_plus=s_g_plus,
            s_g_minus=s_g_minus,
        )

        self.assertEqual(diff(table["nu_shear"], s_f_minus), 0)
        self.assertEqual(diff(table["kappa"], s_g_plus), 0)
        self.assertNotEqual(diff(table["nu_shear"], s_f_plus), 0)
        self.assertNotEqual(diff(table["kappa"], s_g_minus), 0)

    def test_equilibrium_constraints_are_live_inputs_not_target_pde_fixtures(self) -> None:
        perturbed = replace(
            EquilibriumMomentConstraints(),
            scalar_second_pressure=Integer(2),
        )
        table = second_order_residual_table(equilibrium=perturbed)

        self.assertNotEqual(simplify(table["p_grad_T"]), 0)


if __name__ == "__main__":
    unittest.main()
