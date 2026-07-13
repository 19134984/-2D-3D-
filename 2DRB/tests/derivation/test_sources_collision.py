"""Exact source-parity and transformed TRT collision tests."""

from __future__ import annotations

import unittest

from sympy import Matrix, Rational, simplify, symbols, zeros

from tools.derivation.lattice import d2q9, parity_projectors, raw_moment


_SOURCE_COLLISION_IMPORT_ERROR: str | None = None

try:
    from tools.derivation.collision import (
        bgk_collision,
        reconstruct_momentum,
        reconstruct_scalar,
        trt_collision,
    )
    from tools.derivation.sources import (
        flow_source,
        scalar_source,
        source_moment_table,
    )
except (ImportError, ModuleNotFoundError) as exc:
    _SOURCE_COLLISION_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"


class SourceCollisionTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.assertIsNone(
            _SOURCE_COLLISION_IMPORT_ERROR,
            "source/collision API is unavailable: "
            f"{_SOURCE_COLLISION_IMPORT_ERROR}",
        )

        self.lattice = d2q9()
        self.ux, self.uy = symbols("u_x u_y")
        self.fx, self.fy = symbols("F_x F_y")
        self.rho0, self.pressure, self.temperature = symbols("rho_0 p T")
        self.q, self.chi_s, self.chi_b, self.chi_kappa = symbols(
            "Q chi_s chi_b chi_kappa"
        )
        self.grad_tx, self.grad_ty = symbols("gradT_x gradT_y")
        self.sxx, self.sxy, self.syy = symbols("S_xx S_xy S_yy")

        self.velocity = Matrix([self.ux, self.uy])
        self.force = Matrix([self.fx, self.fy])
        self.grad_temperature = Matrix([self.grad_tx, self.grad_ty])
        self.strain = Matrix(
            [
                [self.sxx, self.sxy],
                [self.sxy, self.syy],
            ]
        )

    def flow_terms(self):
        return flow_source(
            self.lattice,
            velocity=self.velocity,
            force=self.force,
            rho0=self.rho0,
            strain=self.strain,
            chi_s=self.chi_s,
            chi_b=self.chi_b,
        )

    def scalar_terms(self):
        return scalar_source(
            self.lattice,
            velocity=self.velocity,
            force=self.force,
            pressure=self.pressure,
            temperature=self.temperature,
            grad_temperature=self.grad_temperature,
            rho0=self.rho0,
            heat_source=self.q,
            chi_kappa=self.chi_kappa,
        )

    def vector_moment(self, populations: Matrix) -> Matrix:
        return Matrix(
            [
                raw_moment(populations, self.lattice.velocities, (1, 0)),
                raw_moment(populations, self.lattice.velocities, (0, 1)),
            ]
        )

    def assert_zero_matrix(self, value: Matrix) -> None:
        self.assertEqual(value.applyfunc(simplify), zeros(*value.shape))

    def assert_expr_equal(self, actual, expected) -> None:
        self.assertEqual(simplify(actual - expected), 0)


class SourceMomentTests(SourceCollisionTestCase):
    def test_sources_are_exact_even_odd_projections(self) -> None:
        for terms in (self.flow_terms(), self.scalar_terms()):
            self.assert_zero_matrix(terms.raw - terms.plus - terms.minus)
            for index, opposite in enumerate(self.lattice.opposite):
                self.assert_expr_equal(terms.plus[index], terms.plus[opposite])
                self.assert_expr_equal(terms.minus[index], -terms.minus[opposite])

    def test_flow_source_required_low_order_moments(self) -> None:
        moments = self.flow_terms().moment_tables

        self.assertEqual(moments["plus"][(0, 0)], 0)
        self.assertEqual(moments["minus"][(0, 0)], 0)
        self.assertEqual(moments["plus"][(1, 0)], 0)
        self.assertEqual(moments["plus"][(0, 1)], 0)
        self.assertEqual(moments["minus"][(1, 0)], self.fx)
        self.assertEqual(moments["minus"][(0, 1)], self.fy)
        for powers in ((2, 0), (1, 1), (0, 2)):
            self.assertEqual(moments["minus"][powers], 0)

    def test_eq13_positive_dimensionless_hermite_contraction(self) -> None:
        moments = self.flow_terms().moment_tables["plus"]
        cs2 = self.lattice.cs2
        trace_strain = self.sxx + self.syy
        a_xx = self.chi_s * self.sxx + Rational(1, 2) * (
            self.chi_b - self.chi_s
        ) * trace_strain
        a_xy = self.chi_s * self.sxy
        a_yy = self.chi_s * self.syy + Rational(1, 2) * (
            self.chi_b - self.chi_s
        ) * trace_strain

        self.assert_expr_equal(
            moments[(2, 0)],
            2 * self.ux * self.fx + 2 * self.rho0 * cs2 * a_xx,
        )
        self.assert_expr_equal(
            moments[(1, 1)],
            self.ux * self.fy
            + self.uy * self.fx
            + 2 * self.rho0 * cs2 * a_xy,
        )
        self.assert_expr_equal(
            moments[(0, 2)],
            2 * self.uy * self.fy + 2 * self.rho0 * cs2 * a_yy,
        )

    def test_scalar_source_required_low_order_moments(self) -> None:
        moments = self.scalar_terms().moment_tables
        cs2 = self.lattice.cs2
        flux = (
            self.pressure * self.grad_temperature
            + self.temperature * self.force
        ) / self.rho0 + self.q * self.velocity + self.chi_kappa * cs2 * (
            self.grad_temperature
        )

        self.assertEqual(moments["plus"][(0, 0)], self.q)
        self.assertEqual(moments["plus"][(1, 0)], 0)
        self.assertEqual(moments["plus"][(0, 1)], 0)
        self.assertEqual(moments["plus"][(2, 0)], cs2 * self.q)
        self.assertEqual(moments["plus"][(1, 1)], 0)
        self.assertEqual(moments["plus"][(0, 2)], cs2 * self.q)
        self.assertEqual(moments["minus"][(0, 0)], 0)
        self.assert_expr_equal(moments["minus"][(1, 0)], flux[0])
        self.assert_expr_equal(moments["minus"][(0, 1)], flux[1])
        for powers in ((2, 0), (1, 1), (0, 2)):
            self.assertEqual(moments["minus"][powers], 0)

    def test_nonzero_third_and_fourth_raw_source_moments_are_retained(self) -> None:
        flow_terms = self.flow_terms()
        flow = flow_terms.moment_tables
        scalar = self.scalar_terms().moment_tables
        cs2 = self.lattice.cs2
        trace_strain = self.sxx + self.syy
        a_trace = self.chi_b * trace_strain
        u_dot_f = self.ux * self.fx + self.uy * self.fy
        scalar_flux = (
            self.pressure * self.grad_temperature
            + self.temperature * self.force
        ) / self.rho0 + self.q * self.velocity + self.chi_kappa * cs2 * (
            self.grad_temperature
        )

        self.assertEqual(
            flow["raw"], source_moment_table(flow_terms.raw, self.lattice)
        )

        self.assertEqual(flow["minus"][(3, 0)], self.fx)
        self.assertEqual(flow["minus"][(2, 1)], cs2 * self.fy)
        self.assertEqual(flow["minus"][(1, 2)], cs2 * self.fx)
        self.assertEqual(flow["minus"][(0, 3)], self.fy)
        self.assertEqual(flow["plus"][(4, 0)], flow["plus"][(2, 0)])
        self.assertEqual(flow["plus"][(3, 1)], flow["plus"][(1, 1)])
        self.assert_expr_equal(
            flow["plus"][(2, 2)],
            2 * cs2 * u_dot_f + 2 * self.rho0 * cs2**2 * a_trace,
        )
        self.assertEqual(flow["plus"][(1, 3)], flow["plus"][(1, 1)])
        self.assertEqual(flow["plus"][(0, 4)], flow["plus"][(0, 2)])

        self.assertEqual(scalar["plus"][(4, 0)], 3 * cs2**2 * self.q)
        self.assertEqual(scalar["plus"][(2, 2)], cs2**2 * self.q)
        self.assertEqual(scalar["plus"][(0, 4)], 3 * cs2**2 * self.q)
        self.assert_expr_equal(scalar["minus"][(3, 0)], scalar_flux[0])
        self.assert_expr_equal(scalar["minus"][(2, 1)], cs2 * scalar_flux[1])
        self.assert_expr_equal(scalar["minus"][(1, 2)], cs2 * scalar_flux[0])
        self.assert_expr_equal(scalar["minus"][(0, 3)], scalar_flux[1])
        for powers in ((4, 0), (3, 1), (2, 2), (1, 3), (0, 4)):
            self.assertEqual(flow["minus"][powers], 0)
            self.assertEqual(scalar["minus"][powers], 0)


class CollisionTests(SourceCollisionTestCase):
    def test_trt_collision_matches_operator_trapezoidal_formula(self) -> None:
        h_tilde = Matrix(symbols("h_tilde_0:9"))
        h_eq = Matrix(symbols("h_eq_0:9"))
        source_plus = Matrix(symbols("source_plus_0:9"))
        source_minus = Matrix(symbols("source_minus_0:9"))
        s_plus, s_minus, dt = symbols("s_plus s_minus dt")
        p_plus, p_minus = parity_projectors(self.lattice)

        actual = trt_collision(
            self.lattice,
            h_tilde=h_tilde,
            h_eq=h_eq,
            source_plus=source_plus,
            source_minus=source_minus,
            s_plus=s_plus,
            s_minus=s_minus,
            dt=dt,
        )
        expected = (
            h_tilde
            - s_plus * p_plus * (h_tilde - h_eq)
            - s_minus * p_minus * (h_tilde - h_eq)
            + dt
            * (
                (1 - s_plus / 2) * source_plus
                + (1 - s_minus / 2) * source_minus
            )
        )

        self.assert_zero_matrix(actual - expected)

    def test_half_source_macroscopic_reconstructions(self) -> None:
        f_tilde = Matrix(symbols("f_tilde_0:9"))
        g_tilde = Matrix(symbols("g_tilde_0:9"))
        dt = symbols("dt")

        expected_momentum = self.vector_moment(f_tilde) + dt * self.force / 2
        self.assert_zero_matrix(
            reconstruct_momentum(
                self.lattice, f_tilde=f_tilde, force=self.force, dt=dt
            )
            - expected_momentum
        )
        self.assertEqual(
            reconstruct_scalar(g_tilde=g_tilde, heat_source=self.q, dt=dt),
            sum(g_tilde) + dt * self.q / 2,
        )

    def test_one_flow_collision_has_exact_net_momentum_source(self) -> None:
        terms = self.flow_terms()
        f_eq = Matrix(symbols("f_eq_0:9"))
        s_plus, s_minus, dt = symbols("s_f_plus s_f_minus dt")
        p_plus, p_minus = parity_projectors(self.lattice)
        f_tilde = f_eq - dt * terms.minus / 2
        nonequilibrium = f_tilde - f_eq
        relaxation = (
            -s_plus * p_plus * nonequilibrium
            - s_minus * p_minus * nonequilibrium
        )
        explicit_source = dt * (
            (1 - s_plus / 2) * terms.plus
            + (1 - s_minus / 2) * terms.minus
        )

        self.assert_zero_matrix(
            self.vector_moment(relaxation) - s_minus * dt * self.force / 2
        )
        self.assert_zero_matrix(
            self.vector_moment(explicit_source)
            - dt * (1 - s_minus / 2) * self.force
        )
        post_collision = trt_collision(
            self.lattice,
            h_tilde=f_tilde,
            h_eq=f_eq,
            source_plus=terms.plus,
            source_minus=terms.minus,
            s_plus=s_plus,
            s_minus=s_minus,
            dt=dt,
        )
        self.assert_zero_matrix(
            self.vector_moment(post_collision - f_tilde) - dt * self.force
        )

    def test_one_scalar_collision_has_exact_net_heat_source(self) -> None:
        terms = self.scalar_terms()
        g_eq = Matrix(symbols("g_eq_0:9"))
        s_plus, s_minus, dt = symbols("s_g_plus s_g_minus dt")
        p_plus, p_minus = parity_projectors(self.lattice)
        g_tilde = g_eq - dt * terms.plus / 2
        nonequilibrium = g_tilde - g_eq
        relaxation = (
            -s_plus * p_plus * nonequilibrium
            - s_minus * p_minus * nonequilibrium
        )
        explicit_source = dt * (
            (1 - s_plus / 2) * terms.plus
            + (1 - s_minus / 2) * terms.minus
        )

        self.assert_expr_equal(sum(relaxation), s_plus * dt * self.q / 2)
        self.assert_expr_equal(
            sum(explicit_source),
            dt * (1 - s_plus / 2) * self.q,
        )
        post_collision = trt_collision(
            self.lattice,
            h_tilde=g_tilde,
            h_eq=g_eq,
            source_plus=terms.plus,
            source_minus=terms.minus,
            s_plus=s_plus,
            s_minus=s_minus,
            dt=dt,
        )
        self.assert_expr_equal(sum(post_collision - g_tilde), dt * self.q)

    def test_equal_rates_recover_bgk_componentwise(self) -> None:
        h_tilde = Matrix(symbols("h_tilde_0:9"))
        h_eq = Matrix(symbols("h_eq_0:9"))
        source_plus = Matrix(symbols("source_plus_0:9"))
        source_minus = Matrix(symbols("source_minus_0:9"))
        rate, dt = symbols("s dt")

        trt = trt_collision(
            self.lattice,
            h_tilde=h_tilde,
            h_eq=h_eq,
            source_plus=source_plus,
            source_minus=source_minus,
            s_plus=rate,
            s_minus=rate,
            dt=dt,
        )
        bgk = bgk_collision(
            h_tilde=h_tilde,
            h_eq=h_eq,
            source=source_plus + source_minus,
            rate=rate,
            dt=dt,
        )

        self.assert_zero_matrix(trt - bgk)


if __name__ == "__main__":
    unittest.main()
