"""Parameterized D2Q5 fourth-order reference-generator tests."""

from __future__ import annotations

import inspect
import unittest
from unittest.mock import patch

import mpmath as mp
from sympy import I, Matrix, N, Rational, diag, exp, eye, simplify, symbols

from tools.derivation.lattice import d2q5
import tools.derivation.d2q5_reference as d2q5_reference_module


_D2Q5_REFERENCE_IMPORT_ERROR: str | None = None
_SERIES_IMPORT_ERROR: str | None = None

try:
    from tools.derivation.d2q5_reference import (
        amplification_matrix,
        amplification_route,
        parameterized_d2q5,
        taylor_moment_route,
    )
except (ImportError, ModuleNotFoundError) as exc:
    _D2Q5_REFERENCE_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"

try:
    from tools.derivation.series import (
        homogeneous_part,
        truncate_in_scale,
        truncate_total_degree,
    )
except (ImportError, ModuleNotFoundError) as exc:
    _SERIES_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"


class D2Q5ReferenceImportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertIsNone(
            _D2Q5_REFERENCE_IMPORT_ERROR,
            "D2Q5 reference API is unavailable: "
            f"{_D2Q5_REFERENCE_IMPORT_ERROR}",
        )

    def test_parameterized_reference_api_is_available(self) -> None:
        self.assertTrue(callable(parameterized_d2q5))


class TruncatedSeriesUtilityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertIsNone(
            _SERIES_IMPORT_ERROR,
            f"truncated-series API is unavailable: {_SERIES_IMPORT_ERROR}",
        )

    def test_total_degree_and_homogeneous_filters_are_exact(self) -> None:
        x, y = symbols("x y")
        expression = 1 + x + y + x**2 + x * y + y**3

        self.assertEqual(
            truncate_total_degree(expression, (x, y), 2),
            1 + x + y + x**2 + x * y,
        )
        self.assertEqual(homogeneous_part(expression, (x, y), 2), x**2 + x * y)

    def test_scale_truncation_drops_only_higher_orders(self) -> None:
        epsilon, coefficient = symbols("epsilon coefficient")
        expression = 1 + coefficient * epsilon + epsilon**2 + epsilon**5

        self.assertEqual(
            truncate_in_scale(expression, epsilon, 2),
            1 + coefficient * epsilon + epsilon**2,
        )


class ParameterizedD2Q5ModelTests(D2Q5ReferenceImportTests):
    def setUp(self) -> None:
        super().setUp()
        self.alpha, self.s1, self.s3, self.s4, self.lattice_speed = symbols(
            "alpha s1 s3 s4 lambda", nonzero=True
        )
        self.model = None
        model_error: str | None = None
        try:
            self.model = parameterized_d2q5(
                self.alpha,
                self.s1,
                self.s3,
                self.s4,
                lattice_speed=self.lattice_speed,
            )
        except (AttributeError, TypeError, ValueError) as exc:
            model_error = f"{type(exc).__name__}: {exc}"
        self.assertIsNone(
            model_error,
            f"parameterized D2Q5 model construction failed: {model_error}",
        )

    def test_eq79_matrix_and_eq39_collision_are_exact(self) -> None:
        expected_m = Matrix(
            [
                [1, 1, 1, 1, 1],
                [0, self.lattice_speed, 0, -self.lattice_speed, 0],
                [0, 0, self.lattice_speed, 0, -self.lattice_speed],
                [-4, 1, 1, 1, 1],
                [0, 1, -1, 1, -1],
            ]
        )
        expected_psi = Matrix(
            [
                [1, 0, 0, 0, 0],
                [0, 1 - self.s1, 0, 0, 0],
                [0, 0, 1 - self.s1, 0, 0],
                [self.alpha * self.s3, 0, 0, 1 - self.s3, 0],
                [0, 0, 0, 0, 1 - self.s4],
            ]
        )

        self.assertEqual(self.model.moment_matrix, expected_m)
        self.assertEqual(self.model.inverse_moment_matrix, expected_m.inv())
        self.assertEqual(
            self.model.moment_matrix * self.model.inverse_moment_matrix,
            eye(5),
        )
        self.assertEqual(self.model.collision_matrix, expected_psi)

    def test_parameterized_equilibrium_populations_are_reconstructed(self) -> None:
        moving_weight = (4 + self.alpha) / 20
        expected_weights = (
            (1 - self.alpha) / 5,
            moving_weight,
            moving_weight,
            moving_weight,
            moving_weight,
        )

        self.assertEqual(
            self.model.equilibrium_moments,
            Matrix([1, 0, 0, self.alpha, 0]),
        )
        self.assertEqual(self.model.equilibrium_weights, expected_weights)
        self.assertEqual(self.model.equilibrium_second_moment, (4 + self.alpha) / 10)
        self.assertEqual(
            self.model.inverse_moment_matrix * self.model.equilibrium_moments,
            Matrix(expected_weights),
        )

    def test_fixed_task1_point_is_only_alpha_minus_two_thirds(self) -> None:
        fixed = parameterized_d2q5(
            Rational(-2, 3),
            self.s1,
            self.s3,
            self.s4,
        )

        self.assertEqual(fixed.equilibrium_weights, d2q5().weights)
        self.assertEqual(fixed.equilibrium_second_moment, Rational(1, 3))


class SecondOrderRouteTests(D2Q5ReferenceImportTests):
    def setUp(self) -> None:
        super().setUp()
        (
            self.alpha,
            self.sigma1,
            self.sigma3,
            self.sigma4,
            self.lattice_speed,
            self.dt,
            self.kx,
            self.ky,
        ) = symbols(
            "alpha sigma1 sigma3 sigma4 lambda dt kx ky", nonzero=True
        )

    def test_amplification_matrix_uses_paper_positive_streaming_phase(self) -> None:
        rates = tuple(
            1 / (sigma + Rational(1, 2))
            for sigma in (self.sigma1, self.sigma3, self.sigma4)
        )
        model = parameterized_d2q5(
            self.alpha,
            *rates,
            lattice_speed=self.lattice_speed,
        )
        phase = diag(
            1,
            exp(I * self.kx * self.lattice_speed * self.dt),
            exp(I * self.ky * self.lattice_speed * self.dt),
            exp(-I * self.kx * self.lattice_speed * self.dt),
            exp(-I * self.ky * self.lattice_speed * self.dt),
        )
        expected = (
            phase
            * model.inverse_moment_matrix
            * model.collision_matrix
            * model.moment_matrix
        )

        self.assertEqual(
            amplification_matrix(
                model,
                self.kx,
                self.ky,
                time_step=self.dt,
            ),
            expected,
        )

    def test_both_routes_generate_the_diffusion_coefficient(self) -> None:
        expected = (
            (4 + self.alpha)
            * self.lattice_speed**2
            * self.dt
            * self.sigma1
            / 10
        )

        route_a = amplification_route(
            self.alpha,
            self.sigma1,
            self.sigma3,
            self.sigma4,
            lattice_speed=self.lattice_speed,
            time_step=self.dt,
            order=2,
        )
        route_b = taylor_moment_route(
            self.alpha,
            self.sigma1,
            self.sigma3,
            self.sigma4,
            lattice_speed=self.lattice_speed,
            time_step=self.dt,
            order=2,
        )

        self.assertEqual(simplify(route_a.diffusion - expected), 0)
        self.assertEqual(simplify(route_b.diffusion - expected), 0)
        self.assertIsNone(route_a.kappa40)
        self.assertIsNone(route_b.kappa40)

    def test_taylor_route_sources_contain_no_eigen_route_vocabulary(self) -> None:
        route_b_callables = {
            "taylor_moment_route": taylor_moment_route,
            "_streaming_taylor_residual": (
                d2q5_reference_module._streaming_taylor_residual
            ),
        }

        for callable_name, callable_object in route_b_callables.items():
            source = inspect.getsource(callable_object).lower()
            for forbidden in (
                "amplification_matrix",
                "amplification_route",
                "characteristic",
                "eigen",
                "gamma",
                "log(",
                "z_h",
            ):
                with self.subTest(callable=callable_name, forbidden=forbidden):
                    self.assertNotIn(forbidden, source)

    def test_taylor_route_does_not_call_amplification_implementation(self) -> None:
        with patch.object(
            d2q5_reference_module,
            "amplification_matrix",
            side_effect=AssertionError("Route A matrix must not be called"),
        ), patch.object(
            d2q5_reference_module,
            "amplification_route",
            side_effect=AssertionError("Route A recursion must not be called"),
        ):
            result = taylor_moment_route(
                Rational(-1, 2),
                Rational(1, 5),
                Rational(2, 7),
                Rational(3, 11),
                order=4,
            )

        self.assertEqual(result.diffusion, Rational(7, 100))
        self.assertEqual(result.kappa40, Rational(26483, 3850))
        self.assertEqual(result.kappa22, Rational(-25317, 1925))


class FourthOrderReferenceTests(D2Q5ReferenceImportTests):
    @classmethod
    def setUpClass(cls) -> None:
        (
            cls.alpha,
            cls.sigma1,
            cls.sigma3,
            cls.sigma4,
            cls.lattice_speed,
            cls.dt,
        ) = symbols("alpha sigma1 sigma3 sigma4 lambda dt", nonzero=True)
        cls.route_error: str | None = None
        cls.route_a = None
        cls.route_b = None
        try:
            cls.route_a = amplification_route(
                cls.alpha,
                cls.sigma1,
                cls.sigma3,
                cls.sigma4,
                lattice_speed=cls.lattice_speed,
                time_step=cls.dt,
                order=4,
            )
            cls.route_b = taylor_moment_route(
                cls.alpha,
                cls.sigma1,
                cls.sigma3,
                cls.sigma4,
                lattice_speed=cls.lattice_speed,
                time_step=cls.dt,
                order=4,
            )
        except (AttributeError, TypeError, ValueError) as exc:
            cls.route_error = f"{type(exc).__name__}: {exc}"

    def setUp(self) -> None:
        super().setUp()
        self.assertIsNone(
            self.route_error,
            f"fourth-order routes are unavailable: {self.route_error}",
        )

    def _paper_kappas(self):
        kappa40 = (
            8
            - 3 * self.alpha
            + 12 * (self.alpha + 4) * self.sigma1**2
            - 12 * (1 - self.alpha) * self.sigma1 * self.sigma3
            - 60 * self.sigma1 * self.sigma4
        )
        kappa22 = (
            -6 * (self.alpha + 4)
            + 24 * (self.alpha + 4) * self.sigma1**2
            - 24 * (1 - self.alpha) * self.sigma1 * self.sigma3
            + 120 * self.sigma1 * self.sigma4
        )
        return kappa40, kappa22

    def test_routes_match_eq41_eq42_symbolically(self) -> None:
        expected40, expected22 = self._paper_kappas()

        for route in (self.route_a, self.route_b):
            with self.subTest(route=route):
                self.assertEqual(simplify(route.kappa40 - expected40), 0)
                self.assertEqual(simplify(route.kappa22 - expected22), 0)
        self.assertEqual(simplify(self.route_a.kappa40 - self.route_b.kappa40), 0)
        self.assertEqual(simplify(self.route_a.kappa22 - self.route_b.kappa22), 0)

    def test_eq40_normalization_and_odd_order_are_exact(self) -> None:
        expected40, expected22 = self._paper_kappas()
        prefactor = (
            self.dt**3
            * self.lattice_speed**4
            * self.sigma1
            * (4 + self.alpha)
            / 1200
        )

        for route in (self.route_a, self.route_b):
            with self.subTest(route=route):
                self.assertEqual(route.cubic, 0)
                self.assertEqual(
                    simplify(route.fourth_axis - prefactor * expected40),
                    0,
                )
                self.assertEqual(
                    simplify(route.fourth_mixed - prefactor * expected22),
                    0,
                )

    def test_isotropy_is_distinct_from_complete_cancellation(self) -> None:
        expected40, expected22 = self._paper_kappas()
        isotropy_residual = simplify(expected22 - 2 * expected40)

        self.assertEqual(
            isotropy_residual,
            40 * (6 * self.sigma1 * self.sigma4 - 1),
        )
        isotropic_residual = simplify(
            isotropy_residual.subs(self.sigma4, 1 / (6 * self.sigma1))
        )
        isotropic_only_kappa40 = expected40.subs(
            {
                self.alpha: Rational(0),
                self.sigma1: Rational(1, 5),
                self.sigma3: Rational(1, 7),
            }
        ).subs(self.sigma4, Rational(5, 6))
        self.assertEqual(isotropic_residual, 0)
        self.assertNotEqual(simplify(isotropic_only_kappa40), 0)

    def test_eq55_and_intermediate_trt_point_cancel_both_coefficients(self) -> None:
        expected40, expected22 = self._paper_kappas()
        sigma3_family = (
            self.sigma1 * (self.alpha + 4) / (1 - self.alpha)
            - (2 + 3 * self.alpha)
            / (12 * self.sigma1 * (1 - self.alpha))
        )
        sigma4_family = 1 / (6 * self.sigma1)

        self.assertEqual(
            simplify(expected40.subs({self.sigma3: sigma3_family, self.sigma4: sigma4_family})),
            0,
        )
        self.assertEqual(
            simplify(expected22.subs({self.sigma3: sigma3_family, self.sigma4: sigma4_family})),
            0,
        )
        trt_point = {
            self.sigma1: 1 / (2 * 3**Rational(1, 2)),
            self.sigma3: 1 / 3**Rational(1, 2),
            self.sigma4: 1 / 3**Rational(1, 2),
        }
        self.assertEqual(simplify(expected40.subs(trt_point)), 0)
        self.assertEqual(simplify(expected22.subs(trt_point)), 0)

    def test_three_generic_rational_points_preserve_exact_route_agreement(self) -> None:
        expected40, expected22 = self._paper_kappas()
        points = (
            (Rational(-1), Rational(1, 5), Rational(2, 7), Rational(3, 11)),
            (Rational(0), Rational(1, 3), Rational(1, 4), Rational(1, 6)),
            (Rational(1, 2), Rational(2, 5), Rational(3, 8), Rational(4, 9)),
        )

        for alpha, sigma1, sigma3, sigma4 in points:
            substitutions = {
                self.alpha: alpha,
                self.sigma1: sigma1,
                self.sigma3: sigma3,
                self.sigma4: sigma4,
                self.lattice_speed: 1,
                self.dt: 1,
            }
            with self.subTest(point=(alpha, sigma1, sigma3, sigma4)):
                self.assertEqual(
                    simplify(self.route_a.kappa40.subs(substitutions)),
                    simplify(expected40.subs(substitutions)),
                )
                self.assertEqual(
                    simplify(self.route_b.kappa22.subs(substitutions)),
                    simplify(expected22.subs(substitutions)),
                )

    def test_high_precision_axis_and_diagonal_eigenvalues_confirm_q2_q4(self) -> None:
        points = (
            (Rational(-1), Rational(1, 5), Rational(2, 7), Rational(3, 11)),
            (Rational(0), Rational(1, 3), Rational(1, 4), Rational(1, 6)),
            (Rational(1, 2), Rational(2, 5), Rational(3, 8), Rational(4, 9)),
        )
        wave_number = mp.mpf("1e-6")

        with mp.workdps(80):
            for point in points:
                alpha, sigma1, sigma3, sigma4 = point
                route = amplification_route(
                    alpha,
                    sigma1,
                    sigma3,
                    sigma4,
                    order=4,
                )
                axis_decay = self._numeric_decay(point, wave_number, mp.mpf("0"))
                diagonal_decay = self._numeric_decay(
                    point,
                    wave_number,
                    wave_number,
                )
                diffusion = _mp(route.diffusion)
                axis_fourth = _mp(route.fourth_axis)
                mixed_fourth = _mp(route.fourth_mixed)
                axis_estimate = (
                    axis_decay - diffusion * wave_number**2
                ) / wave_number**4
                diagonal_estimate = (
                    diagonal_decay - 2 * diffusion * wave_number**2
                ) / wave_number**4

                with self.subTest(point=point, direction="axis"):
                    self.assertLess(abs(axis_estimate - axis_fourth), mp.mpf("1e-10"))
                with self.subTest(point=point, direction="diagonal"):
                    self.assertLess(
                        abs(diagonal_estimate - (2 * axis_fourth + mixed_fourth)),
                        mp.mpf("1e-10"),
                    )

    @staticmethod
    def _numeric_decay(point, kx, ky):
        alpha, sigma1, sigma3, sigma4 = point
        rates = tuple(
            1 / (sigma + Rational(1, 2))
            for sigma in (sigma1, sigma3, sigma4)
        )
        model = parameterized_d2q5(alpha, *rates)
        symbolic = amplification_matrix(model, kx, ky)
        numeric = mp.matrix(
            [
                [_mp(N(symbolic[row, column], 90)) for column in range(5)]
                for row in range(5)
            ]
        )
        values = mp.eig(numeric, left=False, right=False)
        hydrodynamic = min(values, key=lambda value: abs(value - 1))
        return -mp.log(hydrodynamic)


def _mp(value):
    evaluated = N(value, 90)
    real, imaginary = evaluated.as_real_imag()
    return mp.mpc(str(real), str(imaginary))


if __name__ == "__main__":
    unittest.main()
