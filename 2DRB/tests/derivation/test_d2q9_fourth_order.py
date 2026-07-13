"""Independent D2Q9 temperature modified-equation regression tests."""

from __future__ import annotations

import inspect
import unittest
from unittest.mock import patch

from sympy import (
    I,
    Matrix,
    Rational,
    diag,
    exp,
    eye,
    simplify,
    sqrt,
    symbols,
    zeros,
)

import tools.derivation.d2q9_temperature as d2q9_temperature_module
from tools.derivation.lattice import raw_moment


_D2Q9_IMPORT_ERROR: str | None = None

try:
    from tools.derivation.d2q9_temperature import (
        amplification_route,
        build_d2q9_temperature_model,
        d2q9_amplification,
        printed_dubois_coefficients,
        taylor_moment_route,
    )
    from tools.derivation.fourier_verify import high_precision_directional_fit
except (ImportError, ModuleNotFoundError) as exc:
    _D2Q9_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"


class D2Q9TemperatureImportTests(unittest.TestCase):
    def test_d2q9_fourth_order_api_is_available(self) -> None:
        self.assertIsNone(
            _D2Q9_IMPORT_ERROR,
            f"D2Q9 fourth-order API unavailable: {_D2Q9_IMPORT_ERROR}",
        )
        for entrypoint in (
            build_d2q9_temperature_model,
            d2q9_amplification,
            amplification_route,
            taylor_moment_route,
            printed_dubois_coefficients,
            high_precision_directional_fit,
        ):
            with self.subTest(entrypoint=entrypoint):
                self.assertTrue(callable(entrypoint))


class D2Q9FrozenModelTests(D2Q9TemperatureImportTests):
    @classmethod
    def setUpClass(cls) -> None:
        cls.pi, cls.chi = symbols("pi chi")
        cls.sigma_o, cls.sigma_e = symbols(
            "sigma_o sigma_e", nonzero=True
        )
        cls.dt, cls.kx, cls.ky = symbols("dt kx ky", nonzero=True)
        cls.model_error: str | None = None
        cls.model = None
        try:
            cls.model = build_d2q9_temperature_model(
                pi=cls.pi,
                chi_kappa=cls.chi,
                sigma_odd=cls.sigma_o,
                sigma_even=cls.sigma_e,
                dt=cls.dt,
                kx=cls.kx,
                ky=cls.ky,
            )
        except (AttributeError, NotImplementedError, TypeError, ValueError) as exc:
            cls.model_error = f"{type(exc).__name__}: {exc}"

    def setUp(self) -> None:
        super().setUp()
        self.assertIsNone(
            self.model_error,
            f"D2Q9 frozen model unavailable: {self.model_error}",
        )

    def assertMatrixZero(self, matrix: Matrix) -> None:
        self.assertEqual(matrix.applyfunc(simplify), zeros(*matrix.shape))

    def test_equilibrium_conservation_and_projector_algebra(self) -> None:
        model = self.model

        self.assertEqual(simplify((model.ell * model.e)[0]), 1)
        self.assertMatrixZero(model.G * model.G - model.G)
        self.assertMatrixZero(model.ell * model.G - model.ell)
        self.assertMatrixZero(model.p_plus + model.p_minus - eye(9))
        self.assertMatrixZero(model.p_plus * model.p_minus)
        self.assertMatrixZero(model.p_plus * model.p_plus - model.p_plus)
        self.assertMatrixZero(model.p_minus * model.p_minus - model.p_minus)
        self.assertMatrixZero(model.p_minus * model.G)
        self.assertMatrixZero(model.p_plus * model.G - model.G)

        for collision in (model.C0, model.C_ext, model.C_fb):
            with self.subTest(collision=collision):
                self.assertMatrixZero(model.ell * collision - model.ell)

    def test_three_raw_collision_matrices_match_the_frozen_model(self) -> None:
        model = self.model
        s_o = 1 / (self.sigma_o + Rational(1, 2))
        s_e = 1 / (self.sigma_e + Rational(1, 2))
        expected_s = s_e * model.p_plus + s_o * model.p_minus
        expected_base = eye(9) - expected_s * (eye(9) - model.G)
        wavevector = Matrix([self.kx, self.ky])
        expected_ext = (
            expected_base
            + I
            * self.dt
            * (1 - s_o / 2)
            * model.H
            * wavevector
            * model.ell
        )
        expected_fb = (
            expected_base
            - 2
            * (1 - s_o / 2)
            / (model.a + 2 * model.b * self.sigma_o)
            * model.H
            * model.J
            * (eye(9) - model.G)
        )

        self.assertMatrixZero(model.S - expected_s)
        self.assertMatrixZero(model.C0 - expected_base)
        self.assertMatrixZero(model.C_ext - expected_ext)
        self.assertMatrixZero(model.C_fb - expected_fb)

    def test_lambda_and_equilibrium_fourth_raw_moments_are_retained(self) -> None:
        model = self.model
        velocities = model.lattice.velocities

        lambda_expected = {
            (4, 0): Rational(1, 3),
            (0, 4): Rational(1, 3),
            (2, 2): Rational(1, 9),
            (3, 1): 0,
            (1, 3): 0,
        }
        equilibrium_expected = {
            (4, 0): model.a,
            (0, 4): model.a,
            (2, 2): model.lattice.cs2 * model.a,
        }

        for powers, expected in lambda_expected.items():
            with self.subTest(family="lambda", powers=powers):
                self.assertEqual(
                    simplify(raw_moment(model.lattice.lambda_t, velocities, powers)),
                    expected,
                )
        for powers, expected in equilibrium_expected.items():
            with self.subTest(family="equilibrium", powers=powers):
                self.assertEqual(
                    simplify(raw_moment(model.e, velocities, powers) - expected),
                    0,
                )

        dx, dy = symbols("dx dy")
        contraction = simplify(
            sum(
                model.e[index] * (cx * dx + cy * dy) ** 4
                for index, (cx, cy) in enumerate(velocities)
            )
        )
        self.assertEqual(
            simplify(contraction - model.a * (dx**2 + dy**2) ** 2),
            0,
        )

    def test_exact_gradient_source_retains_nonzero_third_raw_moments(self) -> None:
        model = self.model
        velocities = model.lattice.velocities
        source_x = model.H[:, 0]
        source_y = model.H[:, 1]
        cs2 = model.lattice.cs2

        expected = {
            ("x", (3, 0)): model.d,
            ("x", (1, 2)): cs2 * model.d,
            ("x", (2, 1)): 0,
            ("x", (0, 3)): 0,
            ("y", (0, 3)): model.d,
            ("y", (2, 1)): cs2 * model.d,
            ("y", (3, 0)): 0,
            ("y", (1, 2)): 0,
        }
        sources = {"x": source_x, "y": source_y}
        for (axis, powers), value in expected.items():
            with self.subTest(axis=axis, powers=powers):
                self.assertEqual(
                    simplify(raw_moment(sources[axis], velocities, powers) - value),
                    0,
                )

    def test_flux_and_odd_ghost_blocks_are_exact_projectors(self) -> None:
        model = self.model
        cs2 = model.lattice.cs2
        K = Matrix(
            [
                [weight * velocity[axis] / cs2 for axis in range(2)]
                for weight, velocity in zip(
                    model.lattice.weights,
                    model.lattice.velocities,
                    strict=True,
                )
            ]
        )
        p_flux = K * model.J
        p_odd_ghost = model.p_minus - p_flux
        hj = model.H * model.J

        self.assertMatrixZero(model.J * K - eye(2))
        self.assertMatrixZero(model.H - model.d * K)
        self.assertMatrixZero(hj * hj - model.d * hj)
        self.assertMatrixZero(p_flux * p_flux - p_flux)
        self.assertMatrixZero(p_odd_ghost * p_odd_ghost - p_odd_ghost)
        self.assertMatrixZero(model.J * p_odd_ghost)
        self.assertMatrixZero(p_flux * p_odd_ghost)
        self.assertMatrixZero(
            model.C_fb * p_odd_ghost
            - (1 - model.s_odd) * p_odd_ghost
        )

        rational_model = build_d2q9_temperature_model(
            pi=Rational(1, 9),
            chi_kappa=Rational(1, 4),
            sigma_odd=Rational(1, 5),
            sigma_even=Rational(2, 7),
        )
        self.assertEqual((rational_model.H * rational_model.J).rank(), 2)

    def test_feedback_collision_equals_the_effective_three_block_model(self) -> None:
        model = self.model
        cs2 = model.lattice.cs2
        K = Matrix(
            [
                [weight * velocity[axis] / cs2 for axis in range(2)]
                for weight, velocity in zip(
                    model.lattice.weights,
                    model.lattice.velocities,
                    strict=True,
                )
            ]
        )
        p_flux = K * model.J
        p_odd_ghost = model.p_minus - p_flux
        sigma_flux = simplify(model.b * self.sigma_o / model.a)
        s_flux = simplify(1 / (sigma_flux + Rational(1, 2)))
        block_relaxation = (
            s_flux * p_flux
            + model.s_odd * p_odd_ghost
            + model.s_even * model.p_plus
        )
        block_collision = eye(9) - block_relaxation * (eye(9) - model.G)

        self.assertEqual(
            simplify(sigma_flux - model.b * self.sigma_o / model.a),
            0,
        )
        self.assertMatrixZero(model.C_fb - block_collision)

    def test_zero_d_and_pressure_source_limits_are_distinct(self) -> None:
        zero_source = build_d2q9_temperature_model(
            pi=0,
            chi_kappa=0,
            sigma_odd=Rational(1, 5),
            sigma_even=Rational(2, 7),
            kx=Rational(1, 11),
            ky=Rational(1, 13),
        )
        self.assertEqual(zero_source.d, 0)
        self.assertMatrixZero(zero_source.H)
        self.assertMatrixZero(zero_source.C_ext - zero_source.C0)
        self.assertMatrixZero(zero_source.C_fb - zero_source.C0)
        self.assertEqual(
            simplify(zero_source.b * zero_source.sigma_odd / zero_source.a),
            zero_source.sigma_odd,
        )

        pressure_source = build_d2q9_temperature_model(
            pi=Rational(1, 9),
            chi_kappa=0,
            sigma_odd=Rational(1, 5),
            sigma_even=Rational(2, 7),
            kx=Rational(1, 11),
        )
        self.assertNotEqual(pressure_source.d, 0)
        self.assertNotEqual(pressure_source.H, zeros(9, 2))
        self.assertNotEqual(pressure_source.C_ext, pressure_source.C0)

    def test_singular_closures_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "a=0"):
            build_d2q9_temperature_model(
                pi=Rational(-1, 3),
                chi_kappa=0,
                sigma_odd=Rational(1, 5),
                sigma_even=Rational(2, 7),
            )

    def test_feedback_only_singularity_does_not_reject_external_case(self) -> None:
        parameters = {
            "pi": Rational(-2, 3),
            "chi_kappa": 0,
            "sigma_odd_ghost": Rational(1, 2),
            "sigma_even": Rational(1, 3),
            "dt": 1,
            "kx": Rational(1, 7),
            "ky": Rational(1, 11),
        }
        try:
            external = d2q9_amplification(
                case="external",
                sigma_flux=Rational(1, 2),
                **parameters,
            )
        except ValueError as exc:
            self.fail(f"external branch was rejected by feedback-only singularity: {exc}")

        self.assertEqual(external.shape, (9, 9))
        coefficients = amplification_route(
            case="external",
            pi=parameters["pi"],
            chi_kappa=parameters["chi_kappa"],
            sigma_flux=Rational(1, 2),
            sigma_odd_ghost=parameters["sigma_odd_ghost"],
            sigma_even=parameters["sigma_even"],
            order=4,
        )
        self.assertEqual(simplify(coefficients.diffusion - Rational(1, 6)), 0)
        self.assertEqual(simplify(coefficients.pde_c40 - Rational(1, 108)), 0)
        with self.assertRaisesRegex(
            ValueError,
            r"a\+2\*b\*sigma_odd=0",
        ):
            d2q9_amplification(
                case="feedback",
                sigma_flux=Rational(-1, 2),
                **parameters,
            )
        for route in (amplification_route, taylor_moment_route):
            with self.subTest(route=route.__name__):
                try:
                    route(
                        case="feedback",
                        pi=parameters["pi"],
                        chi_kappa=parameters["chi_kappa"],
                        sigma_flux=Rational(-1, 2),
                        sigma_odd_ghost=parameters["sigma_odd_ghost"],
                        sigma_even=parameters["sigma_even"],
                        order=4,
                    )
                except Exception as exc:  # assert the public failure contract
                    self.assertIsInstance(exc, ValueError)
                    self.assertRegex(str(exc), r"a\+2\*b\*sigma_odd=0")
                else:
                    self.fail(f"{route.__name__} accepted a singular feedback closure")
        with self.assertRaisesRegex(ValueError, r"a\+2\*b\*sigma_odd=0"):
            build_d2q9_temperature_model(
                pi=Rational(1, 3),
                chi_kappa=0,
                sigma_odd=-1,
                sigma_even=Rational(2, 7),
            )


class D2Q9AmplificationMatrixTests(D2Q9TemperatureImportTests):
    def assertMatrixZero(self, matrix: Matrix) -> None:
        self.assertEqual(matrix.applyfunc(simplify), zeros(*matrix.shape))

    def _streaming(self, model) -> Matrix:
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

    def test_fixed_negative_phase_is_used_in_all_three_cases(self) -> None:
        kx, ky = symbols("kx ky")
        dt = Rational(3, 2)
        sigma_o = Rational(1, 5)
        sigma_e = Rational(2, 7)
        cases = (
            ("baseline", 0, 0, sigma_o),
            ("external", Rational(1, 9), Rational(1, 4), sigma_o),
            ("feedback", Rational(1, 9), Rational(1, 4), Rational(9, 80)),
        )

        for case, pi, chi, sigma_flux in cases:
            with self.subTest(case=case):
                model = build_d2q9_temperature_model(
                    pi=pi,
                    chi_kappa=chi,
                    sigma_odd=sigma_o,
                    sigma_even=sigma_e,
                    dt=dt,
                    kx=kx,
                    ky=ky,
                )
                try:
                    amplification = d2q9_amplification(
                        case=case,
                        pi=pi,
                        chi_kappa=chi,
                        sigma_flux=sigma_flux,
                        sigma_odd_ghost=sigma_o,
                        sigma_even=sigma_e,
                        dt=dt,
                        kx=kx,
                        ky=ky,
                    )
                except NotImplementedError as exc:
                    self.fail(f"D2Q9 amplification unavailable: {exc}")
                expected_collision = {
                    "baseline": model.C0,
                    "external": model.C_ext,
                    "feedback": model.C_fb,
                }[case]
                self.assertMatrixZero(
                    amplification - self._streaming(model) * expected_collision
                )

    def test_true_baseline_rejects_pressure_or_feedback_parameters(self) -> None:
        for pi, chi in ((Rational(1, 9), 0), (0, Rational(1, 4))):
            with self.subTest(pi=pi, chi=chi):
                try:
                    with self.assertRaisesRegex(ValueError, "true baseline"):
                        d2q9_amplification(
                            case="baseline",
                            pi=pi,
                            chi_kappa=chi,
                            sigma_flux=Rational(1, 5),
                            sigma_odd_ghost=Rational(1, 5),
                            sigma_even=Rational(2, 7),
                            kx=Rational(1, 11),
                        )
                except NotImplementedError as exc:
                    self.fail(f"D2Q9 amplification unavailable: {exc}")

    def test_actual_cases_reject_mismatched_physical_flux_shifts(self) -> None:
        common = {
            "sigma_odd_ghost": Rational(1, 5),
            "sigma_even": Rational(2, 7),
            "kx": Rational(1, 11),
        }
        mismatches = (
            ("baseline", 0, 0, Rational(1, 6)),
            ("external", Rational(1, 9), Rational(1, 4), Rational(1, 6)),
            ("feedback", Rational(1, 9), Rational(1, 4), Rational(1, 5)),
        )

        for case, pi, chi, sigma_flux in mismatches:
            with self.subTest(case=case):
                with self.assertRaisesRegex(ValueError, "sigma_flux"):
                    d2q9_amplification(
                        case=case,
                        pi=pi,
                        chi_kappa=chi,
                        sigma_flux=sigma_flux,
                        **common,
                    )


class D2Q9SecondOrderRouteTests(D2Q9TemperatureImportTests):
    def _routes(self, **kwargs):
        try:
            return (
                amplification_route(order=2, **kwargs),
                taylor_moment_route(order=2, **kwargs),
            )
        except NotImplementedError as exc:
            self.fail(f"D2Q9 second-order route unavailable: {exc}")

    def test_both_routes_recover_k2_transport_in_all_actual_cases(self) -> None:
        dt = Rational(5, 4)
        sigma_o = Rational(2, 7)
        sigma_e = Rational(3, 11)
        a = Rational(4, 9)
        b = Rational(1, 4)
        points = (
            {
                "case": "baseline",
                "pi": 0,
                "chi_kappa": 0,
                "sigma_flux": sigma_o,
                "expected": Rational(1, 3) * dt * sigma_o,
            },
            {
                "case": "external",
                "pi": Rational(1, 9),
                "chi_kappa": Rational(1, 4),
                "sigma_flux": sigma_o,
                "expected": b * dt * sigma_o,
            },
            {
                "case": "feedback",
                "pi": Rational(1, 9),
                "chi_kappa": Rational(1, 4),
                "sigma_flux": b * sigma_o / a,
                "expected": b * dt * sigma_o,
            },
        )

        for point in points:
            with self.subTest(case=point["case"]):
                routes = self._routes(
                    case=point["case"],
                    pi=point["pi"],
                    chi_kappa=point["chi_kappa"],
                    sigma_flux=point["sigma_flux"],
                    sigma_odd_ghost=sigma_o,
                    sigma_even=sigma_e,
                    dt=dt,
                )
                for route in routes:
                    self.assertEqual(
                        simplify(route.diffusion - point["expected"]),
                        0,
                    )
                    self.assertIsNone(route.gamma_axis4)
                    self.assertIsNone(route.pde_c40)

    def test_external_and_feedback_match_at_k2_after_rate_specialization(self) -> None:
        pi = Rational(-1, 12)
        chi = Rational(2, 5)
        sigma_o = Rational(3, 10)
        sigma_e = Rational(4, 13)
        dt = Rational(7, 6)
        a = Rational(1, 4)
        b = Rational(1, 5)
        sigma_feedback = b * sigma_o / a

        external = self._routes(
            case="external",
            pi=pi,
            chi_kappa=chi,
            sigma_flux=sigma_o,
            sigma_odd_ghost=sigma_o,
            sigma_even=sigma_e,
            dt=dt,
        )
        feedback = self._routes(
            case="feedback",
            pi=pi,
            chi_kappa=chi,
            sigma_flux=sigma_feedback,
            sigma_odd_ghost=sigma_o,
            sigma_even=sigma_e,
            dt=dt,
        )

        for external_route, feedback_route in zip(
            external, feedback, strict=True
        ):
            self.assertEqual(
                simplify(external_route.diffusion - feedback_route.diffusion),
                0,
            )
            self.assertEqual(
                external_route.diffusion,
                simplify(b * sigma_o * dt),
            )

    def test_bgk_and_zero_diffusion_edge_limits(self) -> None:
        sigma = Rational(1, 6)
        bgk_routes = self._routes(
            case="baseline",
            pi=0,
            chi_kappa=0,
            sigma_flux=sigma,
            sigma_odd_ghost=sigma,
            sigma_even=sigma,
            dt=1,
        )
        for route in bgk_routes:
            self.assertEqual(route.diffusion, Rational(1, 18))

        zero_routes = self._routes(
            case="external",
            pi=Rational(1, 9),
            chi_kappa=1,
            sigma_flux=sigma,
            sigma_odd_ghost=sigma,
            sigma_even=Rational(2, 7),
            dt=1,
        )
        for route in zero_routes:
            self.assertEqual(route.diffusion, 0)


class D2Q9FourthOrderRouteTests(D2Q9TemperatureImportTests):
    @classmethod
    def setUpClass(cls) -> None:
        cls.points = (
            (
                "baseline_1",
                {
                    "case": "baseline",
                    "pi": 0,
                    "chi_kappa": 0,
                    "sigma_flux": Rational(2, 7),
                    "sigma_odd_ghost": Rational(2, 7),
                    "sigma_even": Rational(3, 11),
                    "dt": Rational(5, 4),
                },
            ),
            (
                "baseline_2",
                {
                    "case": "baseline",
                    "pi": 0,
                    "chi_kappa": 0,
                    "sigma_flux": Rational(1, 5),
                    "sigma_odd_ghost": Rational(1, 5),
                    "sigma_even": Rational(4, 13),
                    "dt": 1,
                },
            ),
            (
                "external_1",
                {
                    "case": "external",
                    "pi": Rational(1, 9),
                    "chi_kappa": Rational(1, 4),
                    "sigma_flux": Rational(2, 7),
                    "sigma_odd_ghost": Rational(2, 7),
                    "sigma_even": Rational(3, 11),
                    "dt": Rational(5, 4),
                },
            ),
            (
                "external_2",
                {
                    "case": "external",
                    "pi": Rational(-1, 12),
                    "chi_kappa": Rational(2, 5),
                    "sigma_flux": Rational(3, 10),
                    "sigma_odd_ghost": Rational(3, 10),
                    "sigma_even": Rational(4, 13),
                    "dt": Rational(7, 6),
                },
            ),
            (
                "feedback_1",
                {
                    "case": "feedback",
                    "pi": Rational(1, 9),
                    "chi_kappa": Rational(1, 4),
                    "sigma_flux": Rational(9, 56),
                    "sigma_odd_ghost": Rational(2, 7),
                    "sigma_even": Rational(3, 11),
                    "dt": Rational(5, 4),
                },
            ),
            (
                "feedback_2",
                {
                    "case": "feedback",
                    "pi": Rational(-1, 12),
                    "chi_kappa": Rational(2, 5),
                    "sigma_flux": Rational(6, 25),
                    "sigma_odd_ghost": Rational(3, 10),
                    "sigma_even": Rational(4, 13),
                    "dt": Rational(7, 6),
                },
            ),
        )
        cls.route_error: str | None = None
        cls.results = {}
        try:
            for name, parameters in cls.points:
                cls.results[name] = (
                    amplification_route(order=4, **parameters),
                    taylor_moment_route(order=4, **parameters),
                )
        except (AttributeError, NotImplementedError, TypeError, ValueError) as exc:
            cls.route_error = f"{type(exc).__name__}: {exc}"

    def setUp(self) -> None:
        super().setUp()
        self.assertIsNone(
            self.route_error,
            f"D2Q9 fourth-order routes unavailable: {self.route_error}",
        )

    def test_a_generic_rational_points_have_zero_exact_route_residual(self) -> None:
        fields = (
            "diffusion",
            "gamma_axis4",
            "gamma_diagonal4",
            "gamma_mixed4",
            "pde_c40",
            "pde_c22",
            "isotropy_residual",
        )
        for name, _parameters in self.points:
            route_a, route_b = self.results[name]
            for field in fields:
                with self.subTest(point=name, field=field):
                    self.assertEqual(
                        simplify(getattr(route_a, field) - getattr(route_b, field)),
                        0,
                    )
            self.assertEqual(route_a.cancellation_residual, route_b.cancellation_residual)

    def test_b_gamma_and_pde_signs_and_direction_mapping_are_distinct(self) -> None:
        parameters = dict(self.points[2][1])
        dt = parameters["dt"]
        for route in self.results["external_1"]:
            self.assertTrue(hasattr(route, "gamma_qq4"))
            self.assertTrue(hasattr(route, "pde_k_equal_diagonal"))
            self.assertEqual(route.gamma_axis3, 0)
            self.assertEqual(route.gamma_diagonal3, 0)
            self.assertEqual(
                simplify(route.gamma_axis4 + dt**3 * route.pde_c40),
                0,
            )
            self.assertEqual(
                simplify(route.gamma_mixed4 + dt**3 * route.pde_c22),
                0,
            )
            self.assertEqual(
                simplify(
                    route.gamma_qq4
                    - 2 * route.gamma_axis4
                    - route.gamma_mixed4
                ),
                0,
            )
            self.assertEqual(
                simplify(
                    route.pde_k_equal_diagonal
                    - route.pde_c40 / 2
                    - route.pde_c22 / 4
                ),
                0,
            )
            self.assertEqual(
                simplify(route.isotropy_residual - route.pde_c22 + 2 * route.pde_c40),
                0,
            )

    def test_c_generated_external_feedback_k4_residual_is_locked_exactly(self) -> None:
        expected_c40 = Rational(-1823, 465696)
        expected_c22 = Rational(-1823, 232848)
        for route_index in range(2):
            external = self.results["external_1"][route_index]
            feedback = self.results["feedback_1"][route_index]
            self.assertEqual(
                simplify(external.pde_c40 - feedback.pde_c40),
                expected_c40,
            )
            self.assertEqual(
                simplify(external.pde_c22 - feedback.pde_c22),
                expected_c22,
            )

    def test_y_route_sources_and_helpers_exclude_the_other_route(self) -> None:
        inspections = (
            (
                "route_a",
                (amplification_route, d2q9_temperature_module._amplification_direction),
                ("taylor_moment_route", "_directional_taylor_coefficients", "_directional_taylor_residual"),
            ),
            (
                "route_b",
                (
                    taylor_moment_route,
                    d2q9_temperature_module._directional_taylor_coefficients,
                    d2q9_temperature_module._directional_taylor_residual,
                ),
                (
                    "d2q9_amplification",
                    "three_block_amplification",
                    "amplification_route",
                    "_amplification_direction",
                    "characteristic",
                    "eigen",
                    "gamma",
                    "log(",
                    "z_h",
                ),
            ),
        )
        for route_name, callables, forbidden_words in inspections:
            for callable_object in callables:
                source = inspect.getsource(callable_object).lower()
                for forbidden in forbidden_words:
                    with self.subTest(
                        route=route_name,
                        callable=callable_object.__name__,
                        forbidden=forbidden,
                    ):
                        self.assertNotIn(forbidden, source)

    def test_z_bidirectional_monkeypatch_keeps_each_route_independent(self) -> None:
        parameters = dict(self.points[1][1])
        expected_a, expected_b = self.results["baseline_2"]

        with patch.object(
            d2q9_temperature_module,
            "taylor_moment_route",
            side_effect=AssertionError("Route B entry must not be called"),
        ), patch.object(
            d2q9_temperature_module,
            "_directional_taylor_residual",
            side_effect=AssertionError("Route B helper must not be called"),
        ):
            actual_a = amplification_route(order=4, **parameters)
        with patch.object(
            d2q9_temperature_module,
            "d2q9_amplification",
            side_effect=AssertionError("actual amplification must not be called"),
        ), patch.object(
            d2q9_temperature_module,
            "_amplification_direction",
            side_effect=AssertionError("Route A helper must not be called"),
        ), patch.object(
            d2q9_temperature_module,
            "amplification_route",
            side_effect=AssertionError("Route A entry must not be called"),
        ), patch.object(
            d2q9_temperature_module,
            "three_block_amplification",
            side_effect=AssertionError("three-block amplification must not be called"),
        ):
            actual_b = taylor_moment_route(order=4, **parameters)

        self.assertEqual(actual_a, expected_a)
        self.assertEqual(actual_b, expected_b)

    def test_zz_zero_source_d2q9_trt_point_cancels_both_coefficients(self) -> None:
        parameters = {
            "case": "baseline",
            "pi": 0,
            "chi_kappa": 0,
            "sigma_flux": 1 / sqrt(12),
            "sigma_odd_ghost": 1 / sqrt(12),
            "sigma_even": 1 / sqrt(3),
            "dt": 1,
        }
        for route in (
            amplification_route(order=4, **parameters),
            taylor_moment_route(order=4, **parameters),
        ):
            self.assertEqual(simplify(route.pde_c40), 0)
            self.assertEqual(simplify(route.pde_c22), 0)
            self.assertEqual(route.cancellation_residual, (0, 0))


class D2Q9QuarticConditionTests(D2Q9TemperatureImportTests):
    def test_general_zero_source_condition_is_generated_and_solved(self) -> None:
        condition_api = getattr(
            d2q9_temperature_module,
            "quartic_condition_system",
            None,
        )
        self.assertTrue(callable(condition_api))

        sigma_o, sigma_e = symbols("sigma_o sigma_e", nonzero=True)
        generated = amplification_route(
            case="baseline",
            pi=0,
            chi_kappa=0,
            sigma_flux=sigma_o,
            sigma_odd_ghost=sigma_o,
            sigma_even=sigma_e,
            order=4,
        )
        system = condition_api(generated, solve_for=(sigma_e,))
        expected_c40 = (
            sigma_o
            * (8 * sigma_e * sigma_o - 4 * sigma_o**2 - 1)
            / 36
        )
        expected_solution = (4 * sigma_o**2 + 1) / (8 * sigma_o)

        self.assertEqual(simplify(system.c40 - expected_c40), 0)
        self.assertEqual(simplify(system.c22 - 2 * expected_c40), 0)
        self.assertEqual(system.isotropy_residual, 0)
        self.assertEqual(len(system.solutions), 1)
        self.assertEqual(
            simplify(system.solutions[0][sigma_e] - expected_solution),
            0,
        )
        nontrt_sigma_o = Rational(2, 7)
        nontrt_sigma_e = simplify(
            expected_solution.subs(sigma_o, nontrt_sigma_o)
        )
        self.assertEqual(
            simplify(
                system.c40.subs(
                    {sigma_o: nontrt_sigma_o, sigma_e: nontrt_sigma_e}
                )
            ),
            0,
        )
        self.assertNotEqual(nontrt_sigma_o, 1 / sqrt(12))

    def test_external_and_feedback_generated_conditions_solve_and_substitute(self) -> None:
        condition_api = getattr(
            d2q9_temperature_module,
            "quartic_condition_system",
            None,
        )
        self.assertTrue(callable(condition_api))
        sigma_e = symbols("sigma_e", nonzero=True)
        points = (
            {
                "case": "external",
                "pi": Rational(1, 9),
                "chi_kappa": Rational(1, 4),
                "sigma_flux": Rational(2, 7),
                "sigma_odd_ghost": Rational(2, 7),
                "sigma_even": sigma_e,
            },
            {
                "case": "feedback",
                "pi": Rational(1, 9),
                "chi_kappa": Rational(1, 4),
                "sigma_flux": Rational(9, 56),
                "sigma_odd_ghost": Rational(2, 7),
                "sigma_even": sigma_e,
            },
        )
        for parameters in points:
            with self.subTest(case=parameters["case"]):
                generated = amplification_route(order=4, **parameters)
                system = condition_api(generated, solve_for=(sigma_e,))
                self.assertEqual(system.isotropy_residual, 0)
                self.assertEqual(len(system.solutions), 1)
                substitution = system.solutions[0]
                self.assertEqual(simplify(system.c40.subs(substitution)), 0)
                self.assertEqual(simplify(system.c22.subs(substitution)), 0)

    def test_degenerate_b_zero_and_a_one_branches_are_classified(self) -> None:
        sigma_e = symbols("sigma_e", nonzero=True)
        external_b_zero = amplification_route(
            case="external",
            pi=Rational(1, 9),
            chi_kappa=1,
            sigma_flux=Rational(1, 5),
            sigma_odd_ghost=Rational(1, 5),
            sigma_even=sigma_e,
            order=4,
        )
        feedback_b_zero = amplification_route(
            case="feedback",
            pi=Rational(1, 9),
            chi_kappa=1,
            sigma_flux=0,
            sigma_odd_ghost=Rational(1, 5),
            sigma_even=sigma_e,
            order=4,
        )
        external_system = d2q9_temperature_module.quartic_condition_system(
            external_b_zero,
            solve_for=(sigma_e,),
        )
        feedback_system = d2q9_temperature_module.quartic_condition_system(
            feedback_b_zero,
            solve_for=(sigma_e,),
        )
        self.assertTrue(hasattr(external_system, "status"))
        self.assertEqual(external_b_zero.pde_c40, Rational(-1, 135))
        self.assertEqual(external_system.status, "incompatible")
        self.assertEqual(feedback_b_zero.cancellation_residual, (0, 0))
        self.assertEqual(feedback_system.status, "identically_satisfied")

        sigma_o_compatible = 1 / sqrt(12)
        compatible_parameters = {
            "pi": Rational(2, 3),
            "chi_kappa": -2,
            "sigma_odd_ghost": sigma_o_compatible,
            "sigma_even": sigma_e,
            "order": 4,
        }
        compatible_external = amplification_route(
            case="external",
            sigma_flux=sigma_o_compatible,
            **compatible_parameters,
        )
        compatible_feedback = amplification_route(
            case="feedback",
            sigma_flux=sigma_o_compatible,
            **compatible_parameters,
        )
        for generated in (compatible_external, compatible_feedback):
            system = d2q9_temperature_module.quartic_condition_system(
                generated,
                solve_for=(sigma_e,),
            )
            self.assertEqual(generated.cancellation_residual, (0, 0))
            self.assertEqual(system.status, "identically_satisfied")

        incompatible_parameters = {
            "pi": Rational(2, 3),
            "chi_kappa": Rational(-1, 2),
            "sigma_odd_ghost": Rational(1, 5),
            "sigma_even": sigma_e,
            "order": 4,
        }
        incompatible_external = amplification_route(
            case="external",
            sigma_flux=Rational(1, 5),
            **incompatible_parameters,
        )
        incompatible_feedback = amplification_route(
            case="feedback",
            sigma_flux=Rational(1, 10),
            **incompatible_parameters,
        )
        for generated in (incompatible_external, incompatible_feedback):
            system = d2q9_temperature_module.quartic_condition_system(
                generated,
                solve_for=(sigma_e,),
            )
            self.assertNotEqual(generated.pde_c40, 0)
            self.assertEqual(system.status, "incompatible")


class DuboisPrintedEquationAuditTests(D2Q9TemperatureImportTests):
    def test_printed_d2q9_trt_point_has_exact_mixed_residual(self) -> None:
        sigma_odd = 1 / sqrt(12)
        sigma_even = 1 / sqrt(3)
        try:
            printed = printed_dubois_coefficients(
                sigma1=sigma_odd,
                sigma3=sigma_even,
                sigma5=sigma_odd,
                sigma7=sigma_even,
                sigma8=sigma_even,
                xi=Rational(1, 3),
                a4=Rational(1, 3),
            )
        except NotImplementedError as exc:
            self.fail(f"printed Dubois audit unavailable: {exc}")

        self.assertEqual(simplify(printed.kappa40), 0)
        self.assertEqual(simplify(printed.kappa22), 1 / sqrt(3))

    def test_printed_expressions_are_encoded_verbatim(self) -> None:
        sigma1, sigma3, sigma5, sigma7, sigma8, xi, a4 = symbols(
            "sigma1 sigma3 sigma5 sigma7 sigma8 xi a4"
        )
        try:
            printed = printed_dubois_coefficients(
                sigma1=sigma1,
                sigma3=sigma3,
                sigma5=sigma5,
                sigma7=sigma7,
                sigma8=sigma8,
                xi=xi,
                a4=a4,
            )
        except NotImplementedError as exc:
            self.fail(f"printed Dubois audit unavailable: {exc}")
        expected40 = sigma1 * (
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
        expected22 = (
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
        self.assertEqual(simplify(printed.kappa40 - expected40), 0)
        self.assertEqual(simplify(printed.kappa22 - expected22), 0)


class D2Q9HighPrecisionFourierTests(D2Q9TemperatureImportTests):
    def test_shrinking_wave_numbers_confirm_q2_sign_and_q4_order(self) -> None:
        try:
            fit = high_precision_directional_fit(
                case="baseline",
                pi=0,
                chi_kappa=0,
                sigma_flux=Rational(1, 5),
                sigma_odd_ghost=Rational(1, 5),
                sigma_even=Rational(4, 13),
                dt=1,
                precision=80,
                wave_numbers=(Rational(1, 50), Rational(1, 100), Rational(1, 200), Rational(1, 400)),
            )
        except NotImplementedError as exc:
            self.fail(f"high-precision Fourier fit unavailable: {exc}")

        self.assertEqual(fit.precision, 80)
        self.assertTrue(all(value.real > 0 for value in fit.axis_gamma))
        self.assertTrue(all(abs(value.imag) < 1e-60 for value in fit.axis_gamma))
        self.assertLess(abs(fit.axis_orders[-1] - 4), Rational(1, 20))
        self.assertLess(abs(fit.equal_diagonal_orders[-1] - 4), Rational(1, 20))
        exact_gamma4 = float(Rational(217, 58500))
        self.assertLess(
            abs(fit.axis_quartic[-1] - exact_gamma4),
            Rational(1, 10) ** 8,
        )
        self.assertLess(
            abs(fit.equal_diagonal_quartic[-1] - exact_gamma4),
            Rational(1, 10) ** 8,
        )
        self.assertLess(
            abs(fit.equal_diagonal_quartic[-1] - fit.axis_quartic[-1]),
            Rational(1, 10) ** 8,
        )

    def test_cancellation_changes_the_fitted_residual_from_q4_to_q6(self) -> None:
        try:
            fit = high_precision_directional_fit(
                case="baseline",
                pi=0,
                chi_kappa=0,
                sigma_flux=1 / sqrt(12),
                sigma_odd_ghost=1 / sqrt(12),
                sigma_even=1 / sqrt(3),
                dt=1,
                precision=80,
                wave_numbers=(Rational(1, 25), Rational(1, 50), Rational(1, 100), Rational(1, 200)),
            )
        except NotImplementedError as exc:
            self.fail(f"high-precision Fourier fit unavailable: {exc}")

        self.assertLess(abs(fit.axis_orders[-1] - 6), Rational(1, 10))
        self.assertLess(
            abs(fit.equal_diagonal_orders[-1] - 6),
            Rational(1, 10),
        )
        self.assertLess(
            abs(fit.axis_quartic[-1]),
            abs(fit.axis_quartic[0]) / 50,
        )


if __name__ == "__main__":
    unittest.main()
