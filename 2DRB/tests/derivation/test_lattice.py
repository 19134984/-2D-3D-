"""Exact algebra tests for the D2Q5 and D2Q9 lattices."""

from __future__ import annotations

import unittest

from sympy import Float, Matrix, Rational, eye, zeros


_LATTICE_IMPORT_ERROR: str | None = None

try:
    from tools.derivation.lattice import (
        d2q5,
        d2q9,
        hermite_moment,
        parity_projectors,
        raw_moment,
    )
except (ImportError, ModuleNotFoundError) as exc:
    _LATTICE_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"


class LatticeTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.assertIsNone(
            _LATTICE_IMPORT_ERROR,
            f"exact lattice API is unavailable: {_LATTICE_IMPORT_ERROR}",
        )


class LatticeDefinitionTests(LatticeTestCase):
    def test_d2q9_exact_definition(self) -> None:
        lattice = d2q9()

        self.assertEqual(
            lattice.velocities,
            (
                (0, 0),
                (1, 0),
                (0, 1),
                (-1, 0),
                (0, -1),
                (1, 1),
                (-1, 1),
                (-1, -1),
                (1, -1),
            ),
        )
        self.assertEqual(
            lattice.weights,
            (
                Rational(4, 9),
                Rational(1, 9),
                Rational(1, 9),
                Rational(1, 9),
                Rational(1, 9),
                Rational(1, 36),
                Rational(1, 36),
                Rational(1, 36),
                Rational(1, 36),
            ),
        )
        self.assertEqual(
            lattice.lambda_t,
            (
                Rational(-5, 9),
                Rational(1, 9),
                Rational(1, 9),
                Rational(1, 9),
                Rational(1, 9),
                Rational(1, 36),
                Rational(1, 36),
                Rational(1, 36),
                Rational(1, 36),
            ),
        )
        self.assertEqual(lattice.opposite, (0, 3, 4, 1, 2, 7, 8, 5, 6))
        self.assertEqual(lattice.cs2, Rational(1, 3))

    def test_d2q5_exact_definition(self) -> None:
        lattice = d2q5()

        self.assertEqual(
            lattice.velocities,
            ((0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)),
        )
        self.assertEqual(
            lattice.weights,
            (
                Rational(1, 3),
                Rational(1, 6),
                Rational(1, 6),
                Rational(1, 6),
                Rational(1, 6),
            ),
        )
        self.assertEqual(lattice.lambda_t, ())
        self.assertEqual(lattice.opposite, (0, 3, 4, 1, 2))
        self.assertEqual(lattice.cs2, Rational(1, 3))

    def test_opposite_maps_are_involutions(self) -> None:
        for lattice in (d2q5(), d2q9()):
            self.assertEqual(
                tuple(
                    lattice.opposite[lattice.opposite[index]]
                    for index in range(len(lattice.velocities))
                ),
                tuple(range(len(lattice.velocities))),
            )

    def test_lattice_constants_contain_no_floating_values(self) -> None:
        for lattice in (d2q5(), d2q9()):
            constants = (*lattice.weights, *lattice.lambda_t, lattice.cs2)
            self.assertFalse(
                any(value.atoms(Float) for value in constants),
                f"floating lattice constant found in {lattice}",
            )


class LatticeMomentTests(LatticeTestCase):
    def test_d2q9_weight_moments(self) -> None:
        lattice = d2q9()
        moment = lambda powers: raw_moment(
            lattice.weights, lattice.velocities, powers
        )

        self.assertEqual(moment((0, 0)), 1)
        self.assertEqual(moment((1, 0)), 0)
        self.assertEqual(moment((0, 1)), 0)
        self.assertEqual(moment((2, 0)), lattice.cs2)
        self.assertEqual(moment((0, 2)), lattice.cs2)
        self.assertEqual(moment((1, 1)), 0)
        self.assertEqual(moment((4, 0)), 3 * lattice.cs2**2)
        self.assertEqual(moment((0, 4)), 3 * lattice.cs2**2)
        self.assertEqual(moment((2, 2)), lattice.cs2**2)
        self.assertEqual(moment((3, 1)), 0)
        self.assertEqual(moment((1, 3)), 0)

    def test_d2q5_weight_moments_and_fourth_order_limitation(self) -> None:
        lattice = d2q5()
        moment = lambda powers: raw_moment(
            lattice.weights, lattice.velocities, powers
        )

        self.assertEqual(moment((0, 0)), 1)
        self.assertEqual(moment((1, 0)), 0)
        self.assertEqual(moment((0, 1)), 0)
        self.assertEqual(moment((2, 0)), lattice.cs2)
        self.assertEqual(moment((0, 2)), lattice.cs2)
        self.assertEqual(moment((1, 1)), 0)
        self.assertEqual(moment((4, 0)), Rational(1, 3))
        self.assertEqual(moment((0, 4)), Rational(1, 3))
        self.assertEqual(moment((2, 2)), 0)

    def test_d2q9_lambda_t_constraints(self) -> None:
        lattice = d2q9()
        moment = lambda powers: raw_moment(
            lattice.lambda_t, lattice.velocities, powers
        )

        self.assertEqual(moment((0, 0)), 0)
        self.assertEqual(moment((1, 0)), 0)
        self.assertEqual(moment((0, 1)), 0)
        self.assertEqual(moment((2, 0)), lattice.cs2)
        self.assertEqual(moment((0, 2)), lattice.cs2)
        self.assertEqual(moment((1, 1)), 0)
        self.assertEqual(moment((3, 0)), 0)
        self.assertEqual(moment((0, 3)), 0)
        self.assertEqual(
            tuple(lattice.lambda_t[lattice.opposite[i]] for i in range(9)),
            lattice.lambda_t,
        )

    def test_d2q9_lambda_t_fourth_raw_moments(self) -> None:
        lattice = d2q9()
        moment = lambda powers: raw_moment(
            lattice.lambda_t, lattice.velocities, powers
        )

        self.assertEqual(
            (
                moment((4, 0)),
                moment((0, 4)),
                moment((2, 2)),
                moment((3, 1)),
                moment((1, 3)),
            ),
            (
                Rational(1, 3),
                Rational(1, 3),
                Rational(1, 9),
                0,
                0,
            ),
        )

    def test_raw_and_hermite_moments_are_exact(self) -> None:
        lattice = d2q9()

        self.assertEqual(
            hermite_moment(
                lattice.weights,
                lattice.velocities,
                (),
                lattice.cs2,
            ),
            1,
        )
        for axes in ((0,), (1,), (0, 0), (0, 1), (1, 1), (0, 0, 0)):
            self.assertEqual(
                hermite_moment(
                    lattice.weights,
                    lattice.velocities,
                    axes,
                    lattice.cs2,
                ),
                0,
            )
        self.assertEqual(
            hermite_moment(
                lattice.lambda_t,
                lattice.velocities,
                (0, 0),
                lattice.cs2,
            ),
            lattice.cs2,
        )


class ParityProjectorTests(LatticeTestCase):
    def test_projector_algebra(self) -> None:
        lattice = d2q9()
        p_plus, p_minus = parity_projectors(lattice)

        self.assertEqual(p_plus * p_plus, p_plus)
        self.assertEqual(p_minus * p_minus, p_minus)
        self.assertEqual(p_plus * p_minus, zeros(9))
        self.assertEqual(p_minus * p_plus, zeros(9))
        self.assertEqual(p_plus + p_minus, eye(9))

    def test_projectors_match_pairwise_even_odd_definition(self) -> None:
        lattice = d2q9()
        p_plus, p_minus = parity_projectors(lattice)
        populations = Matrix(range(9))

        even = p_plus * populations
        odd = p_minus * populations
        for i, opposite in enumerate(lattice.opposite):
            self.assertEqual(even[i], (populations[i] + populations[opposite]) / 2)
            self.assertEqual(odd[i], (populations[i] - populations[opposite]) / 2)


if __name__ == "__main__":
    unittest.main()
