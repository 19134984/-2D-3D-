"""Small exact helpers for truncated formal polynomial series."""

from __future__ import annotations

from collections.abc import Sequence

from sympy import Expr, Poly, expand, series, sympify


def truncate_total_degree(
    expression: Expr,
    variables: Sequence[Expr],
    degree: int,
) -> Expr:
    """Drop polynomial monomials whose total selected-variable degree is larger."""

    polynomial = Poly(expand(sympify(expression)), *variables)
    return expand(
        sum(
            coefficient
            * _monomial(variables, powers)
            for powers, coefficient in polynomial.terms()
            if sum(powers) <= degree
        )
    )


def homogeneous_part(
    expression: Expr,
    variables: Sequence[Expr],
    degree: int,
) -> Expr:
    """Return the exact homogeneous component of the requested total degree."""

    polynomial = Poly(expand(sympify(expression)), *variables)
    return expand(
        sum(
            coefficient
            * _monomial(variables, powers)
            for powers, coefficient in polynomial.terms()
            if sum(powers) == degree
        )
    )


def truncate_in_scale(expression: Expr, scale: Expr, degree: int) -> Expr:
    """Return a formal series in ``scale`` through the requested degree."""

    return expand(series(sympify(expression), scale, 0, degree + 1).removeO())


def _monomial(variables: Sequence[Expr], powers: tuple[int, ...]) -> Expr:
    result = sympify(1)
    for variable, power in zip(variables, powers, strict=True):
        result *= variable**power
    return result
