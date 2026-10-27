"""Exact rational-input inversion for the frozen gamma cumulant formulas.

This is a computational helper, not a proof of a global fiber classification.
No estimation from noisy or floating-point cumulants is implemented.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from typing import Any

import sympy as sp


Z = sp.Symbol("z")


class UnresolvedAdmissibilityError(RuntimeError):
    """The configured exact interval refinement budget was exhausted."""


@dataclass(frozen=True)
class Fiber:
    """One exact positive parameter triple and its normalized root certificate.

    The rational interval is a singleton for a rational root; otherwise it
    isolates exactly one root of the filtered squarefree polynomial.
    """

    z: sp.Expr
    isolating_interval: tuple[sp.Rational, sp.Rational]
    alpha: sp.Expr
    theta: sp.Expr
    beta: sp.Expr


def _rational(value: Any) -> sp.Rational:
    if isinstance(value, bool):
        raise TypeError("Boolean values are not cumulant inputs")
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    if isinstance(value, (int, sp.Rational)):
        return sp.Rational(value)
    if isinstance(value, str):
        try:
            result = sp.Rational(value.strip())
        except (TypeError, ValueError, ZeroDivisionError) as error:
            raise ValueError(f"Not an exact rational input: {value!r}") from error
        if isinstance(result, sp.Rational):
            return result
        raise ValueError(f"Not a finite exact rational input: {value!r}")
    raise TypeError(
        "Cumulants must be str, int, fractions.Fraction, or SymPy Rational; "
        "floating-point and approximate inputs are unsupported"
    )


def _d_bounds(a: sp.Rational, b: sp.Rational, c: sp.Rational):
    """Exact extrema of z**2-z+c on the closed rational interval."""
    endpoints = (a * a - a + c, b * b - b + c)
    upper = max(endpoints)
    lower = c - sp.Rational(1, 4) if a <= sp.Rational(1, 2) <= b else min(endpoints)
    return lower, upper


def _admissible_interval(poly, root, interval, c, max_refinements):
    if isinstance(root, sp.Rational):
        if 0 < root < 1 and root * root - root + c > 0:
            return (root, root)
        return None

    a, b = interval
    for iteration in range(max_refinements + 1):
        if b <= 0 or a >= 1:
            return None
        lower, upper = _d_bounds(a, b, c)
        if upper < 0:
            return None
        if a > 0 and b < 1 and lower > 0:
            return (a, b)
        if iteration == max_refinements:
            break
        a, b = poly.refine_root(a, b, steps=1)
    raise UnresolvedAdmissibilityError(
        f"Could not certify admissibility after {max_refinements} refinements; "
        f"root={root}, interval=({a}, {b}). Increase max_refinements."
    )


def invert_cumulants(*cumulants: Any, max_refinements: int = 256) -> list[Fiber]:
    """Return every positive triple for three or four EXACT rational cumulants.

    Inputs may be positional, or a single list/tuple of length three or four.
    Invalid exact input types raise; nonpositive first/second cumulants and
    incompatible exact cumulants yield an empty list. Floats are rejected.
    A fourth cumulant is imposed by polynomial gcd over QQ, never tolerance.
    Multiplicities do not duplicate returned fibers. Any exhausted interval
    budget raises rather than returning an incomplete answer.
    """
    if len(cumulants) == 1 and isinstance(cumulants[0], (list, tuple)):
        cumulants = tuple(cumulants[0])
    if len(cumulants) not in (3, 4):
        raise ValueError("Exactly three or four rational cumulants are required")
    if isinstance(max_refinements, bool) or not isinstance(max_refinements, int) or max_refinements < 0:
        raise ValueError("max_refinements must be a nonnegative integer")
    values = tuple(_rational(value) for value in cumulants)
    mu, k2, k3 = values[:3]
    if mu <= 0 or k2 <= 0:
        return []
    c, t = k2 / mu**2, k3 / mu**3
    d = Z**2 - Z + c
    f = 2*c**2 + c*Z**2 - c*Z - Z**4 + Z**3
    poly = sp.Poly(Z**4 - Z**3 - c*Z**2 + c*Z + t - 2*c**2, Z, domain=sp.QQ)
    if len(values) == 4:
        j = (Z + 3*d)*f - 2*Z**2*d**2
        fourth = sp.Poly(j - values[3]/mu**4, Z, domain=sp.QQ)
        poly = sp.gcd(poly, fourth)
    if poly.degree() == 0:
        return []
    poly = poly.sqf_part()
    boundary = sp.Poly(Z * (Z - 1) * d, Z, domain=sp.QQ)
    poly = poly.exquo(sp.gcd(poly, boundary)).monic()
    if poly.degree() == 0:
        return []

    # Both lists are ordered increasingly. intervals() uses exact rational
    # Sturm isolation; real_roots(radicals=False) preserves exact CRootOfs.
    intervals = poly.intervals()
    roots = poly.real_roots(radicals=False)
    if len(roots) != len(intervals) or any(multiplicity != 1 for _, multiplicity in intervals):
        raise RuntimeError("Unexpected non-squarefree real-root isolation result")
    result = []
    for root, (interval, _) in zip(roots, intervals):
        certificate = _admissible_interval(poly, root, interval, c, max_refinements)
        if certificate is None:
            continue
        root_d = root**2 - root + c
        result.append(Fiber(
            z=root,
            isolating_interval=certificate,
            alpha=root/root_d,
            theta=mu*root_d,
            beta=(1-root)/root,
        ))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cumulants", nargs="+", help="three or four exact rational strings, e.g. 1 269/720 57241/259200")
    parser.add_argument("--max-refinements", type=int, default=256)
    parser.add_argument("--json", action="store_true", help="serialize exact expressions and interval endpoints as strings")
    args = parser.parse_args()
    try:
        fibers = invert_cumulants(*args.cumulants, max_refinements=args.max_refinements)
    except (TypeError, ValueError, UnresolvedAdmissibilityError) as error:
        parser.error(str(error))
    if args.json:
        print(json.dumps({
            "positive_fibers": len(fibers),
            "fibers": [{
                "z": str(fiber.z),
                "isolating_interval": [str(endpoint) for endpoint in fiber.isolating_interval],
                "alpha": str(fiber.alpha),
                "theta": str(fiber.theta),
                "beta": str(fiber.beta),
            } for fiber in fibers],
        }, indent=2))
        return
    print(f"positive fibers: {len(fibers)}")
    for index, fiber in enumerate(fibers, 1):
        print(f"[{index}] z = {fiber.z}")
        print(f"    isolating interval = {fiber.isolating_interval}")
        print(f"    alpha = {fiber.alpha}")
        print(f"    theta = {fiber.theta}")
        print(f"    beta = {fiber.beta}")


if __name__ == "__main__":
    main()
