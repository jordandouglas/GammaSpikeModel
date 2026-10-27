"""Implementation checks; these do not prove a global classification."""
from fractions import Fraction
import unittest

import sympy as sp

from gamma_cumulants import UnresolvedAdmissibilityError, invert_cumulants


def cumulants(alpha, theta, beta, count=4):
    return tuple(sp.expand(theta**n * (alpha*sp.factorial(n-1) + beta*sp.rf(alpha, n))) for n in range(1, count+1))


class InversionTests(unittest.TestCase):
    def test_exact_collision_and_fourth_separation(self):
        a = (sp.Rational(90, 67), sp.Rational(67, 360), sp.Integer(3))
        b = (sp.Rational(240, 109), sp.Rational(109, 720), sp.Integer(2))
        shared = (1, sp.Rational(269, 720), sp.Rational(57241, 259200))
        self.assertEqual(cumulants(*a)[:3], shared)
        self.assertEqual(cumulants(*b)[:3], shared)
        found = invert_cumulants(*shared)
        self.assertEqual([(x.alpha, x.theta, x.beta) for x in found], [a, b])
        fourths = (sp.Rational(5417707, 31104000), sp.Rational(31505167, 186624000))
        for parameters, fourth in zip((a, b), fourths):
            self.assertEqual(cumulants(*parameters)[3], fourth)
            answer = invert_cumulants(*shared, fourth)
            self.assertEqual([(x.alpha, x.theta, x.beta) for x in answer], [parameters])
        self.assertEqual(invert_cumulants(*shared, 0), [])

    def test_independent_rational_case(self):
        self.assertEqual(cumulants(sp.Integer(1), sp.Integer(1), sp.Integer(1)), (2, 3, 8, 30))
        fibers = invert_cumulants("2", 3, Fraction(8), sp.Rational(30))
        self.assertEqual(len(fibers), 1)
        fiber = fibers[0]
        self.assertEqual((fiber.alpha, fiber.theta, fiber.beta), (1, 1, 1))
        self.assertEqual(fiber.isolating_interval, (sp.Rational(1, 2),)*2)

    def test_independent_irrational_case_and_certificate(self):
        fibers = invert_cumulants(1, "1/3", "2/9")
        self.assertEqual(len(fibers), 1)
        fiber = fibers[0]
        self.assertIsInstance(fiber.z, sp.CRootOf)
        # Independent algebraic certificate: z is the positive root of 3z²-1.
        x = sp.Symbol("x")
        self.assertEqual(sp.Poly(sp.minpoly(fiber.z, x), x).monic(), sp.Poly(x**2-sp.Rational(1,3), x))
        a, b = fiber.isolating_interval
        self.assertTrue(0 < a < b < 1)
        self.assertTrue(b <= sp.Rational(2, 3))
        self.assertTrue(a*a-a+sp.Rational(1,3) > 0)
        self.assertTrue(b*b-b+sp.Rational(1,3) > 0)
        self.assertEqual(fiber.alpha, fiber.z/(fiber.z**2-fiber.z+sp.Rational(1,3)))
        self.assertEqual(fiber.theta, fiber.z**2-fiber.z+sp.Rational(1,3))
        self.assertEqual(fiber.beta, (1-fiber.z)/fiber.z)

    def test_boundary_roots_removed(self):
        # P = z(z-1)(2z-1)(2z+1)/4; z=1/2 has d=0.
        self.assertEqual(invert_cumulants(1, "1/4", "1/8"), [])
        self.assertEqual(invert_cumulants(1, "1/4", "1/8", "1/16"), [])

    def test_repeated_root_is_one_fiber(self):
        # At c=1/4, z=1/4 is a repeated stationary root; d=1/16.
        c = sp.Rational(1, 4)
        z = sp.Rational(1, 4)
        t = 2*c*c + c*z*z-c*z-z**4+z**3
        fibers = invert_cumulants(1, c, t)
        self.assertEqual(sum(fiber.z == z for fiber in fibers), 1)
        self.assertEqual(len(fibers), 1)

    def test_nonpositive_and_invalid_inputs(self):
        self.assertEqual(invert_cumulants(0, 1, 1), [])
        self.assertEqual(invert_cumulants(1, 0, 1), [])
        for bad in (1.0, sp.Float(1), True, complex(1, 0)):
            with self.assertRaises(TypeError):
                invert_cumulants(bad, 1, 1)
        for bad in ("nan", "inf", "sqrt(2)", "1/0"):
            with self.assertRaises((TypeError, ValueError)):
                invert_cumulants(bad, 1, 1)
        with self.assertRaises(ValueError):
            invert_cumulants(1, 2)
        with self.assertRaises(ValueError):
            invert_cumulants(1, 2, 3, max_refinements=-1)
        with self.assertRaises(UnresolvedAdmissibilityError):
            invert_cumulants(1, "1/3", "2/9", max_refinements=0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
