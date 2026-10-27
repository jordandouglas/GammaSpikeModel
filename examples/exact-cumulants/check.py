"""Exact symbolic identities and a rational collision, not a proof by sampling."""

from fractions import Fraction
import json
from pathlib import Path
import sympy as sp


def require(condition, message):
    if not condition:
        raise ValueError(message)


def equal(left, right, message):
    require(sp.cancel(left - right) == 0, message)


def rational(value):
    q = Fraction(value)
    return sp.Rational(q.numerator, q.denominator)


def main():
    data = json.loads(Path(__file__).with_name("witness.json").read_text())
    require(data["schema_version"] == 1, "schema")
    require(data["parameter_order"] == ["alpha", "theta", "beta"], "parameter order")
    require(data["law"] == "B~Poisson(beta); D=sum_{i=1}^{B+1} Gamma_i(alpha,theta)", "model")
    alpha, theta, beta, s = sp.symbols("alpha theta beta s")
    z, w, c, mu, S, T, q = sp.symbols("z w c mu S T q")
    K = -alpha * sp.log(1 - theta * s) + beta * ((1 - theta * s) ** (-alpha) - 1)
    cumulants = [
        sp.factor(sp.diff(K, s, n).subs(s, 0)) for n in range(1, 5)
    ]
    for n, value in enumerate(cumulants, 1):
        equal(value, theta**n * (alpha * sp.factorial(n - 1) + beta * sp.rf(alpha, n)), "CGF derivative")
    h = alpha * theta
    for n in range(1, 4):
        equal(cumulants[n], (h + n * theta) * cumulants[n - 1] - h**2 * sp.factorial(n - 1) * theta**(n - 1), "cumulant recurrence")
    d = z*z - z + c
    f = 2*c*c + c*z*z - c*z - z**4 + z**3
    j = sp.expand((z + 3*d) * f - 2*z*z*d*d)
    sub = {alpha: z/d, theta: mu*d, beta: (1-z)/z}
    for n, target in enumerate([mu, mu**2*c, mu**3*f, mu**4*j]):
        equal(cumulants[n].subs(sub), target, "bijective normalization identity")
    g = sp.diff(f, z)
    equal(g, -4*z**3 + 3*z*z + 2*c*z - c, "fiber derivative")
    equal(sp.diff(g, z), 6*z*(1-2*z) + 2*c, "critical point monotonicity identity")
    equal(g.subs(z, sp.Rational(1, 2)), sp.Rational(1, 4), "central derivative")
    equal(g.subs(z, 1), c-1, "right derivative")
    equal(g.subs(c, z*(1-z)), z*(6*z*(1-z)-1), "admissibility boundary identity")
    equal(f-2*c*c, z*(1-z)*(z*z-c), "baseline level separation")
    equal(g.subs({c: q*q, z: q}), 2*q*q*(1-q), "critical roots versus sqrt(c)")
    equal(rational(data["conditional_injectivity_threshold_c"]), sp.Rational(1, 6), "threshold")

    f_div = sp.cancel((f-f.subs(z, w))/(z-w))
    relation = (S-1)*(c-S*S)+(2*S-1)*T
    equal(f_div, relation.subs({S: z+w, T: z*w}), "three-cumulant divided difference")
    E = (
        2*S**3*(S-1)**3 - S*(S-1)*(8*S*S-15*S+6)*T
        + (3*S-2)*(S-2)*T*T
    )
    c_collision = S*S-(2*S-1)*T/(S-1)
    j_div = sp.cancel((j-j.subs(z, w))/(z-w))
    equal(
        j_div.subs(c, c_collision.subs({S: z+w, T: z*w})),
        (-E/(S-1)).subs({S: z+w, T: z*w}),
        "four-cumulant exclusion identity",
    )
    equal(E.subs(T, 0), 2*S**3*(S-1)**3, "first endpoint")
    equal(E.subs(T, S*S/4), S**3*(S-2)*(3*(S-1)**2+1)/16, "second endpoint")
    equal(E.subs(T, S-1), 2*(S-2)*(S-1)**5, "upper-sum lower endpoint")
    derivative = sp.diff(E, T)
    equal(derivative.subs(T, 0), -S*(S-1)*(8*S*S-15*S+6), "lower-sum derivative")
    equal(derivative.subs(T, S-1), -(S-1)*(8*S**3-21*S*S+22*S-8), "upper-sum derivative")
    equal((8*S**3-21*S*S+22*S-8).subs(S, q+1), 8*q**3+3*q*q+4*q+1, "positive cubic certificate")
    require(all(x > 0 for x in [8, 3, 4, 1]), "positive coefficients")
    poly = 8*S*S-15*S+6
    equal(poly.subs(S, sp.Rational(2, 3)), -sp.Rational(4, 9), "convex endpoint 2/3")
    equal(poly.subs(S, 1), -1, "convex endpoint 1")

    computed = []
    require(len(data["triples"]) == len(data["cumulants"]) == 2, "witness cardinality")
    for triple, expected in zip(data["triples"], data["cumulants"]):
        require(len(triple) == 3 and len(expected) == 4, "witness dimensions")
        aa, tt, bb = map(rational, triple)
        require(aa > 1 and tt > 0 and bb > 0, "witness positive domain")
        values = [sp.factor(k.subs({alpha: aa, theta: tt, beta: bb})) for k in cumulants]
        require(values == list(map(rational, expected)), "independent CGF witness evaluation")
        computed.append(values)
    require(data["triples"][0] != data["triples"][1], "distinct triples")
    require(computed[0][:3] == computed[1][:3], "three-cumulant collision")
    equal(computed[0][3]-computed[1][3], rational(data["fourth_difference"]), "fourth difference")
    require(rational(data["fourth_difference"]) > 0, "nonzero fourth difference")
    require(data["global_minimal_cumulant_order"] == 4, "minimal order constant")
    cc, tt = computed[0][1], computed[0][2]
    P = z**4-z**3-cc*z*z+cc*z+tt-2*cc*cc
    equal(P, (3*z-1)*(4*z-1)*(60*z*z-25*z-42)/720, "complete witness quartic factorization")
    print(json.dumps({
        "status": "exact-symbolic-identities-and-witness-passed",
        "positive_collision_triples": 2,
        "full_domain_root_count_and_sign_proof": "proof.md; requires mathematical review",
        "empirical_inference_or_noise_validation": False,
    }, indent=2))


if __name__ == "__main__":
    main()
