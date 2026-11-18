#!/usr/bin/env python
"""Explicit example, in characteristic coprime to m, showing that the base
change F_q -> F_{q^d} of the Remark in Section 5 of Mazumder-Kapetanakis-Basnet
does not preserve normality over F_q.

q = 5, m = 12, (d_1, d_2) = (4, 6), d = 2, Q = 25, (t_1, t_2) = (2, 3), M = 6.
We exhibit a primitive element of F_{5^12} that is normal over F_25, whose
traces to F_{5^4} and F_{5^6} are normal over F_5 and over F_25, and that is
NOT normal over F_5.  Arithmetic is done directly in F_5[t]/(f); no tables.
Output: ../Data/base_change_example.json
"""
import os, json, random, warnings
warnings.filterwarnings("ignore")
from sympy import factorint, Poly, symbols, factor_list, rem
from verify_intermediate_traces import pmulmod, ppow, find_prim

p, m = 5, 12
f = find_prim(p, m)
n = p ** m - 1
ZERO = [0] * m
ONE = [1] + [0] * (m - 1)
T = [0, 1] + [0] * (m - 2)          # a primitive element


def add(a, b):
    return [(x + y) % p for x, y in zip(a, b)]


def smul(c, a):
    return [(c * x) % p for x in a]


def mul(a, b):
    return pmulmod(a, b, f, p)


def frob(a, j=1):
    return ppow(a, p ** (j % m), f, p)


def L(h, a):
    """L_h(a) = sum h_j a^(5^j), h a list of F_5-coefficients (low first)."""
    r = ZERO
    for j, c in enumerate(h):
        if c % p:
            r = add(r, smul(c % p, frob(a, j)))
    return r


X = symbols("x")


def factors(d):
    return [[int(c) % p for c in g.all_coeffs()[::-1]]
            for g, e in factor_list(Poly(X ** d - 1, X, modulus=p))[1]]


def cof(d, P):
    num = Poly(X ** d - 1, X, modulus=p)
    q_, r_ = divmod(num, Poly(P[::-1], X, modulus=p))
    assert r_.is_zero
    return [int(c) % p for c in q_.all_coeffs()[::-1]]


def normal_Fq(a, d):
    """a in F_{5^d}: normal over F_5."""
    return all(L(cof(d, P), a) != ZERO for P in factors(d))


def normal_FQ(a, t):
    """a in F_{25^t}, t | 6: normal over F_25 (y^t - 1 splits over F_25)."""
    roots = [ppow(T, (n // t) * i, f, p) for i in range(t)]
    for i in range(t):
        b = a
        for j in range(t):
            if j != i:                       # apply (sigma^2 - r_j)
                b = add(frob(b, 2), smul(p - 1, mul(roots[j], b)))
        if b == ZERO:
            return False
    return True


def trace(a, d):
    r = ZERO
    for j in range(m // d):
        r = add(r, frob(a, j * d))
    return r


def primitive(a):
    return a != ZERO and all(ppow(a, n // r, f, p) != ONE for r in factorint(n))


rng = random.Random(7)
while True:
    a0 = [rng.randrange(p) for _ in range(m)]
    if normal_Fq(a0, m):
        break
# an irreducible factor of x^12-1 over F_5 dividing neither x^4-1 nor x^6-1
x4, x6 = Poly(X ** 4 - 1, X, modulus=p), Poly(X ** 6 - 1, X, modulus=p)
new = [P for P in factors(12)
       if not rem(x4, Poly(P[::-1], X, modulus=p)).is_zero
       and not rem(x6, Poly(P[::-1], X, modulus=p)).is_zero]
P1 = new[0]
found = None
while found is None:
    u = [rng.randrange(p) for _ in range(m)]
    a = L(P1, L(u, a0))                      # alpha = L_{P1 u}(alpha_0)
    if primitive(a) and normal_FQ(a, 6):
        found = a
a = found
a1, a2 = trace(a, 4), trace(a, 6)
res = {
    "q": 5, "m": 12, "d": [4, 6], "Q": 25, "t": [2, 3],
    "modulus_low_degree_first": f,
    "factors_of_x^12-1_dividing_neither_x^4-1_nor_x^6-1": new,
    "alpha_coordinates": a,
    "alpha_primitive": primitive(a),
    "alpha_normal_over_F25": normal_FQ(a, 6),
    "alpha_normal_over_F5": normal_Fq(a, 12),
    "a1_in_F625": frob(a1, 4) == a1, "a2_in_F15625": frob(a2, 6) == a2,
    "a1_normal_over_F5": normal_Fq(a1, 4), "a2_normal_over_F5": normal_Fq(a2, 6),
    "a1_normal_over_F25": normal_FQ(a1, 2), "a2_normal_over_F25": normal_FQ(a2, 3),
    # Tr_{4/2}(a1) = Tr_{6/2}(a2)
    "admissible": add(a1, frob(a1, 2)) == add(a2, add(frob(a2, 2), frob(a2, 4))),
}
print(json.dumps(res, indent=1))
assert res["alpha_primitive"] and res["alpha_normal_over_F25"] and not res["alpha_normal_over_F5"]
assert all(res[k] for k in ("a1_in_F625", "a2_in_F15625", "a1_normal_over_F5", "a2_normal_over_F5",
                            "a1_normal_over_F25", "a2_normal_over_F25", "admissible"))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Data",
                       "base_change_example.json"), "w") as fh:
    json.dump(res, fh, indent=1)
