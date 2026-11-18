#!/usr/bin/env python
"""Exhaustive finite-field verification for the audit of
A. C. Mazumder, G. Kapetanakis, D. K. Basnet, "Normal and primitive normal
elements with prescribed traces in intermediate extensions of finite fields",
Finite Fields Appl. 110 (2026) 102745 (arXiv:2510.14840).

Usage:
    python verify_intermediate_traces.py exhaustive [small|all]
    python verify_intermediate_traces.py numeric [c31]

Outputs are written to ../Data/.
"""
import sys, os, json, math, random, time
from math import gcd
import numpy as np
from sympy import factorint, primerange

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "Data")


# ----------------------------------------------------------------------------
# Finite fields F_{p^N} with exp/log tables
# ----------------------------------------------------------------------------
def pmulmod(a, b, f, p):
    N = len(f) - 1
    r = [0] * (2 * N - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i + j] = (r[i + j] + ai * bj) % p
    for i in range(2 * N - 2, N - 1, -1):
        c = r[i]
        if c:
            for j in range(N + 1):
                r[i - N + j] = (r[i - N + j] - c * f[j]) % p
    return r[:N]


def ppow(a, e, f, p):
    N = len(f) - 1
    res = [1] + [0] * (N - 1)
    while e:
        if e & 1:
            res = pmulmod(res, a, f, p)
        a = pmulmod(a, a, f, p)
        e >>= 1
    return res


def find_prim(p, N, seed=1):
    """A primitive polynomial of degree N >= 2 over F_p (low degree first)."""
    rng = random.Random(seed)
    n = p ** N - 1
    prs = list(factorint(n).keys())
    one = [1] + [0] * (N - 1)
    xx = [0, 1] + [0] * (N - 2)
    while True:
        f = [rng.randrange(1, p)] + [rng.randrange(p) for _ in range(N - 1)] + [1]
        if ppow(xx, n, f, p) != one:
            continue
        if all(ppow(xx, n // r, f, p) != one for r in prs):
            return f


class Field:
    """F_{p^N}, N >= 2.  Elements are integer codes (base-p digits of the
    coefficient vector); the constant c in F_p has code c."""

    def __init__(self, p, N):
        assert N >= 2
        self.p, self.N = p, N
        self.size = p ** N
        self.n = n = self.size - 1
        f = find_prim(p, N)
        self.f = f
        if p == 2:
            fm = sum(c << i for i, c in enumerate(f))
            top = 1 << N
            lst = [0] * n
            c = 1
            for i in range(n):
                lst[i] = c
                c <<= 1
                if c & top:
                    c ^= fm
            assert c == 1
            exp = np.array(lst, dtype=np.int64)
        else:
            dig = np.zeros((n, N), dtype=np.int8)
            cur = [1] + [0] * (N - 1)
            for i in range(n):
                dig[i] = cur
                t = cur[-1]
                cur = [0] + cur[:-1]
                if t:
                    for j in range(N):
                        cur[j] = (cur[j] - t * f[j]) % p
            assert cur == [1] + [0] * (N - 1)
            self.dig = dig
            self.pw = np.array([p ** j for j in range(N)], dtype=np.int64)
            exp = dig.astype(np.int64) @ self.pw
        self.exp = exp
        log = np.full(self.size, -1, dtype=np.int64)
        log[exp] = np.arange(n)
        assert (log[1:] >= 0).all()
        self.log = log

    def sum_logs(self, idx_list):
        """codes of sum_j g^{idx_j} (vectorised)."""
        if self.p == 2:
            acc = np.zeros(len(idx_list[0]), dtype=np.int64)
            for ix in idx_list:
                acc ^= self.exp[ix]
            return acc
        acc = np.zeros((len(idx_list[0]), self.N), dtype=np.int16)
        for ix in idx_list:
            acc += self.dig[ix]
            acc %= self.p
        return acc.astype(np.int64) @ self.pw

    def add(self, a, b):
        if self.p == 2:
            return a ^ b
        r = 0
        m = 1
        p = self.p
        for _ in range(self.N):
            r += ((a % p + b % p) % p) * m
            a //= p
            b //= p
            m *= p
        return r

    def neg(self, a):
        if self.p == 2:
            return a
        r = 0
        m = 1
        p = self.p
        for _ in range(self.N):
            r += ((-(a % p)) % p) * m
            a //= p
            m *= p
        return r

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return int(self.exp[(self.log[a] + self.log[b]) % self.n])

    def inv(self, a):
        return int(self.exp[(-self.log[a]) % self.n])


_FIELDS = {}


def get_field(p, N):
    if (p, N) not in _FIELDS:
        _FIELDS[(p, N)] = Field(p, N)
    return _FIELDS[(p, N)]


# ---- polynomials over a Field: lists of codes, low degree first -------------
def polymul_lin(F, P, root):
    """P(x) * (x - root)."""
    res = [0] * (len(P) + 1)
    for i, c in enumerate(P):
        res[i + 1] = F.add(res[i + 1], c)
        res[i] = F.add(res[i], F.neg(F.mul(c, root)))
    return res


def poly_eval(F, P, x):
    r = 0
    for c in reversed(P):
        r = F.add(F.mul(r, x), c)
    return r


def poly_divexact(F, num, den):
    num = num[:]
    dn = len(den) - 1
    q = [0] * (len(num) - dn)
    li = F.inv(den[-1])
    for i in range(len(num) - 1, dn - 1, -1):
        c = F.mul(num[i], li)
        q[i - dn] = c
        if c:
            for j in range(dn + 1):
                num[i - dn + j] = F.add(num[i - dn + j], F.neg(F.mul(c, den[j])))
    assert all(v == 0 for v in num), "division not exact"
    return q


def poly_trim(a):
    while a and a[-1] == 0:
        a = a[:-1]
    return a


def poly_mod(F, a, b):
    a = poly_trim(a[:])
    db = len(b) - 1
    li = F.inv(b[-1])
    while len(a) - 1 >= db and a:
        c = F.mul(a[-1], li)
        sh = len(a) - 1 - db
        for j in range(db + 1):
            a[sh + j] = F.add(a[sh + j], F.neg(F.mul(c, b[j])))
        a = poly_trim(a)
    return a


def poly_gcd_deg(F, a, b):
    a, b = poly_trim(a[:]), poly_trim(b[:])
    while b:
        a, b = b, poly_mod(F, a, b)
    return len(a) - 1


def prime_power(q):
    fa = factorint(q)
    assert len(fa) == 1
    (p, e), = fa.items()
    return p, e


# ----------------------------------------------------------------------------
# Factorisation of x^m - 1 over F_q (q = p^e), with F_q realised inside K
# ----------------------------------------------------------------------------
def embedding(A, K, q, p, e):
    """Field embedding of F_q (inside A) into F_q (inside K), as a dict."""
    if e == 1:
        return {c: c for c in range(p)}
    uA = int(A.exp[A.n // (q - 1)])
    mp_ = [1]
    c = uA
    for _ in range(e):
        mp_ = polymul_lin(A, mp_, c)
        c = int(A.exp[(A.log[c] * p) % A.n])
    assert all(x < p for x in mp_)
    stepK = K.n // (q - 1)
    for j in range(1, q - 1):
        if gcd(j, q - 1) != 1:
            continue
        u = int(K.exp[(j * stepK) % K.n])
        if poly_eval(K, mp_, u) == 0:
            break
    else:
        raise RuntimeError("no embedding found")
    emb = {0: 0}
    stepA = A.n // (q - 1)
    for k in range(q - 1):
        emb[int(A.exp[(k * stepA) % A.n])] = int(K.exp[(k * j * stepK) % K.n])
    return emb


class Factors:
    """Irreducible factors of x^m - 1 over F_q, coefficients as codes of K."""

    def __init__(self, K, q, m):
        p, e = prime_power(q)
        assert p == K.p and K.n % (q - 1) == 0
        self.K, self.q, self.m, self.p = K, q, m, p
        mp = m
        while mp % p == 0:
            mp //= p
        self.mp = mp
        r = 1
        while pow(q, r, mp) != 1 % mp:
            r += 1
        NA = e * r
        if NA < 2:
            NA *= 2
        A = get_field(p, NA)
        assert A.n % mp == 0 and A.n % (q - 1) == 0
        emb = embedding(A, K, q, p, e)
        zl = A.n // mp
        seen = set()
        self.fac = []
        for i in range(mp):
            if i in seen:
                continue
            C = []
            j = i
            while j not in seen:
                seen.add(j)
                C.append(j)
                j = j * q % mp
            P = [1]
            for j in C:
                P = polymul_lin(A, P, int(A.exp[(zl * j) % A.n]))
            self.fac.append({"P": [emb[c] for c in P], "deg": len(C),
                             "ord": mp // gcd(mp, i)})

    def _split(self, d):
        assert self.m % d == 0
        s = 1
        while d % self.p == 0:
            d //= self.p
            s *= self.p
        return s, d

    def mult(self, dlist):
        """multiplicity of each irreducible factor in lcm(x^d - 1 : d in dlist)."""
        out = []
        for f in self.fac:
            mu = 0
            for d in dlist:
                s, dp = self._split(d)
                if dp % f["ord"] == 0:
                    mu = max(mu, s)
            out.append(mu)
        return out

    def Phi(self, dlist):
        r = 1
        for f, mu in zip(self.fac, self.mult(dlist)):
            if mu:
                r *= self.q ** (f["deg"] * (mu - 1)) * (self.q ** f["deg"] - 1)
        return r

    def deg_lcm(self, dlist):
        return sum(f["deg"] * mu for f, mu in zip(self.fac, self.mult(dlist)))

    def omega(self, d):
        return sum(1 for mu in self.mult([d]) if mu)

    def cofactors(self, d):
        """(x^d - 1)/P for every irreducible P | x^d - 1."""
        K = self.K
        xd = [K.neg(1)] + [0] * (d - 1) + [1]
        return [poly_divexact(K, xd, f["P"])
                for f, mu in zip(self.fac, self.mult([d])) if mu]


def L_apply(K, q, coeffs, logs):
    """codes of L_f(g^logs) for f with coefficients (codes) in F_q."""
    terms = []
    for j, c in enumerate(coeffs):
        if c:
            terms.append((int(K.log[c]) + logs * pow(q, j, K.n)) % K.n)
    return K.sum_logs(terms)


def normal_mask(K, fac, d, logs):
    """mask: g^logs (assumed in F_{q^d}) is normal over F_q in F_{q^d}."""
    mask = np.ones(len(logs), dtype=bool)
    for co in fac.cofactors(d):
        mask &= (L_apply(K, fac.q, co, logs) != 0)
    return mask


def trace(K, q, D, e, logs):
    """codes of Tr_{q^D/q^e}(g^logs)."""
    return K.sum_logs([(logs * pow(q, j * e, K.n)) % K.n for j in range(D // e)])


def independent_normal_check(K, q, m):
    """alpha normal over F_q iff gcd(sum alpha^{q^i} x^{m-1-i}, x^m-1) = 1."""
    xm = [K.neg(1)] + [0] * (m - 1) + [1]
    out = np.zeros(K.n, dtype=bool)
    for l in range(K.n):
        poly = [int(K.exp[(l * pow(q, m - 1 - i, K.n)) % K.n]) for i in range(m)]
        out[l] = poly_gcd_deg(K, xm, poly) == 0
    return out


def sub_elems(K, q, d):
    """logs of the nonzero elements of F_{q^d} inside K."""
    return np.arange(q ** d - 1, dtype=np.int64) * (K.n // (q ** d - 1))


def W_int(n):
    return 2 ** len(factorint(n))


def run_case(q, m, tuples, log):
    p, e = prime_power(q)
    t0 = time.time()
    K = get_field(p, e * m)
    n = K.n
    size = K.size
    logs = np.arange(n, dtype=np.int64)
    fac = Factors(K, q, m)
    Phim = fac.Phi([m])
    nm = normal_mask(K, fac, m, logs)
    assert int(nm.sum()) == Phim, "number of normal elements != Phi(x^m-1)"
    prim = np.gcd(logs, n) == 1
    phin = int(prim.sum())
    theta = phin / n
    Wq, Wx = W_int(n), 2 ** fac.omega(m)
    rec = {"q": q, "m": m, "normal": Phim, "primitive": phin,
           "primitive_normal": int((nm & prim).sum()), "W_qm1": Wq, "W_xm1": Wx,
           "tuples": []}
    print(f"== q={q} m={m}: #normal={Phim}=Phi(x^m-1)  #primitive={phin}  "
          f"#primitive normal={rec['primitive_normal']}", file=log)
    if size <= 20000:
        ok = bool((independent_normal_check(K, q, m) == nm).all())
        rec["independent_normality_check"] = ok
        print(f"   independent gcd test of normality agrees on all elements: {ok}", file=log)
        assert ok
    # corollary: normal elements with given trace over F_q
    T1 = trace(K, q, m, 1, logs)
    vals, cnt = np.unique(T1[nm], return_counts=True)
    okc = (len(vals) == q - 1 and set(cnt.tolist()) == {Phim // (q - 1)} and 0 not in vals.tolist())
    rec["corollary_trace_over_Fq_ok"] = bool(okc)
    print(f"   normal elements per nonzero trace over F_q: {sorted(set(cnt.tolist()))} "
          f"(expected {Phim // (q - 1)}): {okc}", file=log)
    for d in tuples:
        assert len(d) == 2 and d[0] < d[1] and d[1] % d[0] != 0
        g = gcd(d[0], d[1])
        L = d[0] * d[1] // g
        lam = d[0] + d[1] - g
        D = d[0] + d[1]
        assert fac.deg_lcm(list(d)) == lam
        Phig = fac.Phi(list(d))
        Ts = [trace(K, q, m, di, logs) for di in d]
        key = Ts[0] * size + Ts[1]
        # sub-field data
        sub = []
        for di in d:
            sl = sub_elems(K, q, di)
            sub.append({"codes": [0] + K.exp[sl].tolist(),
                        "tr": [0] + trace(K, q, di, g, sl).tolist(),
                        "nq": [False] + normal_mask(K, fac, di, sl).tolist(),
                        "logs": sl})
        by = {}
        for c, t, f in zip(sub[1]["codes"], sub[1]["tr"], sub[1]["nq"]):
            by.setdefault(t, []).append((c, f))
        adm, Sq = set(), set()
        for c, t, f in zip(sub[0]["codes"], sub[0]["tr"], sub[0]["nq"]):
            for cc, ff in by.get(t, []):
                k_ = c * size + cc
                adm.add(k_)
                if f and ff:
                    Sq.add(k_)
        assert len(adm) == q ** (D - g)
        # counts per trace tuple
        def cnts(mask):
            kv, kc = np.unique(key[mask] if mask is not None else key, return_counts=True)
            return dict(zip(kv.tolist(), kc.tolist()))
        c_all = cnts(None)
        c_all[0] = c_all.get(0, 0) + 1          # the zero element
        c_nor, c_pri, c_pn = cnts(nm), cnts(prim), cnts(nm & prim)
        reis_ok = set(c_all) == adm and set(c_all.values()) == {q ** (m - lam)}
        thm36_ok = (set(c_nor) == Sq and Phim % Phig == 0
                    and set(c_nor.values()) == {Phim // Phig})
        nonnormal = adm - Sq
        max_normal_on_nonnormal = max([c_nor.get(k_, 0) for k_ in nonnormal], default=0)
        prim_on_nonnormal = [c_pri.get(k_, 0) for k_ in nonnormal]
        pn = [c_pn.get(k_, 0) for k_ in Sq]
        main = theta * Phim / Phig
        err = theta * (Phim / q ** m) * q ** (m / 2) * (Wq - 1) * Wx
        bound_ok = all(abs(v - main) <= err * (1 + 1e-12) for v in pn)
        t = {"d": list(d), "gcd": g, "lcm": L, "lambda": lam, "D": D, "Phi_g": Phig,
             "admissible": len(adm), "normal_admissible": len(Sq),
             "nonnormal_admissible": len(nonnormal),
             "elements_per_admissible_tuple": q ** (m - lam), "reis_count_ok": bool(reis_ok),
             "normal_per_normal_tuple": Phim // Phig, "thm36_ok": bool(thm36_ok),
             "max_normal_on_nonnormal_tuple": int(max_normal_on_nonnormal),
             "min_primitive_on_nonnormal_tuple": int(min(prim_on_nonnormal)),
             "nonnormal_tuples_with_primitive": int(sum(1 for v in prim_on_nonnormal if v > 0)),
             "pn_min": int(min(pn)), "pn_max": int(max(pn)),
             "pn_tuples_with_none": int(sum(1 for v in pn if v == 0)),
             "pn_main_term": main, "pn_error_bound": err, "char_sum_bound_ok": bool(bound_ok),
             "exp_paper": m / 2 - lam - D, "exp_corrected": m / 2 - lam,
             "ineq_paper": bool(q ** (m / 2 - lam - D) >= Wq * Wx),
             "ineq_corrected": bool(q ** (m / 2 - lam) >= Wq * Wx)}
        print(f"   d={d} gcd={g} lcm={L} lambda={lam}: admissible tuples={len(adm)} "
              f"(each hit by {q ** (m - lam)} elements: {reis_ok}); normal admissible={len(Sq)}; "
              f"normal elements per normal tuple={Phim}/{Phig}={Phim // Phig}: {thm36_ok}", file=log)
        print(f"      non-normal admissible tuples={len(nonnormal)}: max #normal elements with these "
              f"traces={max_normal_on_nonnormal}; tuples still hit by primitive elements="
              f"{t['nonnormal_tuples_with_primitive']}", file=log)
        print(f"      primitive normal per normal admissible tuple: min={min(pn)} max={max(pn)} "
              f"tuples with none={t['pn_tuples_with_none']}; main term={main:.2f}, "
              f"|N-main|<=theta*Theta*q^(m/2)(W-1)W'={err:.1f}: {bound_ok}", file=log)
        assert reis_ok and thm36_ok and max_normal_on_nonnormal == 0 and bound_ok
        Sq_arr = np.fromiter(Sq, dtype=np.int64, count=len(Sq))
        member = np.isin(key, Sq_arr)
        # ---- base change of the Remark (gcd > 1): F_Q = F_{q^g}
        if g > 1:
            Q, M = q ** g, m // g
            facQ = Factors(K, Q, M)
            nmQ = normal_mask(K, facQ, M, logs)
            assert int(nmQ.sum()) == facQ.Phi([M])
            nQ = [dict(zip(s["codes"], [False] + normal_mask(K, facQ, di // g, s["logs"]).tolist()))
                  for s, di in zip(sub, d)]
            SQ = set()
            for c, tt in zip(sub[0]["codes"], sub[0]["tr"]):
                if nQ[0][c]:
                    for cc, _ in by.get(tt, []):
                        if nQ[1][cc]:
                            SQ.add(c * size + cc)
            c_Q = cnts(nmQ)
            both_t = Sq & SQ
            both_arr = np.fromiter(both_t, dtype=np.int64, count=len(both_t))
            mem_both = np.isin(key, both_arr)
            b = {"Q": Q, "M": M, "t": [di // g for di in d],
                 "normal_over_Fq": Phim, "normal_over_FQ": int(nmQ.sum()),
                 "normal_over_both": int((nm & nmQ).sum()),
                 "tuples_Fq_normal": len(Sq), "tuples_FQ_normal": len(SQ),
                 "tuples_both": len(both_t),
                 "max_FQ_normal_elements_on_Fq_only_tuples":
                     int(max([c_Q.get(k_, 0) for k_ in Sq - SQ], default=0)),
                 "FQ_route_elements": int((prim & nmQ & mem_both).sum()),
                 "FQ_route_elements_not_Fq_normal": int((prim & nmQ & mem_both & ~nm).sum())}
            t["base_change"] = b
            print(f"      base change Q=q^{g}={Q}, t={b['t']}, M={M}: normal/F_q={Phim}, "
                  f"normal/F_Q={b['normal_over_FQ']}, both={b['normal_over_both']}; "
                  f"F_q-normal tuples={len(Sq)}, F_Q-normal tuples={len(SQ)}, common={len(both_t)}; "
                  f"max #F_Q-normal elements on F_q-normal tuples that are not F_Q-normal="
                  f"{b['max_FQ_normal_elements_on_Fq_only_tuples']}; primitive F_Q-normal elements on "
                  f"common tuples={b['FQ_route_elements']}, of which not F_q-normal="
                  f"{b['FQ_route_elements_not_Fq_normal']}", file=log)
        # ---- route through F_{q^L} when L < m
        if L < m:
            QL, ML = q ** L, m // L
            facL = Factors(K, QL, ML)
            nmL = normal_mask(K, facL, ML, logs)
            assert int(nmL.sum()) == facL.Phi([ML])
            TL = trace(K, q, m, L, logs)
            bn = np.zeros(n, dtype=bool)
            nzl = TL != 0
            bn[nzl] = normal_mask(K, fac, L, K.log[TL[nzl]])
            r = {"L": L, "deg": ML,
                 "route_elements": int((prim & nmL & member).sum()),
                 "route_elements_not_Fq_normal": int((prim & nmL & member & ~nm).sum()),
                 "route_elements_beta_normal": int((prim & nmL & member & bn).sum()),
                 "route_elements_beta_normal_not_Fq_normal":
                     int((prim & nmL & member & bn & ~nm).sum())}
            t["route_lcm"] = r
            print(f"      route through F_(q^L), L={L}: primitive, normal over F_(q^L), traces a "
                  f"normal tuple: {r['route_elements']}, of which not normal over F_q: "
                  f"{r['route_elements_not_Fq_normal']}; with Tr_(m/L) normal over F_q: "
                  f"{r['route_elements_beta_normal']}, of which not normal over F_q: "
                  f"{r['route_elements_beta_normal_not_Fq_normal']}", file=log)
        rec["tuples"].append(t)
        log.flush()
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


SMALL = [
    (2, 6, [(2, 3)]), (3, 6, [(2, 3)]), (4, 6, [(2, 3)]), (5, 6, [(2, 3)]), (7, 6, [(2, 3)]),
    (8, 6, [(2, 3)]), (9, 6, [(2, 3)]),
    (2, 10, [(2, 5)]), (3, 10, [(2, 5)]), (4, 10, [(2, 5)]),
    (2, 12, [(2, 3), (3, 4), (4, 6)]), (3, 12, [(2, 3), (3, 4), (4, 6)]),
    (2, 14, [(2, 7)]), (2, 15, [(3, 5)]),
    (2, 18, [(2, 3), (2, 9), (6, 9)]), (2, 20, [(2, 5), (4, 5), (4, 10)]),
    (2, 21, [(3, 7)]), (2, 22, [(2, 11)]),
]
BIG = [
    (11, 6, [(2, 3)]), (13, 6, [(2, 3)]),
    (2, 24, [(2, 3), (3, 4), (4, 6), (3, 8), (6, 8)]),
    (4, 12, [(2, 3), (3, 4), (4, 6)]),
]


def exhaustive(which):
    cases = SMALL if which == "small" else SMALL + BIG
    recs = []
    with open(os.path.join(DATA, "exhaustive_log.txt"), "w") as log:
        for q, m, tl in cases:
            recs.append(run_case(q, m, tl, log))
            print(f"done q={q} m={m} ({recs[-1]['seconds']} s)", flush=True)
            with open(os.path.join(DATA, "exhaustive.json"), "w") as fh:
                json.dump(recs, fh, indent=1)


# ----------------------------------------------------------------------------
# Numerical part: constants, the inequalities of Theorem 5.1, corrected bounds
# ----------------------------------------------------------------------------
LN2 = math.log(2)
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19]
NUS = list(range(3, 15))


def lnC(nu):
    return sum(LN2 - math.log(p) / nu for p in primerange(2, 2 ** nu + 1))


LC = {nu: lnC(nu) for nu in NUS}


def pi_theta(limit, seg=1 << 24):
    r = int(math.isqrt(limit)) + 1
    base = np.ones(r + 1, dtype=bool)
    base[:2] = False
    for i in range(2, int(r ** 0.5) + 1):
        if base[i]:
            base[i * i::i] = False
    bp = np.nonzero(base)[0]
    cnt, th, lo = 0, 0.0, 2
    while lo <= limit:
        hi = min(lo + seg - 1, limit)
        s = np.ones(hi - lo + 1, dtype=bool)
        for p in bp:
            p = int(p)
            if p * p > hi:
                break
            st = max(p * p, ((lo + p - 1) // p) * p)
            s[st - lo::p] = False
        pr = np.nonzero(s)[0].astype(np.float64) + lo
        cnt += len(pr)
        th += float(np.log(pr).sum())
        lo = hi + 1
    return cnt, th


def tuples_coprime(k, maxm):
    """pairwise coprime 2 <= d_1 < ... < d_k with product <= maxm."""
    out = []

    def rec(pref, prod):
        if len(pref) == k:
            out.append(tuple(pref))
            return
        d = pref[-1] + 1 if pref else 2
        while prod * d ** (k - len(pref)) <= maxm:
            if all(gcd(d, e) == 1 for e in pref):
                rec(pref + [d], prod * d)
            d += 1

    rec([], 1)
    return out


def is_prime_power(n):
    return n >= 2 and len(factorint(n)) == 1


def need_lnq(m, D, k):
    """Smallest ln q such that, for some nu, with W(x^M-1) <= 2^M,
        (M/2 - D + k - 1 - M/nu) ln q >= ln C_nu + M ln 2
    holds at M = m and (the left side growing faster in M) at every multiple
    M of m.  m, D are numpy arrays; the condition is monotone in q."""
    best = np.full(len(m), np.inf)
    for nu in NUS:
        E = m / 2 - D + k - 1 - m / nu
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(E > 0, (LC[nu] + m * LN2) / E, np.inf)
        v = np.maximum(v, LN2 / (0.5 - 1 / nu))
        best = np.minimum(best, v)
    return best


def refined_ok(m, D, k, q):
    """For a fixed prime power q: the same with the Lenstra-Schoof bound
    W(x^M-1) <= 2^((M+gcd(M,q-1))/2), exact at M = m and bounded by
    2^((M+q-1)/2) at the proper multiples M >= 2m (linear in M)."""
    lq = math.log(q)
    ok0 = np.zeros(len(m), dtype=bool)
    okm = np.zeros(len(m), dtype=bool)
    B = (m + np.gcd(m, q - 1)) / 2
    for nu in NUS:
        ok0 |= (m / 2 - D + k - 1 - m / nu) * lq >= LC[nu] + B * LN2
        if (0.5 - 1 / nu) * lq - LN2 / 2 >= 0:
            okm |= (m - D + k - 1 - 2 * m / nu) * lq >= LC[nu] + (2 * m + q - 1) / 2 * LN2
    return ok0 & okm


def tail_ok(k, lnq, M0, c):
    """some nu with G(M0) >= 0 and G'(M0) >= 0, where, for m >= M0,
    G(m) = (m(1/2 - 1/P_{k-1} - 1/nu) - (k-1) sqrt(m/P_{k-2}) + k - 1) ln q - ln C_nu - c m ln 2
    is convex in m and is a lower bound for the logarithm of LHS/RHS."""
    P1 = math.prod(PRIMES[:k - 1])
    P2 = math.prod(PRIMES[:k - 2])
    for nu in NUS:
        a = 0.5 - 1 / P1 - 1 / nu
        G = (M0 * a - (k - 1) * math.sqrt(M0 / P2) + k - 1) * lnq - LC[nu] - c * M0 * LN2
        dG = (a - (k - 1) / (2 * math.sqrt(M0 * P2))) * lnq - c * LN2
        if G >= 0 and dG >= 0:
            return True
    return False


def corrected_threshold(k, M0, out):
    tl = tuples_coprime(k, M0)
    m = np.array([math.prod(d) for d in tl], dtype=np.int64)
    D = np.array([sum(d) for d in tl], dtype=np.int64)
    # (a) bound W(x^m-1) <= 2^m: monotone in q
    need = need_lnq(m, D, k)
    iw = int(np.argmax(need))
    lnq = float(need[iw])
    qa = math.ceil(math.exp(lnq) - 1e-9)
    while not tail_ok(k, math.log(qa), M0, 1.0):
        qa += 1
    # (b) Lenstra-Schoof bound with the exact gcd(m, q-1), prime powers q < qa
    fails = []
    for q in range(2, qa):
        if not is_prime_power(q):
            continue
        ok = M0 >= q and tail_ok(k, math.log(q), M0, 0.75) and bool(refined_ok(m, D, k, q).all())
        if not ok:
            fails.append(q)
    q0 = 2
    if fails:
        q0 = max(fails) + 1
        while not is_prime_power(q0):
            q0 += 1
    print(f"   k={k}: {len(tl)} tuples with product<={M0}; with W(x^m-1)<=2^m: all q>={qa} "
          f"(worst tuple {tl[iw]}, needs q>={math.exp(lnq):.2f}); with Lenstra-Schoof: "
          f"failing prime powers below {qa}: {fails if len(fails) < 40 else str(fails[:5]) + '...' + str(fails[-5:])}"
          f" -> all prime powers q>={q0}", file=out)
    out.flush()
    return {"k": k, "M0": M0, "tuples": len(tl), "q_2m": qa, "worst": list(tl[iw]),
            "worst_need": math.exp(lnq), "n_fails_refined": len(fails),
            "largest_fail_refined": max(fails) if fails else None, "q_refined": q0}


def numeric(c31):
    res = {}
    with open(os.path.join(DATA, "numeric_log.txt"), "w") as out:
        L10 = math.log(10)
        print("-- constants C_nu of Lemma 2.7 (product over all primes <= 2^nu)", file=out)
        for nu in NUS:
            print(f"   log10 C_{nu} = {LC[nu] / L10:.6f}", file=out)
        res["log10_C"] = {nu: LC[nu] / L10 for nu in NUS}
        print(f"   C_11 = {math.exp(LC[11]):.5e} (paper 4.2445e14); C_12 = {math.exp(LC[12]):.5e} "
              f"(paper 1.0573e24)", file=out)
        if c31:
            cnt, th = pi_theta(2 ** 31)
            l31 = cnt * math.log10(2) - th / (31 * L10)
            print(f"   pi(2^31) = {cnt}, theta(2^31) = {th}", file=out)
        else:
            l31 = 1553069.3804745004
        res["log10_C31"] = l31
        print(f"   log10 C_31 = {l31:.6f} (paper: 2.4015e1553069)"
              f"{'' if c31 else '  [stored value; rerun with c31 to recompute]'}", file=out)

        print("\n-- Lemma 2.9: W(t-1) < t^(0.96/log log t), 3 <= t <= 10^6", file=out)
        N = 10 ** 6
        om = np.zeros(N + 1, dtype=np.int64)
        for p in primerange(2, N + 1):
            om[p::p] += 1
        t = np.arange(3, N + 1, dtype=np.float64)
        lhs = om[2:N] * LN2
        ll = np.log(np.log(t))
        ok = (lhs < 0.96 * np.log(t) / ll)
        res["lemma26_ok_to_1e6"] = bool(ok.all())
        print(f"   holds for all t in the range: {bool(ok.all())}", file=out)

        print("\n-- Theorem 5.1(i)(a), k=4: the paper's inequality "
              "q^(m/2-8 floor(m/24)+3) >= C_11 q^(m/11) 2^m", file=out)

        def need4(m):
            E = m / 2 - 8 * (m // 24) + 3 - m / 11
            return math.exp((LC[11] + m * LN2) / E)
        ms = [m for m in range(210, 5000) if len(factorint(m)) >= 4]
        for m in ms[:6]:
            print(f"   m={m}: holds iff q >= {need4(m):.1f}", file=out)
        nf = sum(1 for m in ms if need4(m) > 1334)
        lim4 = 2 ** (1 / (0.5 - 1 / 3 - 1 / 11))
        print(f"   m -> infinity: q >= {lim4:.1f}; among the {len(ms)} values 210<=m<5000 with "
              f"omega(m)>=4 the inequality fails at q=1334 for {nf}", file=out)
        res["k4_paper"] = {"need": {m: need4(m) for m in ms[:6]}, "limit": lim4,
                           "fail_at_1334": nf, "of": len(ms)}

        print("\n-- Theorem 5.1(i)(b)-(d): the paper's inequality with nu in 3..14", file=out)
        res["k567_paper"] = {}
        for k, q, m0 in ((5, 9, 2310), (5, 8, 2310), (6, 7, 30030), (6, 5, 30030),
                         (7, 5, 510510), (7, 4, 510510)):
            good = []
            for nu in list(NUS) + [15, 16, 17, 18]:
                lc = LC[nu] if nu in LC else lnC(nu)
                kf = math.factorial(k)
                # worst case over m >= m0 of a function that is, up to the floor, linear in m
                slope = (0.5 - 2 * k / kf - 1 / nu) * math.log(q) - LN2
                E0 = m0 / 2 - 2 * k * (m0 // kf) + k - 1 - m0 / nu
                if slope > 0 and E0 * math.log(q) >= lc + m0 * LN2:
                    good.append(nu)
            res["k567_paper"][f"{k},{q}"] = good
            print(f"   k={k}, q={q}, m>={m0}: nu for which it holds at m0 with positive slope: {good}",
                  file=out)

        print("\n-- Theorem 5.1(ii), k=3: the paper's inequality q^2 >= C_31 q^(m/31) 2^m", file=out)
        res["k3_paper"] = {}
        for m in (60, 61, 62, 66, 70):
            E = 2 - m / 31
            if E > 1e-12:
                lq = (l31 + m * math.log10(2)) / E
                res["k3_paper"][m] = lq
                print(f"   m={m}: holds iff log10 q >= {lq:.2f} "
                      f"(q >= {10 ** (lq - math.floor(lq)):.4f}e{math.floor(lq)})", file=out)
            else:
                res["k3_paper"][m] = None
                print(f"   m={m}: exponent 2-m/31 = {E:.4f} <= 0, false for every q", file=out)
        tr = tuples_coprime(3, 100000)
        v = [d for d in tr if math.prod(d) >= 60 and sum(d) > math.prod(d) / 4]
        w = [d for d in tr if math.prod(d) >= 70 and d[2] > math.prod(d) / 10]
        print(f"   triples with 60<=m<=10^5 violating D<=m/4: {len(v)}; triples with m>=70 and "
              f"d_3>m/10 (contradicting the bound used): {len(w)}, e.g. {w[:3]}", file=out)
        res["k3_D_le_m4_violations"] = len(v)
        res["k3_d3_gt_m10"] = len(w)

        print("\n-- k=2: exponents m/2-lambda-D (paper) and m/2-lambda (corrected), coprime pairs",
              file=out)
        res["k2"] = []
        for d1 in range(2, 10):
            for d2 in range(d1 + 1, 31):
                if gcd(d1, d2) == 1:
                    m = d1 * d2
                    ep = m / 2 - 2 * (d1 + d2) + 1
                    ec = m / 2 - (d1 + d2) + 1
                    q0 = float(need_lnq(np.array([m]), np.array([d1 + d2]), 2)[0])
                    res["k2"].append({"d": [d1, d2], "exp_paper": ep, "exp_corrected": ec,
                                      "log10_q0": None if math.isinf(q0) else q0 / L10})
        for d1 in range(2, 10):
            row = [r for r in res["k2"] if r["d"][0] == d1]
            print(f"   d1={d1}: paper exponent >0 for d2 in "
                  f"{[r['d'][1] for r in row if r['exp_paper'] > 0]}; corrected exponent >0 for d2 in "
                  f"{[r['d'][1] for r in row if r['exp_corrected'] > 0]}", file=out)
        for d in ((3, 4), (3, 5), (5, 19), (6, 11), (7, 8), (7, 9), (8, 9)):
            r = [x for x in res["k2"] if x["d"] == list(d)][0]
            s = "none with nu<=14" if r["log10_q0"] is None else f"q >= 10^{r['log10_q0']:.2f}"
            print(f"   {d}: paper {r['exp_paper']}, corrected {r['exp_corrected']}; explicit "
                  f"corrected threshold (W(x^m-1)<=2^m, nu<=14): {s}", file=out)

        print("\n-- corrected thresholds: q^(m/2-D+k-1) >= W(q^m-1) W(x^m-1), pairwise coprime, "
              "m = d_1...d_k", file=out)
        res["corrected"] = [corrected_threshold(k, M0, out) for k, M0 in
                            ((3, 100000), (4, 200000), (5, 10 ** 6), (6, 3 * 10 ** 6),
                             (7, 3 * 10 ** 7))]
    with open(os.path.join(DATA, "numeric.json"), "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else "exhaustive"
    if mode == "exhaustive":
        exhaustive(sys.argv[2] if len(sys.argv) > 2 else "all")
    else:
        numeric(len(sys.argv) > 2 and sys.argv[2] == "c31")
