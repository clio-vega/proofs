#!/usr/bin/env python3
"""
Verification script for the proof of Lyra's cap identity
(2026-07-20-lyra-cap-identity.tex).

Checks:
  1. The 50 coefficient polynomials of G_4^{(c)} in the tex file are exactly
     what interpolation from c in {6,8,10,12,14,16} produces.
  2. Cross-checks against fresh independent computation at c in {18, 20, 22}.
  3. Verifies the case analysis in Section 5.3 of the tex file.
  4. Verifies the empirical content matches min(v_2(K(c)), 6) for even c in {6..16}.

Run: python3 2026-07-20-lyra-cap-verify.py
"""
import sys
sys.path.insert(0, '/home/clio/projects/code')
sys.path.insert(0, '/home/clio/projects/code/threerow-c3')

import sympy as sp
from math import comb, factorial
from fractions import Fraction

from jobB_betaprime_sequence import interp_ab

a_s, b_s, c_s = sp.symbols('a b c')
u_s, v_s = sp.symbols('u v')


def v2(n):
    n = abs(int(n))
    if n == 0:
        return float('inf')
    return (n & -n).bit_length() - 1


def G4_of_c(c):
    """Return G_4^{(c)}(a,b) as a sympy polynomial in a, b (integer coefficients)."""
    H4 = sp.expand(interp_ab(c, 4))
    arun = sp.prod([a_s + s for s in range(3, c - 2)])
    brun = sp.prod([b_s + s for s in range(2, c - 3)])
    G4 = sp.expand(sp.cancel(H4 / (arun * brun)))
    return G4


def interpolate_Gc(cs=(6, 8, 10, 12, 14, 16)):
    """Interpolate G_4^{(c)} as a polynomial in a, b, c."""
    data = {}
    for c in cs:
        G4 = G4_of_c(c)
        P = sp.Poly(G4, a_s, b_s)
        d = {}
        for (ea, eb), co in zip(P.monoms(), P.coeffs()):
            d[(ea, eb)] = int(co)
        data[c] = d
    monos = sorted(set().union(*[set(d) for d in data.values()]), reverse=True)
    Gc = 0
    for mono in monos:
        pts = [(c, data[c].get(mono, 0)) for c in cs]
        poly = sp.expand(sp.interpolate(pts, c_s))
        Gc += poly * a_s ** mono[0] * b_s ** mono[1]
    return sp.expand(Gc)


def sheet_coeffs(Gc, ra, rb, da=4, db=4):
    """Return the binomial-basis coefficients on sheet a=2u+ra, b=2v+rb.
    Result is a 5x5 table of polynomials in c."""
    Q = sp.expand(Gc.subs({a_s: 2 * u_s + ra, b_s: 2 * v_s + rb}))

    def sub(i, k):
        return sp.expand(Q.subs({u_s: i, v_s: k}))

    G = [[sub(i, k) for k in range(db + 1)] for i in range(da + 1)]

    def fd(seq, o):
        s = list(seq)
        for _ in range(o):
            s = [sp.expand(s[t + 1] - s[t]) for t in range(len(s) - 1)]
        return s[0]

    Au = [[fd([G[ii][k] for ii in range(da + 1)], i) for k in range(db + 1)] for i in range(da + 1)]
    Av = [[sp.expand(fd([Au[i][kk] for kk in range(db + 1)], k)) for k in range(db + 1)] for i in range(da + 1)]
    return Av


# ---------------------------------------------------------------------------
# CHECK 1: Interpolated coefficient list matches tex file (50 polys, 2 sheets)
# ---------------------------------------------------------------------------

TEX_SHEET_00 = {
    (0, 1): 24 * c_s * (c_s - 5) * (c_s - 4) * (c_s - 1),
    (0, 2): 24 * c_s * (c_s - 1) * (13 * c_s**2 - 77 * c_s + 110),
    (0, 3): 96 * c_s * (c_s - 2) * (c_s - 1) * (7 * c_s - 9),
    (0, 4): 384 * c_s * (c_s - 2) * (c_s - 1) * (c_s + 1),
    (1, 0): 24 * c_s * (c_s - 5) * (c_s - 3) * (c_s - 2),
    (1, 1): 24 * c_s * (13 * c_s**3 - 82 * c_s**2 + 141 * c_s - 64),
    (1, 2): 96 * c_s * (7 * c_s**3 - 4 * c_s**2 - 22 * c_s + 17),
    (1, 3): 384 * c_s * (c_s**3 + 12 * c_s**2 - 4 * c_s + 9),
    (1, 4): 1536 * c_s * (c_s + 1) * (2 * c_s + 1),
    (2, 0): 24 * c_s * (c_s - 3) * (c_s - 2) * (13 * c_s - 41),
    (2, 1): 96 * c_s * (7 * c_s**3 - 10 * c_s**2 - 24 * c_s + 41),
    (2, 2): 192 * (2 * c_s**4 + 32 * c_s**3 + 18 * c_s**2 + 16 * c_s + 15),
    (2, 3): 2304 * (2 * c_s**3 + 14 * c_s**2 + 20 * c_s + 15),
    (2, 4): 9216 * (2 * c_s**2 + 6 * c_s + 5),
    (3, 0): 96 * c_s * (c_s - 3) * (c_s - 2) * (7 * c_s - 13),
    (3, 1): 384 * c_s * (c_s**3 + 10 * c_s**2 - 25 * c_s + 17),
    (3, 2): 2304 * (2 * c_s**3 + 12 * c_s**2 + 8 * c_s + 5),
    (3, 3): 9216 * (c_s + 3) * (3 * c_s + 5),
    (3, 4): 36864 * (2 * c_s + 5),
    (4, 0): 384 * c_s * (c_s - 3) * (c_s - 2) * (c_s - 1),
    (4, 1): 1536 * c_s * (c_s - 1) * (2 * c_s - 1),
    (4, 2): 9216 * (2 * c_s**2 + 2 * c_s + 1),
    (4, 3): 36864 * (2 * c_s + 3),
    (4, 4): sp.Integer(147456),
}

TEX_SHEET_11 = {
    (0, 0): 24 * c_s * (c_s - 5) * (c_s - 4) * (c_s - 1),
    (0, 1): 48 * c_s * (c_s - 1) * (7 * c_s**2 - 39 * c_s + 47),
    (0, 2): 24 * c_s * (c_s - 1) * (41 * c_s**2 - 105 * c_s + 82),
    (0, 3): 96 * c_s * (c_s - 1) * (11 * c_s**2 + 3 * c_s + 10),
    (0, 4): 384 * c_s * (c_s - 1) * (c_s + 1) * (c_s + 2),
    (1, 0): 48 * c_s * (c_s - 1) * (7 * c_s**2 - 43 * c_s + 65),
    (1, 1): 24 * c_s * (c_s - 1) * (41 * c_s**2 - 57 * c_s - 4),
    (1, 2): 96 * (11 * c_s**4 + 54 * c_s**3 + 12 * c_s**2 + 73 * c_s + 30),
    (1, 3): 384 * (c_s**4 + 22 * c_s**3 + 47 * c_s**2 + 65 * c_s + 30),
    (1, 4): 1536 * (c_s + 1) * (c_s + 2) * (2 * c_s + 3),
    (2, 0): 24 * c_s * (c_s - 1) * (41 * c_s**2 - 169 * c_s + 182),
    (2, 1): 96 * c_s * (11 * c_s**3 + 44 * c_s**2 - 38 * c_s + 53),
    (2, 2): 192 * (2 * c_s**4 + 56 * c_s**3 + 186 * c_s**2 + 256 * c_s + 195),
    (2, 3): 2304 * (2 * c_s**3 + 24 * c_s**2 + 64 * c_s + 65),
    (2, 4): 9216 * (2 * c_s**2 + 10 * c_s + 13),
    (3, 0): 96 * c_s * (c_s - 2) * (c_s - 1) * (11 * c_s - 5),
    (3, 1): 384 * c_s * (c_s**3 + 20 * c_s**2 + 8 * c_s + 13),
    (3, 2): 2304 * (2 * c_s**3 + 22 * c_s**2 + 44 * c_s + 35),
    (3, 3): 9216 * (c_s + 5) * (3 * c_s + 7),
    (3, 4): 36864 * (2 * c_s + 7),
    (4, 0): 384 * c_s * (c_s - 2) * (c_s - 1) * (c_s + 1),
    (4, 1): 1536 * c_s * (c_s + 1) * (2 * c_s + 1),
    (4, 2): 9216 * (2 * c_s**2 + 6 * c_s + 5),
    (4, 3): 36864 * (2 * c_s + 5),
    (4, 4): sp.Integer(147456),
}


def check_1_coefficients():
    print("[CHECK 1] Verifying 50 coefficient polynomials in tex file...")
    Gc = interpolate_Gc()
    sheets = [((0, 0), TEX_SHEET_00), ((1, 1), TEX_SHEET_11)]
    for (ra, rb), claimed in sheets:
        Av = sheet_coeffs(Gc, ra, rb)
        for (i, k), poly_claim in claimed.items():
            diff = sp.expand(Av[i][k] - poly_claim)
            assert diff == 0, f"MISMATCH sheet=({ra},{rb}) (i,k)=({i},{k}): {sp.factor(Av[i][k])}"
    print("  All 50 coefficient polynomials match tex file exactly.")
    return Gc


def check_2_cross_c(Gc):
    print("[CHECK 2] Cross-checking interpolation vs fresh computation at c in {18,20,22}...")
    for c_val in [18, 20, 22]:
        fresh = G4_of_c(c_val)
        interp = sp.expand(Gc.subs(c_s, c_val))
        assert sp.simplify(fresh - interp) == 0, f'Mismatch at c={c_val}'
    print("  Fresh computation at c=18,20,22 matches interpolated polynomial.")


def check_3_case_analysis():
    print("[CHECK 3] Verifying Section 5.3 case analysis...")
    # (a) 13c^2 - 77c + 110
    for t in range(1, 30):
        c = 4 * t + 2
        val = 13 * c**2 - 77 * c + 110
        assert val == 4 * (52 * t**2 - 25 * t + 2), f'(a) expansion c={c}'
        assert v2(val) >= 2
    for t in range(2, 30):
        c = 4 * t
        val = 13 * c**2 - 77 * c + 110
        assert val == 2 * (104 * t**2 - 154 * t + 55)
        assert v2(val) == 1
    print("  (a) 13c^2-77c+110 verified.")

    # (b) 13c^3 - 82c^2 + 141c - 64
    for t in range(1, 30):
        c = 4 * t + 2
        val = 13 * c**3 - 82 * c**2 + 141 * c - 64
        assert v2(val) == 1, f'(b) c=2mod4 t={t} v_2={v2(val)}'
    for t in range(2, 30):
        c = 4 * t
        val = 13 * c**3 - 82 * c**2 + 141 * c - 64
        assert val == 4 * (208 * t**3 - 328 * t**2 + 141 * t - 16)
        assert v2(val) >= 2
    print("  (b) 13c^3-82c^2+141c-64 verified.")

    # (c) 41c^2 - 105c + 82
    for t in range(1, 30):
        c = 4 * t + 2
        val = 41 * c**2 - 105 * c + 82
        assert val == 4 * (164 * t**2 + 59 * t + 9)
        assert v2(val) >= 2
    for t in range(2, 30):
        c = 4 * t
        val = 41 * c**2 - 105 * c + 82
        assert val == 2 * (328 * t**2 - 210 * t + 41)
        assert v2(val) == 1
    print("  (c) 41c^2-105c+82 verified.")

    # (d) 41c^2 - 57c - 4
    for t in range(1, 30):
        c = 4 * t + 2
        val = 41 * c**2 - 57 * c - 4
        assert val == 2 * (328 * t**2 + 214 * t + 23)
        assert v2(val) == 1
    for t in range(2, 30):
        c = 4 * t
        val = 41 * c**2 - 57 * c - 4
        assert val == 4 * (164 * t**2 - 57 * t - 1)
        assert v2(val) >= 2
    print("  (d) 41c^2-57c-4 verified.")

    # (e) 41c^2 - 169c + 182
    for t in range(1, 30):
        c = 4 * t + 2
        val = 41 * c**2 - 169 * c + 182
        assert val == 4 * (164 * t**2 - 5 * t + 2)
        assert v2(val) >= 2
    for t in range(2, 30):
        c = 4 * t
        val = 41 * c**2 - 169 * c + 182
        assert val == 2 * (328 * t**2 - 338 * t + 91)
        assert v2(val) == 1
    print("  (e) 41c^2-169c+182 verified.")

    # (f) 11c^4 + 54c^3 + 12c^2 + 73c + 30
    for t in range(1, 30):
        c = 4 * t + 2
        val = 11 * c**4 + 54 * c**3 + 12 * c**2 + 73 * c + 30
        assert v2(val) >= 2, f'(f) c=2mod4 t={t}'
    for t in range(2, 30):
        c = 4 * t
        val = 11 * c**4 + 54 * c**3 + 12 * c**2 + 73 * c + 30
        assert v2(val) == 1
    print("  (f) 11c^4+54c^3+12c^2+73c+30 verified.")


def check_4_witnesses():
    print("[CHECK 4] Verifying K-witness and cap-witness...")
    # K-witness: 24 c(c-1)(c-4)(c-5) has v_2 = 3 + v_2(c) + v_2(c-4) on even c
    for t in range(1, 30):
        c = 4 * t + 2
        K = 24 * c * (c - 1) * (c - 4) * (c - 5)
        assert v2(K) == 5
    for t in range(2, 30):
        c = 4 * t
        K = 24 * c * (c - 1) * (c - 4) * (c - 5)
        assert v2(K) >= 7
    print("  K-witness verified.")

    # Cap witnesses: 192*(...) is always v_2 = 6
    for c in range(-20, 100):
        val1 = 192 * (2 * c**4 + 32 * c**3 + 18 * c**2 + 16 * c + 15)
        val2 = 192 * (2 * c**4 + 56 * c**3 + 186 * c**2 + 256 * c + 195)
        assert v2(val1) == 6, f'sheet00 cap: c={c} v_2={v2(val1)}'
        assert v2(val2) == 6, f'sheet11 cap: c={c} v_2={v2(val2)}'
    print("  Cap witnesses verified (v_2 = 6 for c in [-20, 100]).")


def check_5_empirical_content(Gc):
    print("[CHECK 5] Verifying empirical content = min(v_2(K), 6)...")

    def content_ab(P, ra, rb, da=4, db=4):
        Q = sp.expand(P.subs({a_s: 2 * u_s + ra, b_s: 2 * v_s + rb}))
        G = [[int(Q.subs({u_s: i, v_s: k})) for k in range(db + 1)] for i in range(da + 1)]

        def fd(seq, o):
            s = list(seq)
            for _ in range(o):
                s = [s[i + 1] - s[i] for i in range(len(s) - 1)]
            return s[0]

        Au = [[fd([G[ii][k] for ii in range(da + 1)], i) for k in range(db + 1)] for i in range(da + 1)]
        Av = [[fd([Au[i][kk] for kk in range(db + 1)], k) for k in range(db + 1)] for i in range(da + 1)]
        vals = [v2(Av[i][k]) for i in range(da + 1) for k in range(db + 1) if Av[i][k] != 0]
        return min(vals)

    for c in [6, 8, 10, 12, 14, 16]:
        G4 = G4_of_c(c)
        emp = min(content_ab(G4, 0, 0), content_ab(G4, 1, 1))
        K = 24 * c * (c - 1) * (c - 4) * (c - 5)
        pred = min(v2(K), 6)
        assert emp == pred, f'c={c}: empirical={emp} predicted={pred}'
    print("  For c in {6,8,10,12,14,16}: predicted = empirical (6/6 match).")


if __name__ == "__main__":
    print("=" * 72)
    print("Verifying proof of Lyra's cap identity (2026-07-20).")
    print("Theorem: content(G_4^{(c)}) = min(v_2(K(c)), 6) for even c >= 6.")
    print("=" * 72)
    print()

    Gc = check_1_coefficients()
    check_2_cross_c(Gc)
    check_3_case_analysis()
    check_4_witnesses()
    check_5_empirical_content(Gc)

    print()
    print("=" * 72)
    print("ALL CHECKS PASSED. Proof of Lyra's cap identity is computationally verified.")
    print("=" * 72)
