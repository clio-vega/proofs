#!/usr/bin/env python3
"""
Job B: fit the closed form for c=2 (validate) and c=3.

Ansatz (generalising c=1,c=2):
   M_j = (a-b+1) * C(N, b-j) * P(j) / D(j)
where D(j) is a product of linear factors (a+s-j),(b+t-j),... and a constant,
and P(j) is a polynomial of degree 2c in j (c=1: quadratic, c=2: quartic, c=3: sextic).

Strategy: for each candidate (N, D), extract  P(j) = M_j * D(j) / [(a-b+1) C(N,b-j)]
as exact Fractions over the FULL range j=0..m, and test whether P is a polynomial of
degree exactly 2c in j.  Report the (N,D) that works for ALL tested (a,b).
"""
from fractions import Fraction as Fr
from math import comb, factorial
from itertools import product as iproduct
from job_jstar_engine import Mj_chain_all
from sympy import symbols, interpolate, Poly, factor, simplify, expand

js = symbols('j')

def poly_from_points(pts):
    """Lagrange-interpolate Fraction points {j: value}; return sympy Poly in j and its degree."""
    expr = interpolate([(k, v) for k, v in sorted(pts.items())], js)
    P = Poly(expr, js)
    return P, P.degree()

def extract_P(a, b, c, N, dfactors, dconst):
    """Return dict j-> P(j) over j where the binomial arg b-j in [0,N]; None if denom 0."""
    lam = (a, b, c); m = sum(lam)//2
    M = Mj_chain_all(lam)
    pts = {}
    for j in range(m+1):
        if b-j < 0 or b-j > N:
            continue
        binom = comb(N, b-j)
        if binom == 0:
            continue
        D = dconst
        for (coef_a, coef_b, coef_c, const) in dfactors:
            D *= (coef_a*a + coef_b*b + coef_c*c + const - j)
        if D == 0:
            return None
        num = M[j] * D
        den = (a-b+1) * binom
        if den == 0:
            return None
        pts[j] = Fr(num, den)
    return pts

def test_ansatz(c, N_choices, denom_choices, shapes, target_deg, verbose=False):
    """For each (N, denom) try to fit; require P polynomial of degree=target_deg on all shapes,
    AND identical degree.  Return list of working (Nname, dname)."""
    working = []
    for (Nname, Nfun) in N_choices:
        for (dname, dfactors, dconst) in denom_choices:
            ok_all = True
            degs = set()
            for (a, b) in shapes:
                lam = (a, b, c)
                if sum(lam) % 2 != 0 or not (a >= b >= c >= 1):
                    continue
                N = Nfun(a, b, c, sum(lam)//2)
                pts = extract_P(a, b, c, N, dfactors, dconst)
                if pts is None or len(pts) < target_deg+2:
                    ok_all = False; break
                P, deg = poly_from_points(pts)
                degs.add(deg)
                if deg != target_deg:
                    ok_all = False; break
            if ok_all and degs == {target_deg}:
                working.append((Nname, dname))
                if verbose:
                    print(f"   WORKS: N={Nname}, D={dname}")
    return working

# candidate N functions
N_choices = [
    ('a+2', lambda a,b,c,m: a+2),
    ('a+1', lambda a,b,c,m: a+1),
    ('a+3', lambda a,b,c,m: a+3),
    ('m',   lambda a,b,c,m: m),
    ('2m-3j-no', lambda a,b,c,m: a+2),  # placeholder dup
    ('a',   lambda a,b,c,m: a),
]

def build_denoms(c):
    """denominator = const * product of (row + k - j).  Try combos around the hook structure."""
    # linear factors specified as (coef_a,coef_b,coef_c,const): value = coef*..+const - j
    pool = {
        'a+3': (1,0,0,3), 'a+2': (1,0,0,2), 'a+1': (1,0,0,1),
        'b+2': (0,1,0,2), 'b+1': (0,1,0,1), 'b':   (0,1,0,0),
        'c+1': (0,0,1,1), 'c':   (0,0,1,0),
    }
    choices = []
    # c=2 memory: 2*(a+3-j)(b+2-j)(b+1-j)
    if c == 2:
        for const in [1, 2]:
            choices.append((f'{const}*(a+3)(b+2)(b+1)', [pool['a+3'],pool['b+2'],pool['b+1']], const))
            choices.append((f'{const}*(b+2)(b+1)', [pool['b+2'],pool['b+1']], const))
    if c == 3:
        # try analogues: one a-factor + two/three b/c factors, const 1,2,6,12
        fa = [pool['a+3'], pool['a+2'], pool['a+1']]
        fb_sets = [
            ('(b+2)(b+1)', [pool['b+2'],pool['b+1']]),
            ('(b+2)(b+1)b', [pool['b+2'],pool['b+1'],pool['b']]),
            ('(b+1)b', [pool['b+1'],pool['b']]),
            ('(b+2)(b+1)(c+1)', [pool['b+2'],pool['b+1'],pool['c+1']]),
        ]
        for af in fa:
            an = {1:'a+1',2:'a+2',3:'a+3'}[af[3]]
            for (bn, bs) in fb_sets:
                for const in [1,2,6,12]:
                    choices.append((f'{const}*({an})*{bn}', [af]+bs, const))
    return choices

def run():
    print("="*78); print("VALIDATE METHOD ON c=2 (target degree 4)"); print("="*78)
    shapes2 = [(a,b) for a in range(4,16) for b in range(2,a+1) if (a+b)%2==0 and b>=2][:40]
    w2 = test_ansatz(2, N_choices, build_denoms(2), shapes2, 4, verbose=True)
    print("c=2 working ansaetze (N, D):", w2)

    print("\n"+"="*78); print("FIT c=3 (target degree 6)"); print("="*78)
    shapes3 = [(a,b) for a in range(4,18) for b in range(3,a+1) if (a+b)%2==1 and b>=3]
    w3 = test_ansatz(3, N_choices, build_denoms(3), shapes3, 6, verbose=True)
    print("c=3 working ansaetze (N, D):", w3)
    return w2, w3

if __name__ == "__main__":
    run()
