"""P4: leading-coefficient lemma at w_n -> infinity, EXACT via sympy rational functions.

Psi is computed from eq:Psidef (head-permutation sum) with w_1..w_{n-1} fixed rationals
and w_n = symbol W.  cancel() to P/Q in lowest terms; the expansion at W=infinity has
leading term (lcP/lcQ) W^{degP-degQ}.  Compare to the prediction.
NO nsimplify anywhere (it manufactured 2**(41/58) twice this month).
"""
import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
from itertools import permutations, combinations
from psi import psi_headsum, psi_alt, mvec
from sampler import clean_w
import random

W = sp.Symbol('W')

def psi_headsum_sym(wsym, theta):
    r = len(theta); n = len(wsym)
    total = sp.Integer(0)
    for h in permutations(range(n), r):
        X = [wsym[i] for i in h]
        term = sp.Integer(1)
        for j in range(r):
            term *= (1 - 1/X[j]**2) * X[j]**theta[j]
        for j in range(r):
            for k in range(j+1, r):
                term *= (1 - 1/(X[j]*X[k])) / (1 - X[j]/X[k])
        for t in [t for t in range(n) if t not in h]:
            for j in range(r):
                term *= (1 - 1/(wsym[t]*X[j])) / (1 - wsym[t]/X[j])
        total += term
    return sp.cancel(sp.together(total))

def lead_at_infinity(expr):
    """return (d, c) with expr = c*W^d*(1+O(1/W)) at W=infinity."""
    num, den = sp.fraction(sp.cancel(expr))
    pn = sp.Poly(sp.expand(num), W); pd = sp.Poly(sp.expand(den), W)
    return pn.degree() - pd.degree(), sp.Rational(pn.LC(), pd.LC())

rng = random.Random(55921)
ok = bad = skip = 0; fails = []
caseA = caseB = 0
for trial in range(120):
    r = rng.randint(1,3)
    n = rng.randint(r, r+2)
    wq = clean_w(n, rng, lo=2, hi=12, dlo=1, dhi=5)
    theta = [rng.randint(-5,5) for _ in range(r)]
    m = mvec(theta)
    if len(set(m)) != r:
        skip += 1; continue
    jstar = max(range(r), key=lambda j: m[j])          # 0-indexed argmax, unique
    hat_m = [m[l] for l in range(r) if l != jstar]
    hat_theta = [hat_m[l] - l for l in range(r-1)]     # hat_theta_l = hat_m_l - l (0-indexed)
    wsym = [sp.Rational(x.numerator, x.denominator) for x in wq[:-1]] + [W]
    psi = psi_headsum_sym(wsym, theta)
    if psi == 0:
        fails.append(("PSI IDENTICALLY ZERO IN W", r, n, theta, m)); bad += 1; continue
    d, c = lead_at_infinity(psi)
    # prediction
    pred_d = m[jstar] - r + 1
    Mp = wq[:-1]
    pred_c = sp.Integer(0)
    if r >= 1:
        v = psi_alt(Mp, hat_theta) if r-1 > 0 else F(1)
        pred_c = sp.Rational(v.numerator, v.denominator) * (-1)**(r + (jstar+1))
    if m[jstar] >= 0: caseA += 1
    else: caseB += 1
    if d == pred_d and sp.simplify(c - pred_c) == 0:
        ok += 1
    else:
        bad += 1
        if len(fails) < 6:
            fails.append((r, n, theta, m, jstar+1, "got", d, c, "pred", pred_d, pred_c))
print("P4 leading-coeff lemma: OK=%d BAD=%d  (skipped %d with repeated m)" % (ok, bad, skip))
print("   of the tested cases: max(m)>=0 in %d, max(m)<=-1 in %d" % (caseA, caseB))
for f in fails: print("   FAIL", f)
