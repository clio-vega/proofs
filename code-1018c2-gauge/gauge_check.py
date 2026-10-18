"""PROVE 2026-10-18 c2 -- verification for the GAUGE THEOREM

    Atil(R u D) = ( prod_{y in D} phi(y) ) * Atil(R),      phi(y) = (1-q)/(1-q y^h)

for every INVERSE-CLOSED multiset D of spectral values (h=2 type C, h=1 type B).

The three root families of C_N are the three inversion-orbit types of D:
    e_i - e_j   D = {}            (empty)          -> phi-product 1   (= u-symmetry)
    2 e_m       D = {eps}, eps^2=1 (iota-FIXED pt)  -> phi(eps) = 1
    e_a + e_b   D = {x, 1/x}       (iota-FREE orbit) -> phi(x)phi(1/x) = rho^C

Every ratio is checked SEPARATELY against its closed form (F-ratio, E-ratio, Atil-ratio):
a quotient can be 1 with both halves wrong.  Every control carries an ADMISSIBLE plant
drawn from OUTSIDE the hypothesis class (a D that is not inverse-closed), not from inside.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from itertools import permutations, combinations
from math import factorial
from q484_engine import q, wvars, F, EC, EB, EX, WX

# ---------------------------------------------------------------- generic Atilde(U)
def WU(N, U, w, I, X):
    """W^X with the u-arguments given explicitly as the list U (NOT tied to w_M)."""
    n, r = len(U), len(I)
    if n < r:
        return sp.Integer(0)   # Motegi L2843: Atil = 0 when |M| < |I|.  The engine's
                               # n==0 branch returns 1, which is only right for r==0;
                               # caught by the first smoke test (d=2, r=1, R empty).
    if n == 0:
        return sp.Integer(1)
    return WX(N, n, r, list(U), list(w), tuple(I), X)

def EU(N, U, w, I, X):
    n, r = len(U), len(I)
    return EX(N, n, r, list(U), list(w), X)

def AtilU(N, U, w, I, X):
    """normalised weight function at an arbitrary value multiset U."""
    return sp.cancel(sp.together(WU(N, U, w, I, X) / EU(N, U, w, I, X)))

def FU(N, U, w, I, X):
    """the symmetrisation F only (no 1/E), for the separate F-ratio check."""
    n, r = len(U), len(I)
    if n < r:
        return sp.Integer(0)
    if n == 0:
        return sp.Integer(1)
    beta, gamma = (sp.Integer(1), sp.Integer(0)) if X == 'C' else (sp.Integer(1), sp.Integer(-1))
    return F(N, n, r, list(U), list(w), tuple(I), beta, gamma)

# ---------------------------------------------------------------- closed forms
def phi(y, h):
    return (1 - q) / (1 - q * y**h)

def pred_Atil_ratio(D, h):
    out = sp.Integer(1)
    for y in D:
        out *= phi(y, h)
    return sp.cancel(sp.together(out))

def pred_F_ratio(N, R, D, w):
    """prod_{y in D}[ prod_{v in R}(1-qy/v) prod_{k<=N}(1-q y w_k) ] * P2(D) * (prod y)^{1-|D|}

    (type C; derived in the paper, Step 3).  P2(D)=prod over unordered distinct
    slot pairs of D of (y y' - q).
    """
    d = len(D)
    out = sp.Integer(1)
    for y in D:
        for v in R:
            out *= (1 - q * y / v)
        for k in range(N):
            out *= (1 - q * y * w[k])
    for a in range(d):
        for b in range(a + 1, d):
            out *= (D[a] * D[b] - q)
    py = sp.Integer(1)
    for y in D:
        py *= y
    out *= py**(1 - d)
    return sp.cancel(sp.together(out))

def Q1(D):
    out = sp.Integer(1)
    for a in D:
        for b in D:
            out *= (1 - q * b / a)
    return out

def Q2(D):
    out = sp.Integer(1)
    for a in range(len(D)):
        for b in range(a + 1, len(D)):
            out *= (1 - q * D[a] * D[b])
    return out

def P2(D):
    out = sp.Integer(1)
    for a in range(len(D)):
        for b in range(a + 1, len(D)):
            out *= (D[a] * D[b] - q)
    return out

def pred_E_ratio(N, R, D, w, X):
    """Q1(D)/Q2(D)/(1-q)^d * prod_{y in D}[prod_{v in R}(1-qy/v) prod_k(1-q y w_k)]

    (type C, derived in the paper Step 4, AFTER the iota-cancellation
     prod_{y,v}(1-qv/y) <-> prod_{y,v}(1-q y v)).
    """
    assert X == 'C'
    d = len(D)
    out = Q1(D) / Q2(D) / (1 - q)**d
    for y in D:
        for v in R:
            out *= (1 - q * y / v)
        for k in range(N):
            out *= (1 - q * y * w[k])
    return sp.cancel(sp.together(out))

# ---------------------------------------------------------------- the checks
def eq(a, b):
    return sp.simplify(sp.cancel(sp.together(a - b))) == 0

def adjunction(N, R, D, w, I, X, h, want_F=True, want_E=True):
    """Return dict of verdicts for one (N,R,D,I,X) instance."""
    U = list(R) + list(D)
    out = {}
    aU, aR = AtilU(N, U, w, I, X), AtilU(N, R, w, I, X)
    if sp.simplify(aR) == 0:
        # degenerate stratum: both sides vanish, 0 = lambda*0 carries no information.
        # Report it, never score it -- a vacuous green is the fault of
        # [[a-plant-drawn-from-inside-its-hypothesis-class]].
        return {'VACUOUS': bool(sp.simplify(aU) == 0)}
    out['Atil'] = eq(sp.cancel(aU / aR), pred_Atil_ratio(D, h))
    if want_F:
        fU, fR = FU(N, U, w, I, X), FU(N, R, w, I, X)
        out['F'] = eq(sp.cancel(fU / fR), pred_F_ratio(N, R, D, w)) if fR != 0 else None
    if want_E and X == 'C':
        eU, eR = EU(N, U, w, I, X), EU(N, R, w, I, X)
        out['E'] = eq(sp.cancel(eU / eR), pred_E_ratio(N, R, D, w, X))
    return out
