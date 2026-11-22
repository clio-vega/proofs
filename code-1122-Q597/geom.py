"""PROVE 11-22: shared local geometry (exec'd by qn.py, ccwin.py). Integer tiles ('kind', sorted point tuple), the
configuration table of Theorem AW (code-1120c2-Q597/words.configs, unchanged), plane overlap by SAT (as twoword.py),
4-offsets om (as words.offs)."""
import sys, os
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597')
sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
_cwd = os.getcwd(); os.chdir('/home/clio/projects/proofs/code-1120c2-Q597')
from fractions import Fraction as Fr
from collections import Counter
import words as Wd
import localenum as LE
from psilib import psi
os.chdir(_cwd)
TYPES = ('a0', 'a1', 'b')
OTHER = {'B': 'C', 'C': 'B'}
def it(t): return (t[0], tuple(sorted(tuple(int(x) for x in P) for P in t[1])))
def tr(t, v): return (t[0], tuple(sorted(tuple(a+b for a, b in zip(P, v)) for P in t[1])))
CF = {}
for _T in 'BC':
    for (ck, typ), (Rc, rem, Rn, cre) in Wd.configs(_T).items():
        for P in [q for t in list(rem) + list(cre) + [Rc, Rn] for q in t[1]]: assert all(Fr(x).denominator == 1 for x in P)
        CF[(_T, ck, typ)] = (it(Rc), [it(t) for t in rem], it(Rn), [it(t) for t in cre])
K0 = frozenset({(0, 1), (0, -1)})
K1 = {'B': frozenset({(1, -1), (-1, 1)}), 'C': frozenset({(1, 0), (-1, 0)})}
def classkey(R):
    d = set()
    for A in R[1]:
        for B in R[1]:
            w = tuple(b-a for a, b in zip(A, B))
            if any(w[:2]) and not any(w[2:]): d.add(w[:2])
    return frozenset(d)
def cls(T, R): return 0 if classkey(R) == K0 else (1 if classkey(R) == K1[T] else None)
def base(R): return min(R[1])
def canon(T, c): return CF[(T, K0 if c == 0 else K1[T], 'a0')][0]
def step(T, R, typ):
    """removed tiles (incl. R), new rhombus, created tiles (incl. new rhombus) of the good T-configuration of type typ at R"""
    Rc, rem, Rn, cre = CF[(T, classkey(R), typ)]
    sh = base(R); assert tr(Rc, sh) == R
    return [tr(t, sh) for t in rem] + [R], tr(Rn, sh), [tr(t, sh) for t in cre] + [tr(Rn, sh)]
def nonrh(S): return [t for t in S if t[0] != 'rh']
POLY = {}
def poly(t):
    if t not in POLY: POLY[t] = LE.convex_hull([LE.phi(P) for P in t[1]])
    return POLY[t]
OV = {}
def overlap(s, t):
    """True iff the interiors of the plane polygons intersect (separating axis test)."""
    key = (s, t) if s < t else (t, s)
    if key in OV: return OV[key]
    A, B = poly(s), poly(t); res = True
    for P in (A, B):
        for i in range(len(P)):
            e = P[(i+1) % len(P)] - P[i]; n = complex(-e.imag, e.real)
            pa = [(z*n.conjugate()).real for z in A]; pb = [(z*n.conjugate()).real for z in B]
            if max(pa) <= min(pb) + 1e-9 or max(pb) <= min(pa) + 1e-9: res = False; break
        if not res: break
    OV[key] = res; return res
def meets(t, S): return any(t == s or overlap(t, s) for s in S)
def compatible(t, S): return all(t == s or not overlap(t, s) for s in S)
def offs(t, R, T):
    """om_T(t, R) = (psi_T(t) - psi_T(R), lambda(t) - lambda(R)), lambda = (min a0, min a1) -- as words.offs"""
    f = lambda x: (x[0], frozenset(x[1]))
    pt, pr = psi(f(t), T), psi(f(R), T)
    at = (min(P[0] for P in t[1]), min(P[2] for P in t[1])); ar = (min(P[0] for P in R[1]), min(P[2] for P in R[1]))
    return (int(pt[0]-pr[0]), int(pt[1]-pr[1]), int(at[0]-ar[0]), int(at[1]-ar[1]))
