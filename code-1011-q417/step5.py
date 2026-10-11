"""Q417 STEP 5: two structural claims.

PROP 6 (contradicts the brief's own description of Psi_open).  The brief says Psi_open
"keeps the wrap where it is MEASURED (inside the multiplicities m_h = w-(x_1-x_h)) -- so it
is NOT ordinary psi, the cylinder is still present."  But J_open(theta) is a subset of
[1,h-1], and for i <= h-1 Korff's m_i(x) = x_i - x_{i+1} IS the ordinary multiplicity.
So Psi_open NEVER READS m_h.  CLAIM: Psi_open_T = Macdonald's ordinary psi_T (III (5.11'))
evaluated on the cylindric path; the cylinder enters ONLY through the region R, i.e. as an
ADMISSIBILITY CONSTRAINT on which strips are legal, never as a factor of the weight.
TEST: compute an ordinary psi that never mentions w, and compare.  Firing control: the
same comparison against Korff's CYCLIC Psi, which MUST disagree.

CLAIM 7: the minimum number of J_open-nonempty strips over the alpha=(1^n) fibre is 1,
which is why ord_(t=1) of the signed sum is exactly 1 -- short of the requirement by
exactly ONE order, where the cyclic index set was short by 2M.
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
from korff import paths, J_open, weight_path, mvec, t
import sympy as sp

def thetas(P):
    h = len(P[0])
    return [tuple(P[a+1][i]-P[a][i] for i in range(h)) for a in range(len(P)-1)]

def psi_ordinary(P):
    """Macdonald III (5.11') psi for the horizontal strip x -> x+theta, with NO reference
    to w whatsoever: psi = prod over i in [1,h-1] with theta_i=0, theta_{i+1}=1 of
    (1 - t^{x_i - x_{i+1}}).  The letter w does not appear in this function."""
    h = len(P[0]); out = sp.Integer(1)
    for a in range(len(P)-1):
        x = P[a]; th = tuple(P[a+1][i]-x[i] for i in range(h))
        for i in range(h-1):
            if th[i] == 0 and th[i+1] == 1:
                out *= (1 - t**(x[i]-x[i+1]))
    return sp.expand(out)

print("="*74)
print("PROP 6 -- is Psi_open equal to ordinary psi (wrap-blind)?")
print("="*74)
agree = dis = 0; cyc_agree = cyc_dis = 0
for (k, ell) in [(1,0),(1,1),(1,2),(2,0),(2,1)]:
    h = 2*k; w = 2*ell+2; M = k+ell+1
    for n in [2*M]:
        for alpha in [(1,)*n] + list(itertools.islice(
                itertools.product(range(h+1), repeat=min(n,3)), 40)):
            for P in paths(h, w, alpha):
                a = weight_path(P, w, "Psi_open"); b = psi_ordinary(P)
                if sp.expand(a-b) == 0: agree += 1
                else:
                    dis += 1
                    if dis <= 3: print(f"   DISAGREE k={k} ell={ell} P={P} Psi_open={a} psi={b}")
                c = weight_path(P, w, "Psi")
                if sp.expand(c-b) == 0: cyc_agree += 1
                else: cyc_dis += 1
print(f"   Psi_open vs ordinary psi : {agree} agree, {dis} DISAGREE")
print(f"   FIRING CONTROL, Korff's CYCLIC Psi vs ordinary psi: {cyc_agree} agree, "
      f"{cyc_dis} disagree  (must be >0, else the test is blind)")
print(f"   => Psi_open is wrap-blind: the cylinder enters ONLY via the region R."
      if dis == 0 else "   => claim FALSE")

print()
print("="*74)
print("CLAIM 7 -- minimum number of J_open-nonempty strips over the alpha=(1^n) fibre")
print("="*74)
def cminus(la, k, w):
    m = 2*k
    if all(la[2*i] == la[2*i+1] for i in range(k)): return 1
    if la[0]-la[m-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): return -1
    return 0
def ord_t1(f):
    f = sp.expand(f)
    if f == 0: return None
    u = sp.Symbol('u'); P = sp.Poly(sp.expand(f.subs(t, 1+u)), u)
    return min(m[0] for m, c in zip(P.monoms(), P.coeffs()) if c != 0)
for (k, ell) in [(1,0),(1,1),(1,2),(1,3),(2,1)]:
    h = 2*k; w = 2*ell+2; M = k+ell+1; n = 2*M
    mn = None; wit = None; S = sp.Integer(0)
    for P in paths(h, w, (1,)*n):
        b = sum(1 for th in thetas(P) if J_open(th))
        if mn is None or b < mn: mn, wit = b, P
        c = cminus(P[-1], k, w)
        if c: S += c*weight_path(P, w, "Psi_open")
    print(f"   k={k} ell={ell} M={M} n={n}: min #bad strips = {mn}   "
          f"ord_(t=1) of signed sum = {ord_t1(S)}   (target ord = "
          f"{ord_t1(sp.catalan(M)-(2*M-1)*t**(M-1)+t**(M+1))})")
    print(f"      witness attaining the minimum: {wit}")
    print(f"         its thetas: {thetas(wit)}  Psi_open = {weight_path(wit,w,'Psi_open')}")
