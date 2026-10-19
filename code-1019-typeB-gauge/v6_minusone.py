"""V6.  The eps = -1 tau-fixed point in type B, computed with the ANALYTICALLY
CANCELLED form of W^B.

Why this file exists.  W^B := prod_k (1+u_k^{-1})^{-1} * F|_{gamma=-1}.  The
prefactor has a POLE at u_k = -1 and line 3c has a ZERO there, so substituting
u_k = -1 into the literal expression gives zoo*0 = nan.  V4 did exactly that and
reported 0/18 at eps=-1 -- which is an undefined computation, NOT a refutation.

The cancelled form.  Per permutation sigma, with v_j = u_sigma(j),
  prod_{k=1}^n u_k/(u_k+1) * T_sigma^{(gamma=-1)}
    = (prod_k u_k) * [lines 1,2a,2b,3a,3b] * prod_{j>r}(v_j+1) / prod_{j=1}^n (v_j+1)
    = (prod_k u_k) * [lines] / prod_{j<=r}(v_j+1),
and the surviving 1/(v_j+1) cancels line 1's (v_j^2-1) = (v_j-1)(v_j+1):

   TYPE B HEAD FACTOR IS (v_j - 1)   where   TYPE C's IS (v_j^2 - 1).

That is 1 - e^{-alpha} for the SHORT boundary root e_j versus the LONG one 2e_j,
and it has a consequence for HEAD EXCLUSION:
   (v_j^2-1) vanishes at v_j = +1 AND -1   -> type C excludes both fixed points
   (v_j - 1) vanishes at v_j = +1 ONLY     -> type B excludes only +1
So eps = -1 is predicted to be a tau-fixed point at which type-B head exclusion
FAILS.  This file measures what actually happens to the adjunction law there.
"""
import sympy as sp
from itertools import permutations, combinations
from math import factorial
from q484_engine import EB, WX, q, wvars

phiB = lambda v: (1 - q) / (1 - q * v)


def TB(N, n, r, us, w, I):
    """T_sigma for type B in the cancelled form: head factor (v-1), not (v^2-1)."""
    t = sp.Integer(1)
    for j in range(r):
        ij = I[j]
        t *= (1 - q) * (us[j] - 1)                       # (v^2-1)/(v+1)
        for k in range(ij + 1, N + 1):
            t *= (us[j] - q * w[k - 1])
        for k in range(1, ij):
            t *= (us[j] - w[k - 1])
    for j in range(r):
        for k in range(j + 1, r):
            t *= (q * us[j] - us[k]) / (us[j] - us[k]) * (us[j] * us[k] - 1)
    for j in range(r, n):
        for k in range(N):
            t *= (us[j] - q * w[k])
    for j in range(r, n):
        for k in range(r):
            t *= (q * us[j] - us[k]) / (us[j] - us[k]) * (us[j] * us[k] - 1)
    for j in range(r, n):
        for k in range(j + 1, n):
            t *= (us[j] * us[k] - q)
    return t


def WB_safe(N, n, r, u, w, I):
    tot = sp.Integer(0)
    for sg in permutations(range(n)):
        tot += TB(N, n, r, [u[sg[j]] for j in range(n)], w, I)
    pref = sp.Integer(1)
    for k in range(n):
        pref *= u[k] ** (-(N + n - 1))
    return sp.together(pref * tot / factorial(n - r))


def AtB(N, r, u, w, I):
    n = len(u)
    if n == 0:
        return sp.cancel(1 / EB(N, 0, 0, [], w))
    return sp.cancel(sp.together(WB_safe(N, n, r, u, w, I) / EB(N, n, r, u, w)))


def z(e):
    return sp.simplify(sp.cancel(sp.together(e))) == 0


# ---- step 0: WB_safe agrees with the engine's WX(...,'B') at generic u -------
print("=== V6.0  WB_safe == WX(...,'B') at generic u (consistency of the rewrite)")
okc = tot = 0
for N in (2, 3):
    w = wvars(N)
    for n in (1, 2):
        for r in range(0, n + 1):
            for I in combinations(range(1, N + 1), r):
                u = [sp.Rational(7, 2), sp.Rational(9, 4)][:n]
                tot += 1
                okc += 1 if z(WB_safe(N, n, r, u, w, I) - WX(N, n, r, u, w, I, 'B')) else 0
print(f"    {okc}/{tot}")

# ---- step 1: head exclusion at eps = +1 and eps = -1, measured ---------------
print("\n=== V6.1  head exclusion at a tau-fixed point: does the head term vanish?")
for eps in (sp.Integer(1), sp.Integer(-1)):
    for N in (2,):
        w = wvars(N)
        n, r = 2, 1
        for I in combinations(range(1, N + 1), r):
            # sigma placing eps in the HEAD slot (slot 0)
            head_term = TB(N, n, r, [eps, sp.Rational(7, 2)], w, I)
            print(f"    eps={str(eps):>2}  I={I}  head-slot term == 0 ?  {z(head_term)}")

# ---- step 2: the adjunction law at D = {-1}, cancelled form ------------------
print("\n=== V6.2  the adjunction law at D containing -1 (cancelled form)")
CASES = [("d1 eps=+1", [sp.Integer(1)]),
         ("d1 eps=-1", [sp.Integer(-1)]),
         ("d2 {+1,-1}", [sp.Integer(1), sp.Integer(-1)]),
         ("d2 free {x,1/x}", [sp.Integer(3), sp.Rational(1, 3)]),
         ("d3 {-1,x,1/x}", [sp.Integer(-1), sp.Integer(3), sp.Rational(1, 3)])]
for name, D in CASES:
    rows = []
    meas = []
    for N in (2, 3):
        w = wvars(N)
        for Rl in ([], [sp.Rational(7, 2)]):
            s_ = len(Rl)
            if s_ + len(D) > 3:
                continue
            for r in range(0, s_ + 1):
                for I in combinations(range(1, N + 1), r):
                    base = AtB(N, r, Rl, w, I)
                    if z(base):
                        continue
                    lhs = AtB(N, r, Rl + D, w, I)
                    pred = sp.prod([phiB(v) for v in D])
                    rows.append((N, r, str(I), z(lhs - pred * base)))
                    if not rows[-1][3]:
                        meas.append((N, r, I, sp.cancel(sp.together(lhs / base))))
    g = sum(1 for a in rows if a[3])
    print(f"    {name:18s} predicted phi={sp.cancel(sp.prod([phiB(v) for v in D]))}"
          f"   ->  {g}/{len(rows)}")
    for m in meas[:4]:
        print(f"        MEASURED at N={m[0]},r={m[1]},I={m[2]}:  {sp.nsimplify(m[3])}")
