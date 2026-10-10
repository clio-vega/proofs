"""The t=1 evaluation test.  Korff's Psi_T(1) = 0 unless EVERY step has J = empty,
and J = empty iff theta is constant (all-0 or all-1) -- because a non-constant cyclic
0/1 word must contain a 0->1 ascent.  So t^{stat(T)} cannot help: it is 1 at t=1.
Prediction: c_alpha(1) = +1 if every alpha_a in {0,h}, else 0."""
import sys, itertools
sys.path.insert(0,'.')
from korff import *
from walks import AB_from_walks
import sympy as sp
t = sp.Symbol('t')

def trivial_signed_count(k, ell, alpha):
    """sum of HKKO signs over T with Psi_T(1) != 0, computed by BRUTE FORCE over all
    cylindric tableaux (no use of the theta-constant characterisation)."""
    h = 2*k; w = 2*ell+2
    tot = 0
    for P in paths(h, w, alpha):
        if sp.expand(weight_path(P, w)).subs(t, 1) == 0: continue
        la = P[-1]
        if all(la[2*i] == la[2*i+1] for i in range(k)): c = 1
        elif la[0]-la[h-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): c = -1
        else: c = 0
        tot += c
    return tot

def predicted(k, alpha):
    h = 2*k
    return 1 if all(a in (0, h) for a in alpha) else 0

for (k, ell, n) in [(1,0,4),(1,1,6),(1,0,3),(1,1,5),(2,0,4),(2,1,6),(1,2,8)]:
    h = 2*k; w = 2*ell+2
    viol = []; checked = 0; pred_ok = 0
    for alpha in itertools.product(range(h+1), repeat=n):
        A,B,TOT = AB_from_walks(k, w, alpha)
        c1 = A - B                      # = c_alpha(1), HKKO Thm 3.3 / my Lemma t1
        tsc = trivial_signed_count(k, ell, alpha)
        checked += 1
        if tsc == predicted(k, alpha): pred_ok += 1
        if c1 != tsc: viol.append((alpha, c1, tsc))
    print(f"k={k} ell={ell} n={n} (h={h},w={w}): {checked} contents; "
          f"brute-force trivial-weight signed count matches the theta-constant prediction "
          f"in {pred_ok}/{checked}; "
          f"{len(viol)} contents VIOLATE  c_alpha(1) = sum_{{Psi_T(1)!=0}} eps(T)")
    for v in viol[:4]:
        print(f"      alpha={v[0]}  c_alpha(1)={v[1]}  but trivial-weight signed count={v[2]}")
