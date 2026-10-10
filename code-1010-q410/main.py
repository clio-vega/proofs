"""Q410 main test: does Korff's cylindric HL weight carry Warnaar's generalised
affine Jacobi-Trudi coefficient?"""
import sys, itertools
sys.path.insert(0,'.')
from korff import *
from core import Engine
from warn import G_times_prodx
import sympy as sp
t = sp.Symbol('t')

def cminus(la, k, w):
    """HKKO eq (eq:wt_C) with h=k: c^-_{2k,w}(la) on la in Par(2k,w), as a tuple of length 2k."""
    m = 2*k
    if all(la[2*i] == la[2*i+1] for i in range(k)): return 1
    if la[0]-la[m-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): return -1
    return 0

def korff_side(k, ell, nvars, kind="Psi"):
    """sum_{la in Par(2k,w)} c^-(la) * sum_{T in CRST(la;2k,w)} Psi_T(t) x^T
    i.e. HKKO Thm 3.3 with the cylindric Schur function replaced by Korff's
    cylindric Hall-Littlewood function, same signs, same tableaux."""
    h = 2*k; w = 2*ell+2
    xs = sp.symbols(f'x1:{nvars+1}', positive=True)
    tot = {}
    for alpha in itertools.product(range(h+1), repeat=nvars):
        for P in paths(h, w, alpha):
            c = cminus(P[-1], k, w)
            if c == 0: continue
            mon = tuple(alpha)
            tot[mon] = tot.get(mon, 0) + c*weight_path(P, w, kind)
    return {a: sp.expand(v) for a, v in tot.items() if sp.expand(v) != 0}

def warnaar_side(k, ell, nvars):
    E = Engine(nvars)
    F = G_times_prodx(E, k, ell)
    P = sp.Poly(F, *E.x)
    out = {}
    for mon, co in P.terms():
        out[tuple(mon)] = sp.expand(co)
    return out

if __name__ == "__main__":
    for (k, ell, nvars) in [(1,0,4), (1,1,6)]:
        w = 2*ell+2; M = k+ell+1
        print(f"\n===== k={k} ell={ell} n_C={nvars}  (h={2*k}, w={w}, M={M}, N={2*M}) =====")
        print(f"      Korff dictionary: n_K = h = {2*k} cyclic slots,  k_K = w = {w} level")
        W = warnaar_side(k, ell, nvars)
        K = korff_side(k, ell, nvars)
        alpha = tuple([1]*nvars)
        ca = W.get(alpha, sp.Integer(0))
        ka = K.get(alpha, sp.Integer(0))
        print(f"  WITNESS alpha=(1^{nvars}):")
        print(f"    Warnaar c_alpha(t)          = {sp.factor(ca)}   = {sp.expand(ca)}")
        print(f"    Korff  sum_T c^- Psi_T      = {sp.factor(ka)}")
        print(f"    c_alpha(0) = {ca.subs(t,0)}   c_alpha(1) = {ca.subs(t,1)}")
        print(f"    Korff(0)   = {sp.expand(ka).subs(t,0)}   Korff(1)   = {sp.expand(ka).subs(t,1)}")
        # full comparison over all monomials
        keys = set(W) | set(K)
        agree = sum(1 for a in keys if sp.expand(W.get(a,0)-K.get(a,0)) == 0)
        print(f"  FULL: {len(keys)} monomials, {agree} agree, {len(keys)-agree} disagree")
        # t -> 1/t variant, allowing a per-monomial power of t
        ok = 0
        for a in keys:
            f = sp.expand(W.get(a,0)); g = sp.expand(K.get(a,0))
            if f == 0 and g == 0: ok += 1; continue
            gr = sp.expand(sp.simplify(g.subs(t, 1/t)*t**sp.Poly(g,t).total_degree())) if g != 0 else sp.Integer(0)
            if f == 0 or g == 0: continue
            q = sp.cancel(f/gr)
            if q.is_Pow or q == 1 or (q.free_symbols <= {t} and sp.simplify(q - t**sp.log(q,t)) == 0):
                ok += 1
        print(f"  under t->1/t with a free monomial prefactor: {ok}/{len(keys)} monomials reconcilable")
