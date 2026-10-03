"""DUAL FORM.  Summing over y first:

   k(a,b) = sum_{r in Box_g, sum r = a + sigma}  N(r),
   N(r)   = # { y : max(P_i, r_{i-1}-g_{i-1}) <= y_i <= min(Q_i, r_i),  sum y_i = sigma }

So k(.,b) is the SLICE-SUM TRANSFORM of ONE FIXED nonnegative function N on the
box Box_g = prod_i [0,g_i].  Condition (A) = "that transform is log-concave".

This file (i) calibrates the dual form against the primal, and (ii) tests which
discrete-concavity property N has.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen, diffsys
from itertools import product
from collections import Counter

def N_fun(P, Q, g, sig, r):
    m = len(g)
    lo = [max(P[i], r[(i-1) % m] - g[(i-1) % m]) for i in range(m)]
    hi = [min(Q[i], r[i]) for i in range(m)]
    if any(lo[i] > hi[i] for i in range(m)): return 0
    # coefficient of z^sig in prod (z^lo_i + ... + z^hi_i)
    cur = {0: 1}
    for i in range(m):
        nxt = {}
        for s, c in cur.items():
            for t in range(lo[i], hi[i]+1):
                nxt[s+t] = nxt.get(s+t, 0) + c
        cur = nxt
    return cur.get(sig, 0)

def k_dual(mu, lam, n, m, b):
    a_, g, P, Q = diffsys.data(mu, lam, n, m)
    sig = sum(mu) + b - sum(a_)
    out = Counter()
    for r in product(*[range(0, g[i]+1) for i in range(m)]):
        v = N_fun(P, Q, g, sig, r)
        if v: out[sum(r) - sig] += v
    return dict(out)

# ---------- M^natural-concavity test ----------
def Mnat_concave(f, box, report=False):
    """f: dict r -> value (>=0, 0 = outside dom).  box: list of (lo,hi).
       M^nat-concave: for all x,y in dom, all i with x_i>y_i, EITHER
         f(x)+f(y) <= f(x-e_i)+f(y+e_i)   OR
         exists j with x_j<y_j: f(x)+f(y) <= f(x-e_i+e_j)+f(y+e_i-e_j)
       (values outside dom treated as -inf, i.e. the inequality fails there)"""
    import math
    dom = [r for r, v in f.items() if v > 0]
    m = len(box)
    LG = {r: math.log(v) for r, v in f.items() if v > 0}
    def val(r): return LG.get(tuple(r), None)
    bad = []
    for x in dom:
        for y in dom:
            for i in range(m):
                if x[i] <= y[i]: continue
                lhs = LG[x] + LG[y]
                xm = list(x); xm[i] -= 1
                yp = list(y); yp[i] += 1
                a, bb = val(xm), val(yp)
                if a is not None and bb is not None and lhs <= a + bb + 1e-12:
                    continue
                ok = False
                for j in range(m):
                    if x[j] >= y[j]: continue
                    xx = list(x); xx[i] -= 1; xx[j] += 1
                    yy = list(y); yy[i] += 1; yy[j] -= 1
                    a2, b2 = val(xx), val(yy)
                    if a2 is not None and b2 is not None and lhs <= a2 + b2 + 1e-12:
                        ok = True; break
                if not ok:
                    bad.append((x, y, i))
                    if len(bad) > 3: return bad
    return bad

if __name__ == '__main__':
    st = Counter()
    for m in range(2, 5):
        for n in range(m, m+4):
            for (mu, lam) in gen.pairs(n, m, 6):
                for b in range(0, 2*n+1):
                    p = {k: v for k, v in diffsys.k_diffsys(mu, lam, n, m, b).items() if v}
                    d = {k: v for k, v in k_dual(mu, lam, n, m, b).items() if v}
                    st['agree' if p == d else 'DISAGREE'] += 1
    print("dual form vs primal:", dict(st))

    print("\n--- is N  M^natural-concave?  (real shapes) ---")
    st2 = Counter(); wit = []
    for m in range(2, 5):
        for n in range(m, m+4):
            for (mu, lam) in gen.pairs(n, m, 6):
                a_, g, P, Q = diffsys.data(mu, lam, n, m)
                for b in range(0, 2*n+1):
                    sig = sum(mu) + b - sum(a_)
                    f = {}
                    for r in product(*[range(0, g[i]+1) for i in range(m)]):
                        v = N_fun(P, Q, g, sig, r)
                        if v: f[r] = v
                    if len(f) < 2: st2['trivial'] += 1; continue
                    bad = Mnat_concave(f, [(0, g[i]) for i in range(m)])
                    st2['Mnat_ok' if not bad else 'Mnat_FAIL'] += 1
                    st2['domsize'] = max(st2.get('domsize', 0), len(f))
                    if bad and len(wit) < 3:
                        wit.append((m, n, mu, lam, b, g, P, Q, bad[:2], dict(sorted(f.items()))))
    print(dict(st2))
    for w in wit:
        print(f"   m={w[0]} n={w[1]} mu={w[2]} lam={w[3]} b={w[4]} g={w[5]} P={w[6]} Q={w[7]}")
        print(f"      bad triples {w[8]}")
        print(f"      N = {w[9]}")
