r"""ROUTE (iv) PROBE.  Test Conjecture R:  for m notin A,
      Atilde^X_{A,B}(w_1..w_N)|_{w_m = eps}  =?=  c * Atilde^X_{A',B'}(w')
where w' = (w_j)_{j != m} relabelled to [N-1], A' = image of A, and
B' = image of B\{m}.   c is whatever rational function the ratio turns out to be.

We do NOT guess c.  We compute the RATIO symbolically and factor it.
DOMAIN: N = 2,3;  X in {B,C};  every m;  every admissible eps;
        every (A,B) with A >= B and m notin A.  Fully symbolic in q and the
        remaining w's -- no sampling.

GATE (brief instrument fault 3): before substituting w_m = eps, check that the
substitution point is in the domain, i.e. the reduced denominator of Atilde does
not vanish identically at w_m = eps.
"""
import sys, itertools
import sympy as sp
sys.path.insert(0, '.')
import q484_engine as E
q = E.q

def relabel(S, m):
    """[N]\{m} -> [N-1] order-preserving; S must avoid m."""
    return tuple(j if j < m else j - 1 for j in S)

def restrict(expr, wm, eps):
    """Substitute w_m = eps after checking the point is in the domain."""
    f = sp.cancel(sp.together(expr))
    num, den = sp.fraction(f)
    d0 = sp.simplify(den.subs(wm, eps))
    if d0 == 0:
        return None, 'DENOM VANISHES at w_m=eps'
    return sp.cancel(num.subs(wm, eps) / d0), 'ok'

rows = []
for N in (2, 3):
    w = E.wvars(N)
    wp = E.wvars(N - 1) if N > 1 else []
    P = E.subsets(N)
    for X in ('C', 'B'):
        h = 2 if X == 'C' else 1
        epss = [sp.Integer(1), sp.Integer(-1)] if h == 2 else [sp.Integer(1)]
        for m in range(1, N + 1):
            for eps in epss:
                for A in P:
                    if m in A: continue
                    for B in P:
                        if not E.geq(A, B): continue
                        val, st = restrict(E.Atilde(N, w, A, B, X=X), w[m-1], eps)
                        if val is None:
                            rows.append((N,X,m,int(eps),A,B,'DOMAIN-FAIL',None)); continue
                        zero = sp.simplify(val) == 0
                        # the rank-(N-1) comparison object
                        Ap, Bp = relabel(A, m), relabel(tuple(b for b in B if b != m), m)
                        # map the surviving w's: w_j (j!=m) -> wp
                        sub = {w[j-1]: wp[relabel((j,), m)[0]-1] for j in range(1, N+1) if j != m}
                        vsub = val.subs(sub, simultaneous=True)
                        low = E.Atilde(N-1, wp, Ap, Bp, X=X) if N-1 >= 0 else sp.Integer(1)
                        if sp.simplify(low) == 0:
                            ratio = 'LOW=0'
                        else:
                            ratio = sp.factor(sp.cancel(sp.together(vsub / low)))
                        rows.append((N,X,m,int(eps),A,B,'ZERO' if zero else 'nonzero',ratio))

print("N X m eps   A        B        status    ratio  Atilde|_S / Atilde^{N-1}_{A',B'}")
print("-"*100)
for r in rows:
    N,X,m,eps,A,B,st,ratio = r
    cl = 'c3' if A==B else ('c2' if m in B else 'c1')
    print(f"{N} {X} {m} {eps:+d}  {str(A):9s}{str(B):9s}{cl} {st:8s} {ratio}")
nz = sum(1 for r in rows if r[6]=='nonzero')
z  = sum(1 for r in rows if r[6]=='ZERO')
df = sum(1 for r in rows if r[6]=='DOMAIN-FAIL')
print(f"\nSUMMARY  nonzero {nz}   ZERO {z}   DOMAIN-FAIL {df}   total {len(rows)}")
