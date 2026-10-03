"""(c) the m=2 CONCAVITY mechanism, and exactly where it stops.
   (f) how many REAL slices are one-sided (zero horizontal-strip defect throughout)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen, diffsys
from onesided import F
from gen import is_pf2
from itertools import product
from collections import Counter

def concave_interior(f):
    """f(a-1)+f(a+1) <= 2 f(a) for every a with f(a-1),f(a+1) > 0."""
    for a in range(1, len(f)-1):
        if f[a-1] > 0 and f[a+1] > 0 and f[a-1] + f[a+1] > 2*f[a]: return False
    return True

print("--- (c) is F CONCAVE on its support interior? ---", flush=True)
for m in (2, 3, 4):
    st = Counter(); wit = []
    Um = {2: 8, 3: 5, 4: 4}[m]
    for L in product(range(1, 3), repeat=m):
        for U in product(*[range(L[i], Um+1) for i in range(m)]):
            for T in range(sum(L), sum(U)+1):
                f, nW = F(list(L), list(U), T)
                if f is None or len(f) < 3: continue
                st['cases'] += 1
                if concave_interior(f): st['concave'] += 1
                else:
                    st['NOT_concave'] += 1
                    if len(wit) < 2: wit.append((L, U, T, f))
    print(f"  m={m}: {dict(st)}", flush=True)
    for w in wit: print("     not concave:", w, flush=True)

print("--- (f) real slices with ZERO horizontal-strip defect throughout ---", flush=True)
st = Counter()
for m in range(2, 6):
    for n in range(m, m+5):
        for (mu, lam) in gen.pairs(n, m, 7):
            a_, g, P, Q = diffsys.data(mu, lam, n, m)
            for b in range(0, 2*n+1):
                sig = sum(mu) + b - sum(a_)
                Y = [y for y in product(*[range(P[i], Q[i]+1) for i in range(m)])
                     if sum(y) == sig]
                Yeff = [y for y in Y
                        if all(g[i]+1-max(y[i],0)-max(-y[(i+1) % m],0) >= 1 for i in range(m))]
                if not Yeff: continue
                st['slices'] += 1
                st['nontrivial'] += (1 if len(Yeff) >= 2 else 0)
                if all(min(y) >= 0 for y in Yeff):
                    st['one_sided'] += 1
                    if len(Yeff) >= 2: st['one_sided_nontrivial'] += 1
                    st['one_sided_maxsize'] = max(st.get('one_sided_maxsize', 0), len(Yeff))
                hw = {sum(g[i]+1-max(y[i],0)-max(-y[(i+1) % m],0) for i in range(m)) for y in Yeff}
                if len(hw) == 1:
                    st['single_halfwidth'] += 1
                    if all(min(y) >= 0 for y in Yeff): st['single_hw_and_onesided'] += 1
print(dict(st), flush=True)
