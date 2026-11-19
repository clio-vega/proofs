"""(1) validate stepset.step against recorded steps (infer d among 6 unit directions);
(2) candidate-level covariance along the clean middles."""
from replay import *
from stepset import *
import cmath
from collections import Counter
DIRS = [cmath.exp(1j*cmath.pi*k/12) for k in range(24)]
def T(S, v): return frozenset((k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S)
st = Counter(); n = 0
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    n += 1
    if n > 200: break
    side, A, B, Jo = clean(r)
    mA, mB, p, q = middles(A, B)
    M0 = board0(r); MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
    # infer d from the middle in M
    dd = None
    for d in DIRS:
        Bm = set(MP); good = True
        for s in mA:
            Rm = next(u for u in s[1] if u[0] == 'rh'); got, _ = step(Bm, Rm, d)
            if got is None or (got[0], frozenset(got[1]), got[2]) != (s[0], s[1], s[2]): good = False; break
            Bm = apply(Bm, [s])
        if good: dd = d; break
    st[(side, 'reimpl reproduces M-middle for some d', dd is not None)] += 1
    if dd is None: continue
    a, C = mA[0], mA[1:]
    alpha = tuple(y-x for x, y in zip(cen(rh(a[1])), cen(rh(a[2])))); ma = tuple(-v for v in alpha)
    gamma = tuple(y-x for x, y in zip(cen(rh(C[0][1])), cen(rh(C[-1][2]))))
    Bm = apply(MP, [a]); Bx = set(MxP); ok = True; ncand = Counter()
    for s, t in zip(C, mB[:-1]):
        Rm = next(u for u in s[1] if u[0] == 'rh'); Rx = next(u for u in t[1] if u[0] == 'rh')
        cm = {T(S, ma) for _, S, _ in candidates(Bm, Rm, dd)}; cx = {S for _, S, _ in candidates(Bx, Rx, dd)}
        ok &= (cm == cx); ncand[len(cm)] += 1
        Bm = apply(Bm, [s]); Bx = apply(Bx, [t])
    st[(side, 'C candidates translate by -alpha', ok)] += 1
    st[(side, 'candidates per C step', tuple(sorted(ncand.items())))] += 1
    Ra = next(u for u in a[1] if u[0] == 'rh'); Rl = next(u for u in mB[-1][1] if u[0] == 'rh')
    ca = {T(S, gamma) for _, S, _ in candidates(MP, Ra, dd)}; cl = {S for _, S, _ in candidates(Bx, Rl, dd)}
    st[(side, 'last candidates = a-candidates + gamma', ca == cl)] += 1
    st[(side, '#cands at a', len(ca))] += 1
for k in sorted(st, key=str): print(st[k], k)
