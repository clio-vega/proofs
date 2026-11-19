"""Full run of candidate covariance (all 521) + half-wake structure:
rem = x's wake tiles before (M_P \\ Mx_P), add = after.  H = {t in rem : t - alpha in add} (the translating half).
Are all tiles of all C-candidates (in M_1..) and a-candidates that lie in rem, inside H?  |H|/|rem|?"""
from replay import *
from stepset import *
import cmath
from collections import Counter
DIRS = [cmath.exp(1j*cmath.pi*k/12) for k in range(24)]
def T(S, v): return frozenset((k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S)
st = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side, A, B, Jo = clean(r)
    mA, mB, p, q = middles(A, B)
    M0 = board0(r); MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
    rem, add = MP - MxP, MxP - MP
    a, C = mA[0], mA[1:]
    alpha = tuple(y-x for x, y in zip(cen(rh(a[1])), cen(rh(a[2])))); ma = tuple(-v for v in alpha)
    gamma = tuple(y-x for x, y in zip(cen(rh(C[0][1])), cen(rh(C[-1][2]))))
    H = {t for t in rem if T([t], ma) <= add}
    st[(side, '|H|/|rem| bucket', round(len(H)/len(rem), 1))] += 1
    dd = None
    for d in DIRS:
        Bm = set(MP); good = True
        for s in mA:
            Rm = next(u for u in s[1] if u[0] == 'rh'); got, _ = step(Bm, Rm, d)
            if got is None or (got[0], got[1], got[2]) != (s[0], s[1], s[2]): good = False; break
            Bm = apply(Bm, [s])
        if good: dd = d; break
    st[(side, 'direction found', dd is not None)] += 1
    if dd is None: continue
    Bm = apply(MP, [a]); Bx = set(MxP); ok = True; inH = True; outside = Counter()
    for s, t in zip(C, mB[:-1]):
        Rm = next(u for u in s[1] if u[0] == 'rh'); Rx = next(u for u in t[1] if u[0] == 'rh')
        cs = candidates(Bm, Rm, dd)
        ok &= ({T(S, ma) for _, S, _ in cs} == {S for _, S, _ in candidates(Bx, Rx, dd)})
        for _, S, _ in cs:
            for u in S:
                if u[0] == 'rh': continue
                outside['H' if u in H else ('rem-notH' if u in rem else ('y-made' if u not in MP else 'common'))] += 1
        Bm = apply(Bm, [s]); Bx = apply(Bx, [t])
    st[(side, 'C cands translate', ok)] += 1
    st[(side, 'C-candidate non-rh tiles: kinds', tuple(sorted(outside)))] += 1
    Ra = next(u for u in a[1] if u[0] == 'rh'); Rl = next(u for u in mB[-1][1] if u[0] == 'rh')
    st[(side, 'last cands = a cands + gamma', {T(S, gamma) for _, S, _ in candidates(MP, Ra, dd)} == {S for _, S, _ in candidates(Bx, Rl, dd)})] += 1
for k in sorted(st, key=str): print(st[k], k)
