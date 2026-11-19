"""Windows (Lemma window): N(R,T) = R + non-rh tiles within 3 hops of R through non-rh tiles.
For each step of the clean middle in M (boards M_1.. after a) and in Mx: is N(R_i - alpha, Mx_i) == N(R_i, M_i) - alpha (C part),
and for the last Mx step: N == N(a in M_P) + gamma?  Where do the window tiles of C come from (wake D / a-new / common-outside)?"""
from replay import *
from collections import Counter, defaultdict
from itertools import combinations
UNIT = {(1, 0), (0, 1), (-1, 0), (0, -1), (1, -1), (-1, 1)}
def is_edge(P, Q):
    d = tuple(b-a for a, b in zip(P, Q))
    return (d[2:] == (0, 0) and d[:2] in UNIT) or (d[:2] == (0, 0) and d[2:] in UNIT)
def edges(t):
    return {frozenset(e) for e in combinations(t[1], 2) if is_edge(*e)}
def index(Bd):
    E = defaultdict(set)
    for t in Bd:
        for e in edges(t): E[e].add(t)
    return E
def window(Bd, R, E=None):
    E = E or index(Bd); N = {R}; fr = {R}
    for _ in range(3):
        nf = set()
        for t in fr:
            for e in edges(t):
                for u in E[e]:
                    if u not in N and u[0] != 'rh': nf.add(u)
        N |= nf; fr = nf
    return N
def T(S, v): return {(k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S}
st = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side, A, B, Jo = clean(r)
    mA, mB, p, q = middles(A, B)
    M0 = board0(r); MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
    rem = MP - MxP
    a, C = mA[0], mA[1:]
    alpha = tuple(y-x for x, y in zip(cen(rh(a[1])), cen(rh(a[2])))); ma = tuple(-v for v in alpha)
    gamma = tuple(y-x for x, y in zip(cen(rh(C[0][1])), cen(rh(C[-1][2]))))
    Bm = apply(MP, [a]); Bx = set(MxP); allok = True; orig = Counter()
    for s, t in zip(C, mB[:-1]):
        Rm = next(u for u in s[1] if u[0] == 'rh'); Rx = next(u for u in t[1] if u[0] == 'rh')
        Wm = window(Bm, Rm); Wx = window(Bx, Rx)
        allok &= (T(Wm, ma) == Wx)
        for u in Wm: orig['D' if u in rem else ('y-made' if u not in MP else 'common')] += 1
        Bm = apply(Bm, [s]); Bx = apply(Bx, [t])
    st[(side, 'C windows translate by -alpha', allok)] += 1
    st[(side, 'C windows contain common (non-wake, not y-made) tiles', orig['common'] > 0)] += 1
    Ra = next(u for u in a[1] if u[0] == 'rh'); Rl = next(u for u in mB[-1][1] if u[0] == 'rh')
    Wa = window(MP, Ra); Wl = window(Bx, Rl)
    st[(side, 'last window = a-window + gamma', T(Wa, gamma) == Wl)] += 1
    st[(side, 'a-window: #D, #common', sum(u in rem for u in Wa), sum(u not in rem for u in Wa))] += 1
for k in sorted(st, key=str): print(st[k], k)
