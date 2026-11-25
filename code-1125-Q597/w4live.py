"""(W4) liveness + window: for every middle word (r <= RMAX), every translation vector v in {0} u 12 unit vectors u {-alpha}:
does S = H_p+gamma meet (A's step p+i) + v for some 1<=i<=m?  The claim is the v = -alpha row (must be 0).
Also the window: for v = -alpha, at which i is the 4-offset test (om(S+alpha tile, R_{p+i}) in Q_T) passed at all (must be i in {m-1, m})."""
import sys, itertools
from collections import Counter
from wmodel import *
RMAX = int(sys.argv[1]); st = Counter()
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
VS = [(0, 0, 0, 0)] + [u + (0, 0) for u in UN] + [(0, 0) + u for u in UN]
for T in 'BC':
    R_p = canon(T); aT = AT[T]
    remA, _, Rp1 = step(T, R_p, aT); alpha = tuple(int(x) for x in vsub(W.base(Rp1), W.base(R_p)))
    QT = set()
    for (ck, typ), (Rc, rem, Rn, cre) in CF[T].items():
        for t in rem: QT.add(W.offs(t, Rc, T))
    for r in range(RMAX + 1):
        for xs in itertools.product(('a0', 'a1'), repeat=r):
            word = [(1, aT), (1, 'b')] + [(0, x) for x in xs] + [(0, 'b')]
            steps, Rend = run_word(T, R_p, word); m = len(word) - 1
            gamma = vsub(W.base(Rend), W.base(steps[1][0]))
            S = frozenset(W.trs(t, gamma) for t in remA)
            for v in VS:
                hit = any(S & frozenset(W.trs(t, v) for t in steps[i][1]) for i in range(1, m + 1))
                st[(T, 'v', v, 'is -alpha' if v == tuple(-x for x in alpha) else '', 'some C+v step meets S', hit)] += 1
            Sa = [W.trs(t, alpha) for t in S]
            feas = tuple(m - i for i in range(1, m + 1) if any(W.offs(t, steps[i][0], T) in QT for t in Sa))
            st[(T, 'window: m-i at which some S+alpha tile has offset in Q_T', feas, 'last two letters', word[-2][1] if m >= 2 else None)] += 1
for k in sorted(st, key=str): print(st[k], k)
