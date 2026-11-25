"""Lesion of the mosaic hypothesis at its point of use in (T#) (PROVE 2026-11-25 c2).
(T#) uses one consequence of "M0 = M0(ell) is a Purbhoo mosaic": nest T's frame origin X lies on the side line of Delta_k
(Cor img + Thm FB frame).  Lesion: keep the board and nest T's packed cells, translate the frame origin X by v, for v in the
12 unit vectors +-e_i of Z^4 and v=0 (control).  Then evaluate the engine's stop test packed_at(T, R_p+gamma) +
addability (journey's condition).  (T#) predicts: never fires at v=0.  If it fires for some v, the re-run would stop at
R_p+gamma on that (non-mosaic) frame, CT2 and BS(L) fail there.
usage: python3 framelesion.py SIZE   (case-b records from code-1123c2-Q597/nest{SIZE}.pkl)"""
import sys, pickle
from fractions import Fraction as Fr
from collections import Counter
SIZE = int(sys.argv[1]); sys.argv = [sys.argv[0], str(SIZE), '1', '1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1123c2-Q597/nest%d.pkl' % SIZE, 'rb'))
def shift(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
UB, UC = (-1, 1, 0, 0), (1, 0, 0, 0)
def fires(N, X, R, pos):
    u, w = N['u'], N['w']; S = set(R[1])
    for p in S:
        if {p, tuple(a+b for a, b in zip(p, u)), tuple(a+b for a, b in zip(p, w)), tuple(a+b+c for a, b, c in zip(p, u, w))} == S:
            D = tuple(a-b for a, b in zip(p, X)); Du = D[:2] + (0, 0); Dw = (0, 0) + D[2:]
            i = migrate.Board.coef(Du, u) if any(Du) else Fr(0); j = migrate.Board.coef(Dw, w) if any(Dw) else Fr(0)
            if i is None or j is None: continue
            if i >= 0 and j >= 0 and i.denominator == 1 and j.denominator == 1:
                c = (int(i), int(j))
                if c not in pos and migrate.Board.is_ideal(pos | {c}): return True
    return False
VS = [(0, 0, 0, 0)] + [tuple(s if k == i else 0 for k in range(4)) for i in range(4) for s in (1, -1)]
aT = ('a1', 'a0'); st = Counter(); nrec = Counter()
for P in PAIRS:
    roles = P['roles']
    if len(roles) != 2: continue
    try: T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in roles]
    except Exception: continue
    Ls = [i for i in (0, 1) if T_[i] == (1, aT[i])]
    if len(Ls) != 1: continue
    L = Ls[0]; t, A, Bq, Jo = roles[L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    TT = 'B' if L == 0 else 'C'; nrec[TT] += 1
    N = next(N for N in P['nests'] if N['u'] == (UB if TT == 'B' else UC))
    pos = set(map(tuple, N['pos']))
    gam = disp(mA[1], mA[-1]); al = disp(mA[0], mA[0]); Rp = ('rh', frozenset(rh(A[p][1])))
    Rpg = shift(Rp, gam)
    for v in VS:
        X = tuple(a+b for a, b in zip(N['X'], v))
        if fires(N, X, Rpg, pos): st[(TT, 'frame shift', v, '-alpha' if tuple(-a for a in al) == v else '')] += 1
print('case-b records', dict(nrec))
print('records where the stop test FIRES at R_p+gamma, per frame shift v (absent = 0):')
for k in sorted(st, key=str): print(st[k], k)
