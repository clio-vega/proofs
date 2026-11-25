"""(T#) probe: the re-run's last-middle rhombus R_p+gamma.  Frame coordinates (i,j) in nest T's frame (packed_at logic, i,j
any integers), membership of (i,j) in nest T (pos), addability.  Same for A's rhombus R_{p+m+1}.  usage: python3 tsharp.py SIZE"""
import sys, pickle
from fractions import Fraction as Fr
from collections import Counter
SIZE = int(sys.argv[1]); sys.argv = [sys.argv[0], str(SIZE), '1', '1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1123c2-Q597/nest%d.pkl' % SIZE, 'rb'))
def shift(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
UB, UC = (-1, 1, 0, 0), (1, 0, 0, 0)
def fcoord(N, R):
    X, u, w = N['X'], N['u'], N['w']; S = set(R[1])
    for p in S:
        if {p, tuple(a+b for a, b in zip(p, u)), tuple(a+b for a, b in zip(p, w)), tuple(a+b+c for a, b, c in zip(p, u, w))} == S:
            D = tuple(a-b for a, b in zip(p, X)); Du = D[:2] + (0, 0); Dw = (0, 0) + D[2:]
            i = migrate.Board.coef(Du, u) if any(Du) else Fr(0); j = migrate.Board.coef(Dw, w) if any(Dw) else Fr(0)
            if i is None or j is None: return 'off-lattice'
            if i.denominator != 1 or j.denominator != 1: return 'non-integral'
            return (int(i), int(j))
    return 'wrong shape'
def ideal(S): return all(((i-1, j) in S or i == 0) and ((i, j-1) in S or j == 0) for (i, j) in S)
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    if P['kind'] if 'kind' in P else False: pass
    roles = P['roles']
    if len(roles) != 2: continue
    try: T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in roles]
    except Exception: continue
    Ls = [i for i in (0, 1) if T_[i] == (1, aT[i])]
    if len(Ls) != 1: continue
    L = Ls[0]; t, A, Bq, Jo = roles[L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    TT = 'B' if L == 0 else 'C'
    N = next(N for N in P['nests'] if N['u'] == (UB if TT == 'B' else UC))
    pos = set(map(tuple, N['pos']))
    gam = disp(mA[1], mA[-1]); al = disp(mA[0], mA[0]); Rp = ('rh', frozenset(rh(A[p][1])))
    Rpg = shift(Rp, gam); Rsh = shift(Rpg, al)
    for name, R in (('R_p+gamma', Rpg), ('R_{p+m+1}', Rsh), ('R_p', Rp)):
        c = fcoord(N, R)
        if isinstance(c, tuple):
            add = c not in pos and ideal(pos | {c}) and c[0] >= 0 and c[1] >= 0
            st[(TT, name, 'frame cell', 'in nest' if c in pos else ('ADDABLE (packed)' if add else 'not addable'), 'i<0 or j<0' if min(c) < 0 else '')] += 1
        else: st[(TT, name, c)] += 1
    st[(TT, 'A ends at R_{p+m+1}', len(A) == p + len(mA))] += 1
for k in sorted(st, key=str): print(st[k], k)
