"""PROVE 11-23 c2, (Ea) on the nest census: in case a, does the engine's termination test fire at R_{p+1} (A) or R'_{p+1}
(re-run), and where are R_p, R_{p+1}, R'_{p+1}, R_{p+2} in nest T's cell grid?  The test is re-implemented from
migrate.py (packed_at + is_ideal over nests X != src).
CONTROL: the test must fire at the final rhombus of every journey (A, re-run, J_o), and at no earlier rhombus.
usage: python3 eaA.py SIZE"""
import sys, os, pickle
from fractions import Fraction as Fr
from collections import Counter
SIZE = int(sys.argv[1]); sys.argv = [sys.argv[0], str(SIZE), '1', '1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1123c2-Q597/nest%d%s.pkl' % (SIZE, os.environ.get('TAG', '')), 'rb'))
def coef(D, v):
    k = None
    for a, b in zip(D, v):
        if b == 0:
            if a != 0: return None
        else:
            kk = Fr(a) / b
            if k is None: k = kk
            elif kk != k: return None
    return k
def padd(a, b): return tuple(x + y for x, y in zip(a, b))
def cellpos(info, R):
    """(i, j) as Fractions if R is the u,w-parallelogram at X + i u + j w (i, j any sign), else None"""
    X, u, w = info['X'], info['u'], info['w']; S = set(R[1]) if isinstance(R, tuple) else set(R)
    for q in S:
        if {q, padd(q, u), padd(q, w), padd(padd(q, u), w)} == S:
            D = tuple(a - b for a, b in zip(q, X))
            Du = (D[0], D[1], 0, 0) if u[2:] == (0, 0) else (0, 0, D[2], D[3]); Dw = tuple(a - b for a, b in zip(D, Du))
            i = coef(Du, u) if any(Du) else Fr(0); j = coef(Dw, w) if any(Dw) else Fr(0)
            if i is None or j is None: continue
            return (i, j)
    return None
def is_ideal(S): return all(((i-1, j) in S or i == 0) and ((i, j-1) in S or j == 0) for (i, j) in S)
def stops(nests, R, src=0):
    for X, info in enumerate(nests):
        if X == src: continue
        c = cellpos(info, R)
        if c is None or c[0] < 0 or c[1] < 0 or c[0].denominator != 1 or c[1].denominator != 1: continue
        ij = (int(c[0]), int(c[1]))
        if ij not in info['pos'] and is_ideal(set(info['pos']) | {ij}): return X
    return None
def rhs(J): return [rh(s[1]) for s in J] + [rh(J[-1][2])]
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    N0 = P['nests']
    N1 = [dict(N0[i], pos=P['nests2' if True else ''][i]['pos']) for i in range(3)]   # nests after J_y (board b2)
    N2 = [dict(N0[i], pos=P['nests1'][i]['pos']) for i in range(3)]                   # nests after J_x (board b1)
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    for r, (t, A, Bq, Jo) in enumerate(P['roles']):
        Nre = N1 if r == 0 else N2     # role x: re-run on board after J_y; role y: after J_x
        for name, J, NN in (('A', A, N0), ('rerun', Bq, Nre), ('Jo', Jo, N0)):
            R_ = rhs(J)
            st[('CONTROL', name, 'fires at end', stops(NN, R_[-1]) is not None,
                'fires earlier', any(stops(NN, R) is not None for R in R_[:-1]))] += 1
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] == 'b': continue
    TT = 1 if L == 0 else 2; Nre = N1 if L == 0 else N2
    Rp, Rp1, Rq1 = rh(A[p][1]), rh(A[p][2]), rh(Bq[p][2]); Rp2 = rh(A[p+1][2])
    cp = [cellpos(N0[TT], R) for R in (Rp, Rp1, Rq1, Rp2)]
    st[('T', 'BC'[L], 'stop at R_p+1', stops(N0, Rp1), 'stop at R\'_p+1 (rerun board)', stops(Nre, Rq1))] += 1
    st[('T', 'BC'[L], 'cell coords (R_p, R_p+1, R\'_p+1, R_p+2) rel. to', tuple(None if c is None else ((c[0] - cp[0][0], c[1] - cp[0][1]) if cp[0] else c) for c in cp),
        'R_p coords integral', cp[0] is not None and all(x.denominator == 1 for x in cp[0]))] += 1
    if cp[0] is not None:
        i, j = cp[0]; pos = set(N0[TT]['pos'])
        st[('T', 'BC'[L], 'R_p cell (i,j) in pos', (i, j) in pos, 'i<0', i < 0, 'j<0', j < 0,
            'R_p+1 cell in pos', (cp[1][0], cp[1][1]) in pos if cp[1] else None,
            'R\'_p+1 cell in pos', (cp[2][0], cp[2][1]) in pos if cp[2] else None,
            'R_p+2 cell in pos', (cp[3][0], cp[3][1]) in pos if cp[3] else None)] += 1
for k in sorted(st, key=str): print(st[k], k)
