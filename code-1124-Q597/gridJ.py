"""Is J_o already on its own stop-grid line at R^o_{e+1}?  F_eps^{T'}(R) = f_eps^{T'}(min R) (Lemma grid, 11-23 c2):
f0B=a0+b0, f1B=a1+b1, f0C=b0, f1C=b1.  s' = coordinate moved by J_o's straight step: T'=C straight (1,a0) moves F0;
T'=B straight (1,a1) moves F1.  Compare F(R^o_{e+1}) with F at J_o's final rhombus (= nest value by Lemma grid).
Also at R* = R^o_e and the same for A at R# vs A's final.  CONTROL: generic (0,b) steps of J_o followed by straight.
usage: python3 gridJ.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE = int(sys.argv[1]); sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
FG = {'B': (lambda P: P[0] + P[1], lambda P: P[2] + P[3]), 'C': (lambda P: P[1], lambda P: P[3])}
SP = {'C': 0, 'B': 1}; STRT = {'B': (1, 'a1'), 'C': (1, 'a0')}; aT = ('a1', 'a0'); OTH = {'B': 'C', 'C': 'B'}
def rhof(S): return next(u for u in S if u[0] == "rh")
def Fv(T, R): m = min(R[1]); return (FG[T][0](m), FG[T][1](m))
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    for role in (0, 1):
        t, A, Bq, Jo = P['roles'][role]; T = 'BC'[role]; Tp = OTH[T]; TO = typed(Jo)
        cfin = Fv(Tp, rhof(Jo[-1][2]))
        mA, mB, p, q = middles(A, Bq)
        split = None
        if role == L and typed(Bq)[p][1] == 'b':
            M0 = set(P['M0']); W = M0 - set(apply(M0, Jo)); rem = {}
            for j, s in enumerate(Jo):
                for u in s[1]: rem.setdefault(u, j)
            split = min(rem[u] for u in A[p][1] if u in W) + len(mA) - 1
        for j, s in enumerate(Jo):
            if TO[j] != (0, 'b') or j + 1 >= len(Jo): continue
            v = Fv(Tp, rhof(s[2])); d = (v[0] - cfin[0], v[1] - cfin[1])
            lab = 'SPLIT' if j == split else ('generic next=straight' if TO[j+1] == STRT[Tp] else 'generic next=other')
            st[(lab, 'Tp', Tp, "F_s' - final", d[SP[Tp]], 'other F - final', d[1 - SP[Tp]])] += 1
for k in sorted(st, key=str): print(st[k], k)
