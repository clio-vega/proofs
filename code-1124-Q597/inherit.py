"""Band-boundary inheritance.  For split records and every shadow step j in [c0+2, e]: the M0 tiles at v_j that are NOT in
H^o_j and NOT J_o-created ('lower-sector tiles', provenance X or W>0) -- does each contain v_{j-1} (was already in the star
at the previous shadow vertex)?  Also the displacement v_j - v_{j-1} and the types.   usage: python3 inherit.py SIZE"""
import sys, os, pickle, math, cmath
from collections import Counter
SIZE = int(sys.argv[1])
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
from fractions import Fraction as Fr
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import phi
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
aT = ('a1', 'a0')
def nz(t): return (t[0], frozenset(tuple(Fr(x) for x in P) for P in t[1]))
st = Counter()
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; Tp = OTHER[T]; TO = typed(Jo)
    M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); m = len(mA) - 1; e = c0 + m
    def vtx(j):
        o = min(nz(next(u for u in Jo[j][1] if u[0] == 'rh'))[1])
        return o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
    bd = set(M0)
    for j, s in enumerate(Jo):
        if c0 + 2 <= j <= e:
            v = vtx(j); vp = vtx(j - 1)
            low = [u for u in bd if u in M0 and u not in s[1] and v in nz(u)[1]]
            inh = sum(1 for u in low if vp in nz(u)[1])
            st[('T', T, 'e' if j == e else 'sh', TO[j - 1], '->', TO[j], 'dv', tuple(int(a - b) for a, b in zip(v, vp)), '#low', len(low), '#contain v_prev', inh)] += 1
        bd = (bd - s[1]) | s[2]
for k in sorted(st, key=str): print(st[k], k)
