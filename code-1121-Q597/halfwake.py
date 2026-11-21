"""PROVE 11-21: half-wake shift. On the 11-20 journey dumps (not regenerated), for every step: (target, class index rel. start,
type) -> multiset of triangle shift vectors (0-triangles: shift of q in Lambda_1; 1-triangles: shift of p in Lambda_0),
read via bs6.sdelta (cells matched between removed and added sets).  Frame: u_A=(0,-1;0,0), w_A=(0,0;0,1)."""
import pickle, glob, sys
from collections import Counter
from fractions import Fraction as Fr
src = open('/home/clio/projects/proofs/code-1118c2-Q597/bs6.py').read()
# only the pure helper functions (no census, no engine import side effects beyond what they need)
ns = {}
exec(src[src.index('def pre(A, B):'):src.index('import migrate')], ns)
rh, sdelta = ns['rh'], ns['sdelta']
def l0(s):
    vs = list(rh(s[1]))
    for P in vs:
        for Q in vs:
            d = tuple(Fr(b) - Fr(a) for a, b in zip(P, Q))
            if any(d[:2]) and not any(d[2:]):
                return d[:2] if (d[0], d[1]) > (0, 0) else tuple(-x for x in d[:2])
def cen(vs): return tuple(sum(Fr(v[i]) for v in vs)/len(vs) for i in range(4))
st = Counter(); nj = Counter()
for f in sorted(glob.glob('/home/clio/projects/proofs/code-1120-Q597/dumpj5-*.pkl')):
    d = pickle.load(open(f, 'rb'))
    for kind in ('BC', 'BB2'):
        for TT, M0, J in d[kind]:
            nj[kind] += 1; c0 = l0(J[0])
            for s in J:
                cl = 0 if l0(s) == c0 else 1
                if s[0] == 'b': t = 'b'
                else:
                    a = tuple(y - x for x, y in zip(cen(rh(s[1])), cen(rh(s[2]))))
                    t = 'a1' if not any(a[:2]) else 'a0'
                D = sdelta(s)
                if D is None: st[(kind, TT, cl, t, 'sdelta None')] += 1; continue
                for (side, cell), v in D.items():
                    st[(kind, TT, cl, t, '%d-tri shift' % side, tuple(int(x) for x in v))] += 1
print('journeys', dict(nj))
for k in sorted(st, key=str): print(st[k], k)
