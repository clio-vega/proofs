"""Trunk lemma check on the uncapped size-5 closure, EVERY nest-A corner x of every state:
K = x->B, H = x->C step logs (IBC.STEPLOG records: typ, R_before, R_after, removed, created).
 T1 every step before the first type-b step of K is type a, and K, H agree on all of them (same removed AND created);
 T2 the first type-b step of K and of H has the same index f and the same removed set (hexagon);
 T3 at f the created sets differ, and the new rhombus classes are B-class for K and C-class for H;
 T4 no landing before f (journey length > f).
usage: python3 trunk.py S0 S1"""
import sys
exec(open('/home/clio/projects/proofs/code-1128c2-Q585-SW/swfull.py').read().split("if __name__")[0])
from fractions import Fraction as Fr
def jl(b, x, dv):
    c = copy.deepcopy(b); IBC.STEPLOG.clear(); t = c.nest[0]['pos'].pop(x); c.journey(t, 0, unit_dir(dv)); return list(IBC.STEPLOG)
def ucls(vs):
    ps = sorted(set((P[0]+P[2], P[1]+P[3]) for P in vs)); d = (ps[-1][0]-ps[0][0], ps[-1][1]-ps[0][1])
    d = (Fr(d[0], 2), Fr(d[1], 2)); return d if d > (0, 0) else (-d[0], -d[1])
st = Counter(); S0, S1 = int(sys.argv[1]), int(sys.argv[2])
for k, (key, b0) in enumerate(census()):
    if k < S0 or k >= S1: continue
    for b, dirs in closure(b0):
        for x in corners(set(b.nest[0]['pos'])):
            K = jl(b, x, dirs['B']); H = jl(b, x, dirs['C']); st['corners'] += 1
            fK = next(i for i, s in enumerate(K) if s[0] == 'b'); fH = next(i for i, s in enumerate(H) if s[0] == 'b')
            st[('T1', all(K[i][0] == 'a' and K[i][3] == H[i][3] and K[i][4] == H[i][4] for i in range(min(fK, fH))) and fK == fH)] += 1
            st[('T2 same hexagon at f', K[fK][3] == H[fH][3])] += 1
            st[('T3 created differ', K[fK][4] != H[fH][4], 'classes', ucls(K[fK][2]), ucls(H[fH][2]))] += 1
            st[('f', fK)] += 1
            st[('class before f', tuple(sorted(set(ucls(K[i][1]) for i in range(fK+1)))))] += 1
for k in sorted(st, key=str): print(st[k], k)
