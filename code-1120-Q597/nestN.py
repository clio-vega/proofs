"""PROVE 11-20: fact (N) on census START boards.  For every occupied nest-A cell (i,j) (nest[0]):
the back edges E_u = {X0+i u} x s1_j and E_w = s0_i x {X1+j w} of its rhombus have, on the far side,
the nest rhombus (i-1,j) / (i,j-1) if i>0 / j>0, and NO tile (region boundary) if i=0 / j=0.
usage: python3 nestN.py SIZE"""
import sys
SZ = sys.argv[1]; sys.argv = [sys.argv[0], SZ, '1', '1']
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
from stepset import index
from collections import Counter
st = Counter(); nb = 0
for key, b in census_sized():
    nb += 1; Bd = bare(b); E = index(Bd); A = b.nest[0]
    X, u, w = A['X'], A['u'], A['w']
    add = lambda P, Q: tuple(x+y for x, y in zip(P, Q)); mul = lambda k, v: tuple(k*x for x in v)
    for (i, j), tid in A['pos'].items():
        R = (b.tiles[tid][0], frozenset(b.tiles[tid][1]))
        base = add(X, add(mul(i, u), mul(j, w)))
        Eu = frozenset({base, add(base, w)}); Ew = frozenset({base, add(base, u)})
        for nm, e, back in (('E_u', Eu, i), ('E_w', Ew, j)):
            assert e <= R[1], 'edge not on R'
            far = E.get(e, set()) - {R}
            kinds = tuple(sorted(t[0] for t in far))
            st[(nm, 'index>0' if back > 0 else 'index=0', kinds)] += 1
print('SIZE', SIZE, 'start boards', nb)
for k in sorted(st, key=str): print(st[k], k)
