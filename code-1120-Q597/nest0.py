"""PROVE 11-20: k=0.  For the first step of every journey, what lies across the two back edges E_u, E_w of R_0
(canonical R_0=[0,v]x[0,v], E_u={v}x[0,v], E_w=[0,v]x{0}), split by whether the nest rhombus across exists."""
import sys, pickle, glob
from collections import Counter
sys.argv = [sys.argv[0], '5', '1', '1']
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
from stepset import index
import os; os.chdir('/home/clio/projects/proofs/code-1120-Q597')
st = Counter()
for fn in sorted(glob.glob('dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in D[fam]:
            Bd = set(M0); E = index(Bd)
            R = next(t for t in J[0][1] if t[0] == 'rh')
            pts0 = sorted({P[:2] for P in R[1]}); pts1 = sorted({P[2:] for P in R[1]})
            v = tuple(b-a for a, b in zip(*pts0)); b0, b1 = pts0[0], pts1[0]
            add = lambda P, Q: tuple(x+y for x, y in zip(P, Q))
            Eu = frozenset({add(b0, v) + b1, add(b0, v) + add(b1, v)})
            Ew = frozenset({b0 + b1, add(b0, v) + b1})
            for nm, e in (('E_u', Eu), ('E_w', Ew)):
                across = E.get(e, set()) - {R}
                st[(fam, TT, 'v', v, nm, tuple(sorted(t[0] for t in across)) or 'NOTHING (region boundary)')] += 1
for k in sorted(st, key=str): print(st[k], k)
