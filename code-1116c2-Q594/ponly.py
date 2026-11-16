"""Cayley-trick rigidity test (pre-registered).  Pool every state reachable (full BFS, single A->B / A->C moves, as rigid2.py)
from every census start, grouped by nu (same canonical hexagon AND same nest-C content; cross-start comparison is meaningful
because coordinates come from canonical(nt,...)).  Keys tested for 'determines the tiling':
  K_qp = (q,p)  [11-16: injective on 860 states]   K_p = p only   K_q = q only   K_nest = nest data (control: RIGID fails)
Report #states, #classes, #classes with >1 tiling.  Non-vacuity: #K_p classes containing states from >1 start / >1 nest data."""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1116-Q585-MC')
from qp import *
from collections import defaultdict, Counter
pool = defaultdict(dict)   # nt -> {tiling: (q,p,nestdata,startkey)}
for key, b0 in census():
    lam, mu, nu, bb = key; nt = tuple(nu)
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    seen = set(); stack = [b0]
    while stack:
        b = stack.pop(); T = tuple(tiles_of(b))
        if T in seen: continue
        seen.add(T)
        q, p = qp(b)
        pool[nt][T] = (frozenset(q.items()), frozenset(p.items()), nestdata(b), key)
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
st = Counter()
for nt, D in pool.items():
    for name, f in (('qp', lambda v: (v[0], v[1])), ('p', lambda v: v[1]), ('q', lambda v: v[0]), ('nest', lambda v: v[2])):
        g = defaultdict(set); src = defaultdict(set)
        for T, v in D.items(): g[f(v)].add(T); src[f(v)].add(v[3])
        st[(name, 'classes')] += len(g); st[(name, 'classes >1 tiling')] += sum(len(s) > 1 for s in g.values())
        st[(name, 'classes from >1 start')] += sum(len(s) > 1 for s in src.values())
    st['states'] += len(D); st['pools'] += 1
for k in sorted(st, key=str): print(k, st[k])
