"""Local swap laws on all reachable mixed states (cap 60/cell), h = A-corner, g = cell directly BELOW h in its column (row coord = j).
(S1) landing:  C-landing of h in (h->C)  ==  C-landing of g in (h->B ; g->C)
(S2) tiling :  (h->B ; g->C)  ==  (h->C ; g->B)   [same tiles, same landing positions per nest]
Control (S2'): the same with g = cell directly LEFT of h (same row) -- should be a different statement."""
import sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1114c2-Q590')
from common import *
from collections import Counter
def key_of(b): return (frozenset(b.nest[0]['pos']), tuple(tiles_of(b)))
st = Counter(); bad = {}
for key, b0 in census():
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    stack = [b0]; seen = set()
    while stack and len(seen) < 60:
        b = stack.pop(); k = key_of(b)
        if k in seen: continue
        seen.add(k); S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
        for h in cs:
            for nm, g in (('below', (h[0], h[1]-1)), ('left', (h[0]-1, h[1]))):
                if g not in S or g not in corners(S - {h}): continue
                p = copy.deepcopy(b); _, Lh, _ = mig1(p, 0, h, DV['C'])
                q = copy.deepcopy(b); _, Xh, _ = mig1(q, 0, h, DV['B']); _, Lg, _ = mig1(q, 0, g, DV['C'])
                r = copy.deepcopy(p); _, Xg, _ = mig1(r, 0, g, DV['B'])
                s1 = Lh == Lg; s2 = tiles_of(q) == tiles_of(r) and Xh == Xg
                st[(nm, 'n')] += 1; st[(nm, 'S1 landing')] += s1; st[(nm, 'S2 tiling')] += s2
                if not s1 and (nm, 'S1') not in bad: bad[(nm, 'S1')] = (key, sorted(S), h, g, Lh, Lg)
                if not s2 and (nm, 'S2') not in bad: bad[(nm, 'S2')] = (key, sorted(S), h, g)
    st['cells'] += 1
for x in sorted(st, key=str): print(x, st[x])
for k, v in bad.items(): print(k, v)
