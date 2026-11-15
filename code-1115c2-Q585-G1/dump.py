"""Dump every diamond (BC and BB) with full hexagon logs, pickled for classification."""
from rec import *
import pickle
out = []
for key, b0 in census():
    for b, DV in states(b0):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind in ('BC', 'BB'):
                    if kind == 'BB' and x > y: continue
                    b1 = copy.deepcopy(b); X1, l1x, h1x = mig_rec(b1, 0, x, DV[kind[0]]); Y1, l1y, h1y = mig_rec(b1, 0, y, DV[kind[1]])
                    b2 = copy.deepcopy(b); Y2, l2y, h2y = mig_rec(b2, 0, y, DV[kind[1]]); X2, l2x, h2x = mig_rec(b2, 0, x, DV[kind[0]])
                    ok = tiles_of(b1) == tiles_of(b2) and (l1x, l1y) == (l2x, l2y)
                    out.append(dict(key=key, kind=kind, x=x, y=y, ok=ok, l=(l1x, l1y, l2x, l2y),
                                    h1x=h1x, h1y=h1y, h2y=h2y, h2x=h2x, A=frozenset(b.nest[0]['pos'])))
pickle.dump(out, open('/tmp/c2-diamonds.pkl', 'wb'))
from collections import Counter
c = Counter((d['kind'], d['ok'], [s[0] for s in d['h1x']] != [s[0] for s in d['h2x']] or [s[0] for s in d['h1y']] != [s[0] for s in d['h2y']]) for d in out)
print(c)
