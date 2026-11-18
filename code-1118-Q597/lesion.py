"""Lesion control: replace s by b or c in (E). Must fail. Uses the q/p step vectors of 11-17c2."""
from collections import Counter
from q597lib import *
V = {0: {'s': (0, -1), 'b': (1, -1), 'c': (-1, 0)}, 1: {'s': (0, 1), 'b': (-1, 1), 'c': (1, 0)}}
st = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side = r['Us'][0][0]; I = set(addd(r['Yx'], neg(r['Y'])))
    for nm, vec in V[side].items():
        ok = all(tuple(a+b for a, b in zip(r['X'].get(k, (0, 0)), r['Yx'].get(k, (0, 0)))) == vec for k in I)
        st[(nm, ok)] += 1
print(dict(st))
