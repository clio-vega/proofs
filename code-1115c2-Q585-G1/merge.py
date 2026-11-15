"""SW mechanism: for column pairs (h above g): compare g's C-journey after h->B with h's SOLO C-journey (and g's B-journey after
h->C with h's solo B-journey).  'merge index' = smallest i such that g's steps from i on equal a suffix of h's solo steps
(hexagon tile-sets).  Small states only (first 40 cells, 30 states)."""
from rec import *
from collections import Counter
st = Counter(); n = 0
def merge(G, H):
    for i in range(len(G)+1):
        tail = G[i:]
        if len(tail) <= len(H) and H[len(H)-len(tail):] == tail: return i, len(H) - len(tail)
for key, b0 in census():
    n += 1
    if n > 40: break
    for b, DV in states(b0, 30):
        S = set(b.nest[0]['pos'])
        for h in corners(S):
            g = (h[0], h[1]-1)
            if g not in S or g not in corners(S - {h}): continue
            for first, second in (('B', 'C'), ('C', 'B')):
                p = copy.deepcopy(b); _, _, Hs = mig_rec(p, 0, h, DV[second])          # h solo to `second`
                q = copy.deepcopy(b); mig_rec(q, 0, h, DV[first]); _, _, Gs = mig_rec(q, 0, g, DV[second])
                Hh = [s[1] for s in Hs]; Gg = [s[1] for s in Gs]
                i, j = merge(Gg, Hh)
                st[(second, 'g-prefix before merge', i, 'h-prefix skipped', j, 'merged tail empty', i == len(Gg))] += 1
for k in sorted(st, key=str): print(k, st[k])
