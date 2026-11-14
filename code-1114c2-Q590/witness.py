"""Smallest explicit diamond failure (->B, direction N_A) for the write-up."""
from common import *
for key, b0 in census():
    A = b0.nest[0]; cells = set(A['pos']); cs = corners(cells)
    if len(cs) < 2: continue
    found = False
    for i in range(len(cs)):
        for j in range(i+1, len(cs)):
            x, y = cs[i], cs[j]; res = {}
            for o in ((x, y), (y, x)):
                b = copy.deepcopy(b0); land = []
                for z in o: Y, ij, s = mig1(b, 0, z, A['w']); land.append((z, ij, dict(s)))
                res[o] = (sorted(b.nest[1]['pos']), land, tiles_of(b))
            if res[(x, y)][2] != res[(y, x)][2]:
                print('cell', key, 'alpha cells', sorted(cells), 'beta', sorted(b0.nest[1]['pos']))
                for o, r in res.items(): print('  order', o, '-> B nest', r[0], ' landings', r[1])
                found = True; break
        if found: break
    if found: break
