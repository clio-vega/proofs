"""Locate the rigidity failures: BFS with parent pointers (move words) and print the two words reaching each bad class."""
from qp import *
from collections import Counter, defaultdict
for key, b0 in census():
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    seen = {}; stack = [(b0, ())]
    while stack:
        b, w = stack.pop(); k = fullkey(b)
        if k in seen: continue
        seen[k] = w
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); Y, ij, s = mig1(bx, 0, x, DV[t]); stack.append((bx, w + ((x, t, ij),)))
    g = defaultdict(list)
    for k, w in seen.items(): g[k[0]].append(w)
    for D, ws in g.items():
        if len(ws) > 1:
            print(key, 'alpha', sorted(A['pos']), '\n nestdata', [sorted(s) for s in D])
            for w in ws: print('  word', w)
