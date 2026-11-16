"""RIGID': within one orbit (full BFS from a census start), is the state a function of (nest data, set of A-positions sent to B)?
Key built from the move word.  Control: plain nest data (RIGID) is known to fail in 2 classes; RIGID' must be tested on the
same BFS.  Non-vacuity: count classes of RIGID' reached by >1 distinct move word."""
from qp import *
from collections import Counter, defaultdict
import time
st = Counter(); t0 = time.time()
for key, b0 in census():
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    seen = {}; stack = [(b0, frozenset(), ())]
    words = defaultdict(set)
    while stack:
        b, toB, w = stack.pop(); k = fullkey(b)
        words[(k[0], toB)].add(w)
        if (k, toB) in seen: continue
        seen[(k, toB)] = 1
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append((bx, toB | {x} if t == 'B' else toB, w + ((x, t),)))
    g = defaultdict(set); g0 = defaultdict(set)
    for (k, toB) in seen: g[(k[0], toB)].add(k[1]); g0[k[0]].add(k[1])
    st['starts'] += 1; st["RIGID' classes"] += len(g); st["RIGID' classes with >1 tiling"] += sum(len(v) > 1 for v in g.values())
    st["RIGID' classes reached by >1 move word"] += sum(len(words[c]) > 1 for c in g)
    st['RIGID classes with >1 tiling (control: expect 2 total)'] += sum(len(v) > 1 for v in g0.values())
for k in sorted(st, key=str): print(k, st[k], flush=True)
print('time', round(time.time()-t0))
