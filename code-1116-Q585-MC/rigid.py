"""RIGIDITY: within the set of states reachable from one census start (single A->B / A->C migrations), is the state a
function of its nest data (the three position sets)?  Full BFS per start; reports group sizes.  Control: across DIFFERENT
starts of the same cell (lam,mu,nu), the same nest data must occur with different tilings (else rigidity is a property of the
cell, i.e. c=1, and says nothing)."""
from qp import *
from collections import Counter, defaultdict
import time
st = Counter(); t0 = time.time(); percell = defaultdict(dict)
for key, b0 in census():
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    seen = {}; stack = [b0]
    while stack:
        b = stack.pop(); k = fullkey(b)
        if k in seen: continue
        seen[k] = 1
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
    g = Counter(k[0] for k in seen)
    st['starts'] += 1; st['states'] += len(seen); st['nestdata classes'] += len(g); st['classes size>1'] += sum(v > 1 for v in g.values())
    cell = key[:3]
    for k in seen:
        d = percell[cell].setdefault(k[0], set()); d.add(k[1])
    st['max |alpha|'] = max(st['max |alpha|'], len(A['pos']))
for cell, D in percell.items():
    st['cells'] += 1; st['cells where some nestdata has >1 tiling across starts'] += any(len(v) > 1 for v in D.values())
for k in sorted(st, key=str): print(k, st[k])
print('time', round(time.time()-t0))
