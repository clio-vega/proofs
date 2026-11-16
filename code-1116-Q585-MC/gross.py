"""Per journey: does any triangle move in more than one step (gross), and is the net move = the single step's move?
And the per-step direction set by migration type and flip type.  Exact."""
from qp import *
import migrate
from collections import Counter
st = Counter(); CUR = {}
_orig = migrate.Board.step
def step3(self, R, d):
    b4 = qp(self); R2 = _orig(self, R, d)
    if R2 is None: return R2
    af = qp(self)
    for i in (0, 1):
        for t in b4[i]:
            if af[i][t] != b4[i][t]:
                CUR.setdefault((i, t), []).append(tuple(a-c for a, c in zip(af[i][t], b4[i][t])))
    return R2
migrate.Board.step = step3
n = 0
for key, b0 in census():
    for b, DV in states(b0, cap=12):
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); CUR.clear(); mig1(bx, 0, x, DV[t])
                st[(t, 'journeys')] += 1
                st[(t, 'some triangle moved twice')] += any(len(v) > 1 for v in CUR.values())
                for (i, tt), v in CUR.items(): st[(t, 'part', i, 'step vec', v[0])] += 1
    n += 1
    if n >= 50: break
for k in sorted(st, key=str): print(k, st[k])
