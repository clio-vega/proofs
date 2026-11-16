"""Classify every rotated hexagon by the parts of its two triangles.  Claim (to prove): type a = two triangles of the SAME
type (both 0 or both 1) forming a centrally symmetric lozenge, the other part a unit segment; type b = one 0-tri + one 1-tri.
Also test the per-step displacement law on (q,p)."""
from qp import *
import migrate
from collections import Counter
st = Counter(); n = 0
_orig = migrate.Board.step
def step2(self, R, d):
    before = qp(self); bt = dict(self.tiles); R2 = _orig(self, R, d)
    if R2 is None: return R2
    gone = [bt[t] for t in bt if t not in self.tiles]
    tris = [vs for k, vs in gone if k == 'tri']
    kinds = []
    for vs in tris:
        p1 = {v[2:] for v in vs}; kinds.append(0 if len(p1) == 1 else 1)
    after = qp(self)
    dq = Counter(tuple(a-c for a, c in zip(after[0][t], before[0][t])) for t in before[0] if after[0][t] != before[0][t])
    dp = Counter(tuple(a-c for a, c in zip(after[1][t], before[1][t])) for t in before[1] if after[1][t] != before[1][t])
    st[(self.last, tuple(sorted(kinds)), 'moved q', sum(dq.values()), 'moved p', sum(dp.values()))] += 1
    if self.last == 'a' and kinds[0] == kinds[1]:
        i = kinds[0]; sl = (slice(0, 2) if i == 0 else slice(2, 4)); ot = (slice(2, 4) if i == 0 else slice(0, 2))
        A0 = {v[sl] for v in tris[0]}; A1 = {v[sl] for v in tris[1]}
        st[('a', 'lozenge: share an edge')] += len(A0 & A1) == 2
        o0 = {v[ot] for v in tris[0]}.pop(); o1 = {v[ot] for v in tris[1]}.pop()
        e = tuple(a-c for a, c in zip(o1, o0))
        D = dq if i == 0 else dp
        st[("a", "offset", e, "D", tuple(sorted(D.items())))] += 1
        st[('a', 'two triangles swap (dq = +-e)')] += sorted(D) == sorted([e, tuple(-a for a in e)])
    return R2
migrate.Board.step = step2
for key, b0 in census():
    for b, DV in states(b0, cap=10):
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t])
    n += 1
    if n >= 6: break
for k in sorted(st, key=str): print(k, st[k])
