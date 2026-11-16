"""Step-2 instrument (pre-registered control: B/B).  Census, states() cap 60 (as c1/c2).  For each ordered corner pair:
 BC  : x->B ; y->C  vs  y->C ; x->B          (MC)
 BB  : x->B ; y->B  vs  y->B ; x->B, x<y     (control: must FAIL the invariant wherever the tilings differ)
 SW  : x->B ; y->C  vs  x->C ; y->B, y=(c,r-1) below x=(c,r), y a corner after x leaves
Invariants (exact integers):
 U  = union law: moved-triangle sets' union equal in both composites
 N  = no-repeat: no triangle moves twice inside one journey
 D  = direction law: every 0-tri step in {(0,-1),(1,-1)} for B, {(-1,0),(0,-1)} for C; 1-tri {(-1,1),(0,1)} / {(1,0),(0,1)}
 Q  = final (q,p) equal
 T  = final tilings equal (ground truth)
Coverage: C-after-B second journeys that change route and have EMPTY common hexagon suffix with solo ('land inside wake')."""
from qp import *
import migrate
from collections import Counter
DIRS = {('B', 0): {(0, -1), (1, -1)}, ('C', 0): {(-1, 0), (0, -1)}, ('B', 1): {(-1, 1), (0, 1)}, ('C', 1): {(1, 0), (0, 1)}}
_orig = migrate.Board.step; CUR = {}
def stepI(self, R, d):
    b4 = qp(self); R2 = _orig(self, R, d)
    if R2 is None: return R2
    af = qp(self)
    for i in (0, 1):
        for t in b4[i]:
            if af[i][t] != b4[i][t]: CUR.setdefault((i, t), []).append(tuple(a-c for a, c in zip(af[i][t], b4[i][t])))
    return R2
migrate.Board.step = stepI
def move(b, x, T, DV):
    CUR.clear(); Y, ij, h = mig_rec(b, 0, x, DV[T]); mv = dict(CUR)
    N = all(len(v) == 1 for v in mv.values()); D = all(v[0] in DIRS[(T, i)] for (i, t), v in mv.items())
    return dict(land=(Y, ij), moved=set(mv), N=N, D=D, hex=[s[1] for s in h])
st = Counter()
def run(kind, b, x, y, DV, s1, s2):
    b1 = copy.deepcopy(b); m1 = [move(b1, x, s1[0], DV), move(b1, y, s1[1], DV)]
    b2 = copy.deepcopy(b)
    if kind == 'SW': m2 = [move(b2, x, s2[0], DV), move(b2, y, s2[1], DV)]
    else: m2 = [move(b2, y, s1[1], DV), move(b2, x, s1[0], DV)]
    T = tiles_of(b1) == tiles_of(b2); Q = qp(b1) == qp(b2)
    U = (m1[0]['moved'] | m1[1]['moved']) == (m2[0]['moved'] | m2[1]['moved'])
    L = sorted(m['land'] for m in m1) == sorted(m['land'] for m in m2)
    st[(kind, 'n')] += 1; st[(kind, 'T tiling equal')] += T; st[(kind, 'Q (q,p) equal')] += Q; st[(kind, 'Q==T agree')] += Q == T
    st[(kind, 'U union law')] += U; st[(kind, 'U holds where T fails')] += U and not T
    st[(kind, 'landing multiset equal')] += L
    st[(kind, 'N,D all four journeys')] += all(m['N'] and m['D'] for m in m1 + m2)
    return b1, b2, m1, m2
for key, b0 in census():
    for b, DV in states(b0, cap=60):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                b1, b2, m1, m2 = run('BC', b, x, y, DV, 'BC', None)
                # coverage: C-after-B (second journey of order 1 is y->C) vs solo y->C at b
                bs = copy.deepcopy(b); solo = move(bs, y, 'C', DV)
                H2, HA = m1[1]['hex'], solo['hex']
                if H2 != HA:
                    st[('BC', 'C-after-B route changed')] += 1
                    st[('BC', 'C-after-B lands inside wake (empty common suffix)')] += H2[-1] != HA[-1]
                if x < y: run('BB', b, x, y, DV, 'BB', None)
            c, r = x
            if r >= 1 and (c, r-1) in cs or (r >= 1 and (c, r-1) in corners(set(b.nest[0]['pos']) - {x}) and (c+1, r-1) not in b.nest[0]['pos']):
                y = (c, r-1)
                if y in corners(set(b.nest[0]['pos']) - {x}): run('SW', b, x, y, DV, 'BC', 'CB')
    st['cells'] += 1
for k in sorted(st, key=str): print(k, st[k], flush=True)
