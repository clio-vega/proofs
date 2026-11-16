"""Q594 /assumptions suspects A7, A8, A16, A17 on the 1558 B/C (MC) pairs and 779 B/B pairs (same sampler as cactus.py).
For B/C with solo wakes W_B = X(x->B @ s), W_C = X(y->C @ s), and composite wakes W'_C = X(y->C @ after x->B), W'_B = X(x->B @ after y->C):
  I = W_B & W_C (size distribution);  dC = W'_C ^ W_C, dB = W'_B ^ W_B (symmetric differences);
  A8: dC, dB  subset of I?   A7: |dC|,|dB| distribution.
Also the per-triangle displacement in each order (exact integer vectors) on I: which journeys push t, with which vector."""
import sys, copy
sys.path.insert(0, '/home/clio/projects/proofs/code-1116-Q585-MC')
from qp import *
import migrate
from collections import Counter
_orig = migrate.Board.step
CUR = {}
def stepI(self, R, d):
    b4 = qp(self); R2 = _orig(self, R, d)
    if R2 is None: return R2
    af = qp(self)
    for i in (0, 1):
        for t in b4[i]:
            if af[i][t] != b4[i][t]: CUR[(i, t)] = CUR.get((i, t), ()) + (tuple(a-c for a, c in zip(af[i][t], b4[i][t])),)
    return R2
migrate.Board.step = stepI
def mv(b, x, T, DV):
    bb = copy.deepcopy(b); CUR.clear(); mig_rec(bb, 0, x, DV[T]); return dict(CUR), bb
st = Counter(); pat = Counter()
for key, b0 in census():
    for b, DV in states(b0, cap=60):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind, (tx, ty) in (('BC', 'BC'), ('BB', 'BB')):
                    if kind == 'BB' and not x < y: continue
                    WB, b1 = mv(b, x, tx, DV); WC, b2 = mv(b, y, ty, DV)
                    WCp, b12 = mv(b1, y, ty, DV); WBp, b21 = mv(b2, x, tx, DV)
                    com = tiles_of(b12) == tiles_of(b21)
                    I = set(WB) & set(WC)
                    if not I:
                        st[(kind, 'disjoint; composite wakes == solo')] += (WCp == WC and WBp == WB)
                        st[(kind, 'disjoint n')] += 1
                        continue
                    dC = set(WCp) ^ set(WC); dB = set(WBp) ^ set(WB)
                    st[(kind, com, '|I|', len(I))] += 1
                    st[(kind, com, '|dC|,|dB|', len(dC), len(dB))] += 1
                    st[(kind, com, 'dC,dB subset I')] += (dC <= I and dB <= I)
                    st[(kind, com, 'vectors agree off I')] += all(WCp.get(t) == WC.get(t) for t in set(WCp) | set(WC) if t not in I) and \
                        all(WBp.get(t) == WB.get(t) for t in set(WBp) | set(WB) if t not in I)
                    if kind == 'BC':
                        for t in I:
                            # order 1: x->B then y->C ; order 2: y->C then x->B
                            pat[(t[0], 'o1', WB.get(t), WCp.get(t), 'o2', WC.get(t), WBp.get(t))] += 1
for k in sorted(st, key=str): print(k, st[k])
print('--- B/C patterns on I: (type, o1: B-vec, C-vec after B, o2: C-vec, B-vec after C) ---')
for k, v in sorted(pat.items(), key=lambda kv: -kv[1]): print(v, k)
