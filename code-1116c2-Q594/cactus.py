"""Q594 step 1: the cactus-shape test on B/B pairs (pre-registered; defined BEFORE running).

Population: census (134 cells), states() cap 60 -- identical sampler to code-1116-Q585-MC/verify.py, which gave 779 B/B pairs
(x<y distinct A-corners), 339 commuting (tilings equal), 440 not.

Window index (defined here, before the run).  For a move m = (x -> B) at state sigma, run SOLO from sigma:
  W_tri(m) = set of triangles moved by the journey, a triangle named by its conserved shape (Cor. cons, 11-16):
             ('0', pi_0 tau) or ('1', pi_1 sigma).  This is the 11-16 'moved set' X(m).
  W_hex(m) = set of hexagons rotated by the journey, a hexagon named by the frozenset of its Z^4 corner points
             (union of the vertex sets of the four tiles that rotate).  Purbhoo's 'route'.
Shape of a pair (W(x), W(y)): 'disjoint' if W(x) & W(y) empty; 'nested' if one contains the other (and they meet);
'overlap' otherwise.  Also a third index, purely tableau-side, cell-interval:
  W_row(m) = the interval [min, max] of A-nest rows (second coordinate) hit by ... -- NOT defined: a single corner has a
  one-point row set, so the row index is trivial for single migrations.  Recorded as not used.

Falsifier (pre-registered): one NON-commuting B/B pair with a DISJOINT window pair kills 'cactus shape' for B/B.
Also print the same table for B/C (MC, all commute) for comparison."""
import sys, copy
sys.path.insert(0, '/home/clio/projects/proofs/code-1116-Q585-MC')
from qp import *
import migrate
from collections import Counter

_orig = migrate.Board.step
CUR = {}; HEX = []
def stepI(self, R, d):
    b4 = qp(self); before = {t: self.tiles[t] for t in self.tiles}
    R2 = _orig(self, R, d)
    if R2 is None: return R2
    af = qp(self)
    gone = set(before) - set(self.tiles)
    corners_ = frozenset(v for t in gone for v in before[t][1])
    HEX.append(corners_)
    for i in (0, 1):
        for t in b4[i]:
            if af[i][t] != b4[i][t]: CUR.setdefault((i, t), []).append(1)
    return R2
migrate.Board.step = stepI

def solo(b, x, T, DV):
    bb = copy.deepcopy(b); CUR.clear(); HEX.clear()
    mig_rec(bb, 0, x, DV[T])
    return set(CUR), set(HEX), bb

def shape(U, V):
    if not (U & V): return 'disjoint'
    if U <= V or V <= U: return 'nested'
    return 'overlap'

def run_pair(b, x, y, DV, tx, ty):
    _, _, b1 = solo(b, x, tx, DV); _, _, b1 = solo(b1, y, ty, DV)
    _, _, b2 = solo(b, y, ty, DV); _, _, b2 = solo(b2, x, tx, DV)
    return tiles_of(b1) == tiles_of(b2)

tab = Counter(); n = Counter(); hexsize = Counter()
for key, b0 in census():
    for b, DV in states(b0, cap=60):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind in ('BB', 'BC'):
                    if kind == 'BB' and not x < y: continue
                    tx, ty = ('B', 'B') if kind == 'BB' else ('B', 'C')
                    Tx, Hx, _ = solo(b, x, tx, DV); Ty, Hy, _ = solo(b, y, ty, DV)
                    com = run_pair(b, x, y, DV, tx, ty)
                    n[kind] += 1
                    tab[(kind, 'tri', shape(Tx, Ty), com)] += 1
                    tab[(kind, 'hex', shape(Hx, Hy), com)] += 1
                    if not Tx or not Ty: tab[(kind, 'EMPTY tri window', com)] += 1

print('n', dict(n))
for kind in ('BB', 'BC'):
    for idx in ('tri', 'hex'):
        print(f'== {kind} index={idx}:  shape | commute | not')
        for s in ('disjoint', 'nested', 'overlap'):
            print(f'   {s:9s} | {tab[(kind, idx, s, True)]:5d} | {tab[(kind, idx, s, False)]:5d}')
    print(kind, 'empty tri windows (commute,not):', tab[(kind, 'EMPTY tri window', True)], tab[(kind, 'EMPTY tri window', False)])
fals = tab[('BB', 'tri', 'disjoint', False)] + tab[('BB', 'hex', 'disjoint', False)]
print('FALSIFIER (non-commuting B/B on disjoint window, either index):', fals)
