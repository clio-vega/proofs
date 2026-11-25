"""Step 0 (PROVE 2026-11-25 c2): the commutor-loop / monodromy test on B/B pairs, with B/C control.

Definitions (fixed BEFORE running):
  move x->T : pop nest-A corner x, journey toward T (dvec DV[T]), lands at nest-T cell c_x.
  undo(T, c): pop nest-T cell c (must be a corner of nest T's ideal), journey with dvec -DV[T]; record where it lands.
  INV control : M -> Mx -> undo(T, c_x).  Returns M?  (single-migration invertibility, Purbhoo Thm 2 for one rhombus)
  LIFO loop   : M -> Mx -> Mxy -> undo(y's cell) -> undo(x's cell).  Returns M?
  COMMUTOR loop (the monodromy test): M -> Mx -> Mxy, then undo x's cell FIRST, then y's cell.  Returns M?
     For a commuting pair the two orders are interchangeable, so this is a genuine loop around the square.
Population: census(), states(cap=60) = the 5194-state size-5 census; ordered pairs of distinct nest-A corners.
"""
import sys, copy
sys.path.insert(0, '/home/clio/projects/proofs/code-1115c2-Q585-G1')
from rec import *
from collections import Counter

def key(b): return (tuple(frozenset(b.nest[i]['pos']) for i in range(3)), tuple(tiles_of(b)))

def move(b, X, ij, dvec):
    t = b.nest[X]['pos'].pop(ij)
    R, Y, ij2, s = b.journey(t, X, unit_dir(dvec))
    return Y, ij2

def undo(b, X, c, dvec):
    if c not in corners(set(b.nest[X]['pos'])): return None
    try:
        return move(b, X, c, tuple(-v for v in dvec))
    except RuntimeError as e:
        return ('err', str(e)[:30])

NEST = {'B': 1, 'C': 2}
tab = Counter()
for kk, b0 in census():
    for b, DV in states(b0, cap=60):
        k0 = key(b)
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for (tx, ty) in (('B', 'B'), ('B', 'C')):
                    if tx == ty == 'B' and not x < y: continue
                    cls = tx + ty
                    # INV control on x alone
                    bx = copy.deepcopy(b); Yx, cx = move(bx, 0, x, DV[tx])
                    bi = copy.deepcopy(bx); r = undo(bi, NEST[tx], cx, DV[tx])
                    tab[(cls, 'INV', 'ok' if (r is not None and key(bi) == k0) else ('notcorner' if r is None else 'FAIL'))] += 1
                    bxy = copy.deepcopy(bx); Yy, cy = move(bxy, 0, y, DV[ty])
                    by = copy.deepcopy(b); move(by, 0, y, DV[ty]); byx = copy.deepcopy(by); move(byx, 0, x, DV[tx])
                    comm = key(bxy) == key(byx)
                    tag = cls + ('-comm' if comm else '-noncomm')
                    tab[(tag, 'count')] += 1
                    # LIFO loop
                    bl = copy.deepcopy(bxy); r1 = undo(bl, NEST[ty], cy, DV[ty]); r2 = undo(bl, NEST[tx], cx, DV[tx]) if r1 is not None and r1[0] != 'err' else None
                    tab[(tag, 'LIFO', 'ok' if (r2 is not None and r2[0] != 'err' and key(bl) == k0) else ('undefined' if r2 is None else 'FAIL'))] += 1
                    # COMMUTOR loop: undo x's cell first
                    bc = copy.deepcopy(bxy); r1 = undo(bc, NEST[tx], cx, DV[tx])
                    if r1 is None or r1[0] == 'err':
                        tab[(tag, 'COMMUTOR', 'x-cell not a corner' if r1 is None else 'err')] += 1; continue
                    r2 = undo(bc, NEST[ty], cy, DV[ty])
                    if r2 is None or r2[0] == 'err':
                        tab[(tag, 'COMMUTOR', 'y-cell not a corner after' if r2 is None else 'err')] += 1; continue
                    tab[(tag, 'COMMUTOR', 'closes' if key(bc) == k0 else 'OPEN')] += 1
for k in sorted(tab): print(k, tab[k])
