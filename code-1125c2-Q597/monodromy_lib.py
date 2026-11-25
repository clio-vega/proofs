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

