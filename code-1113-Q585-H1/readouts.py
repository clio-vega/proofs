"""Table-1 / Prop 5.1 readouts in the engine.
readout(board, X, sgn, dirn): give nest X the orientation sgn*(E_X,N_X) (E_X,N_X = nest (u,w)), form the unique flock,
migrate in compass direction dirn (relative to that orientation). Returns (dst nest, tableau dict in induced frame)."""
import sys, copy; sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from flock import *

def unique_flock(board, X, E_, N_):
    info = board.nest[X]
    bx = {box_of(board, X, board.tiles[t][1], E_, N_): t for t in info['pos'].values()}
    if not bx: return {}
    shape = tuple(sorted(Counter(i for i, j in info['pos']).values(), reverse=True))  # cols or rows: try both
    for sh in {shape, tuple(sorted(Counter(j for i, j in info['pos']).values(), reverse=True))}:
        L = lr_tableaux(set(bx), sh)
        if len(L) == 1: return {bx[b]: v for b, v in L[0].items()}
    raise RuntimeError('no unique flock')

def readout(board, X, sgn, dirn):
    b = copy.deepcopy(board); u, w = b.nest[X]['u'], b.nest[X]['w']
    E_ = tuple(sgn*a for a in u); N_ = tuple(sgn*a for a in w)
    fl = unique_flock(b, X, E_, N_)
    if not fl: return None, {}
    dst, out, fr, boxes, st = migrate_flock(b, X, fl, E_, N_, dirn)
    return dst, boxes

def norm(T):
    if not T: return ()
    c0 = min(c for c, r in T); r0 = min(r for c, r in T)
    return tuple(sorted(((c-c0, r-r0), v) for (c, r), v in T.items()))
