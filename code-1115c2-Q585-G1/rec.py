"""Instrumented migration: record each hexagon rotation of a journey.  Engine unchanged (monkeypatch wraps Board.step)."""
import sys, copy
sys.path.insert(0,'/home/clio/projects/proofs/code-1114c2-Q590')
from common import *
import migrate
from fractions import Fraction as Fr
_orig_step = migrate.Board.step
LOG = []
def step_rec(self, R, d):
    before = {t: self.tiles[t] for t in self.tiles}
    R2 = _orig_step(self, R, d)
    if R2 is not None:
        gone = set(before) - set(self.tiles)
        S = frozenset((before[t][0], frozenset(before[t][1])) for t in gone)
        LOG.append((self.last, S, frozenset(self.tiles[R2][1]), frozenset(before[R][1])))
    return R2
migrate.Board.step = step_rec
def mig_rec(b, X, ij, dvec):
    LOG.clear(); t = b.nest[X]['pos'].pop(ij); R, Y, ij2, s = b.journey(t, X, unit_dir(dvec))
    return Y, ij2, list(LOG)
def key_of(b): return (frozenset(b.nest[0]['pos']), tuple(tiles_of(b)))
def states(b0, cap=60):
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    stack = [b0]; seen = set()
    while stack and len(seen) < cap:
        b = stack.pop(); k = key_of(b)
        if k in seen: continue
        seen.add(k); yield b, DV
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
