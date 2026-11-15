"""Instrumented single migrations: log each hexagon step (centre c, type, rotation k, tile kinds+vertex sets before/after)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1114c2-Q590')
from common import *
import migrate
from fractions import Fraction as Fr
LOG = []
_orig_step = migrate.Board.step
def step_logged(self, R, d):
    before = {t: self.tiles[t] for t in self.tiles}
    newR = _orig_step(self, R, d)
    if newR is not None:
        after = {t: self.tiles[t] for t in self.tiles}
        gone = {t: before[t] for t in before if t not in after}
        new = {t: after[t] for t in after if t not in before}
        LOG.append(dict(typ=self.last, gone=sorted(gone.values()), new=sorted(new.values()), R_before=before[R], R_after=after[newR]))
    return newR
migrate.Board.step = step_logged
def mig_logged(b, X, ij, dvec):
    LOG.clear(); t = b.nest[X]['pos'].pop(ij); R, Y, ij2, s = b.journey(t, X, unit_dir(dvec))
    return Y, ij2, list(LOG)
def hexkey(step):
    """region of a hexagon step = frozenset of its 6 boundary vertices (via all tile vertices minus interior)"""
    return frozenset(v for k, vs in step['gone'] for v in vs)
