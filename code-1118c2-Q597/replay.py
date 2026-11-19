"""Replay helpers on align.pkl (<=5): boards as sets of (kind, frozenset(vertices)); apply steps."""
import sys; sys.path.insert(0, '/home/clio/projects/proofs/code-1118-Q597')
from q597lib import *
def board0(r): return {(k, frozenset(vs)) for k, vs in r['M'].values()}
def apply(Bd, steps):
    Bd = set(Bd)
    for _, g, n in steps:
        assert g <= Bd, 'step not applicable'
        Bd -= g; Bd |= n
    return Bd
def clean(r):
    side = r['Us'][0][0]
    return (side, r['Jy'], r['Jyx'], r['Jx']) if side == 0 else (side, r['Jx'], r['Jxy'], r['Jy'])
