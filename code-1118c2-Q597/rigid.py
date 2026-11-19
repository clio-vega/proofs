"""Control for hexa: within each census entry, is a state determined by its nest data (the three position SETS)?
If yes, 'closes after a third move' reduces to 'the nest shapes agree', which is cell combinatorics, not journeys."""
import sys, os
os.chdir('/home/clio/projects/proofs/code-1118c2-Q597')
sys.argv = [sys.argv[0], '60']
exec(open('/home/clio/projects/proofs/code-1118-Q597/widen.py').read().split('st = Counter()')[0])
from collections import defaultdict
st = Counter()
for key, b0 in census():
    byn = defaultdict(set)
    for b, DV in states(b0, cap=CAP):
        byn[tuple(frozenset(b.nest[i]['pos']) for i in range(3))].add(tuple(tiles_of(b)))
    for k, v in byn.items(): st[('#states per nest datum', min(len(v), 5))] += 1
for k in sorted(st, key=str): print(st[k], k)
