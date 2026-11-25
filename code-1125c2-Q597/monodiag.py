"""Diagnose the OPEN commutor loops of monodromy.py (B/B non-commuting).
For each OPEN loop record: (a) does M' have the same nest-A cell set as M?  (b) same tiling?
(c) where do the two undone rhombi land (which nest-A cells), vs the original corners x,y?
(d) is y's solo landing cell in M equal to y's landing cell in Mx (c_y)?  Same for x in My."""
import sys, copy
sys.path.insert(0, '/home/clio/projects/proofs/code-1115c2-Q585-G1')
from rec import *
from collections import Counter
from monodromy_lib import key, move, undo
tab = Counter()
for kk, b0 in census():
    for b, DV in states(b0, cap=60):
        k0 = key(b); cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if not x < y: continue
                bx = copy.deepcopy(b); _, cx = move(bx, 0, x, DV['B'])
                bxy = copy.deepcopy(bx); _, cy = move(bxy, 0, y, DV['B'])
                by = copy.deepcopy(b); _, cy0 = move(by, 0, y, DV['B'])
                byx = copy.deepcopy(by); _, cx1 = move(byx, 0, x, DV['B'])
                if key(bxy) == key(byx): continue
                tab['noncomm'] += 1
                tab[('cells', 'cy==cy0', cy == cy0, 'cx==cx1', cx == cx1, 'samecellset', {cx, cy} == {cy0, cx1})] += 1
                bc = copy.deepcopy(bxy); r1 = undo(bc, 1, cx, DV['B'])
                if r1 is None: continue
                r2 = undo(bc, 1, cy, DV['B'])
                if r2 is None: continue
                k1 = key(bc)
                if k1 == k0: tab['closes'] += 1; continue
                tab[('OPEN', 'nestA same', k1[0][0] == k0[0][0], 'tiles same', k1[1] == k0[1], 'land', (r1[1] == x, r2[1] == y), 'swapped', (r1[1] == y, r2[1] == x))] += 1
for k in sorted(tab, key=str): print(k, tab[k])
