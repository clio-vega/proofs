"""Swap test.  From every reachable state (mixed single migrations), for each corner x=(i,j) of A:
 column swap: y=(i,j-1) in A  :  x->B ; y->C   vs   x->C ; y->B     (does the pair of landings / tiling agree?)
 row swap:    y=(i-1,j) in A  :  same two paths.
Report: same tiling; same SET of landings; landing of the C-move equal; landing of the B-move equal."""
from instr import *
from collections import Counter
st = Counter(); ex = []
def key_of(b): return (frozenset(b.nest[0]['pos']), tuple(tiles_of(b)))
for key, b0 in census():
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    stack=[b0]; seen=set(); ns=0
    while stack and ns < 60:
        b = stack.pop(); k = key_of(b)
        if k in seen: continue
        seen.add(k); ns += 1; S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
        for x in cs:
            i, j = x
            for kind, y in (('col', (i, j-1)), ('row', (i-1, j))):
                if y not in S: continue
                if y not in corners(S - {x}): st[kind+' skipped: y not a corner after x leaves'] += 1; continue
                try:
                  b1 = copy.deepcopy(b); _, lxB, _ = mig1(b1,0,x,DV['B']); _, lyC, _ = mig1(b1,0,y,DV['C'])
                  b2 = copy.deepcopy(b); _, lxC, _ = mig1(b2,0,x,DV['C']); _, lyB, _ = mig1(b2,0,y,DV['B'])
                except RuntimeError as e:
                  st[kind+' ERR '+str(e)[:25]] += 1; continue
                st[kind+' n'] += 1
                st[kind+' same tiling'] += tiles_of(b1) == tiles_of(b2)
                st[kind+' C-landing equal'] += lyC == lxC
                st[kind+' B-landing equal'] += lxB == lyB
                if tiles_of(b1) != tiles_of(b2) and len(ex) < 4: ex.append((kind, key, x, y, lxB, lyC, lxC, lyB))
for k in sorted(st): print(k, st[k])
for e in ex: print(e)
