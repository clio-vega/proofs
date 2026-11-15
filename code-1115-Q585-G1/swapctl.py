"""Control for (SW): for distinct outer corners x,y of the A-nest (non-adjacent, both corners now), compare
   x->B ; y->C   with   x->C ; y->B   (assignment swap for a NON-adjacent pair).
If this always holds too, (SW) is a special case of a general 'assignment-irrelevance' and the adjacency is not what matters."""
from instr import *
from local import key_of
from collections import Counter
st = Counter(); ex = []
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
            for y in cs:
                if not x < y: continue
                b1 = copy.deepcopy(b); _, lxB, _ = mig1(b1,0,x,DV['B']); _, lyC, _ = mig1(b1,0,y,DV['C'])
                b2 = copy.deepcopy(b); _, lxC, _ = mig1(b2,0,x,DV['C']); _, lyB, _ = mig1(b2,0,y,DV['B'])
                st['n'] += 1; same = tiles_of(b1) == tiles_of(b2); st['same tiling'] += same
                st['same landing set'] += {lxB, lyC} == {lxC, lyB}
                st['B box equal and C box equal'] += (lxB == lyB and lyC == lxC)
                if not same and len(ex) < 2: ex.append((key, x, y, lxB, lyC, lxC, lyB))
    st['cells'] += 1
for k in sorted(st): print(k, st[k])
for e in ex: print(e)
