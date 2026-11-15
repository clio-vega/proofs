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
        seen.add(k); ns += 1; cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for t in 'BC':
                bx = copy.deepcopy(b); mig_logged(bx, 0, x, DV[t]); stack.append(bx)
        for x in cs:
            for y in cs:
                if x == y: continue
                b1 = copy.deepcopy(b); _, l1x, s1x = mig_logged(b1,0,x,DV['B']); _, l1y, s1y = mig_logged(b1,0,y,DV['C'])
                b2 = copy.deepcopy(b); _, l2y, s2y = mig_logged(b2,0,y,DV['C']); _, l2x, s2x = mig_logged(b2,0,x,DV['B'])
                hx1=[hexkey(s) for s in s1x]; hx2=[hexkey(s) for s in s2x]; hy1=[hexkey(s) for s in s1y]; hy2=[hexkey(s) for s in s2y]
                inter = [ (a['typ'],) for a in s1x] != [(a['typ'],) for a in s2x] or [a['typ'] for a in s1y] != [a['typ'] for a in s2y]
                if not inter: continue
                st['inter']+=1
                st[('lens', len(s1x),len(s2x),len(s1y),len(s2y))]+=1
                if len(ex)<3: ex.append((key,x,y,[a['typ'] for a in s1x],[a['typ'] for a in s2x],[a['typ'] for a in s1y],[a['typ'] for a in s2y]))
print(st); [print(e) for e in ex]
