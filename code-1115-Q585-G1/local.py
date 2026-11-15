"""Locality of the B/C interaction.  For each diamond (x->B, y->C), compare x's journey in P1 (x first) and P2 (y first)
as lists of (hexagon region, type).  Measure: common prefix, common suffix, and whether every differing x-step region
meets some region touched by y's journey (in either path).  Same for y.  Also the B/B control with the same measure."""
from instr import *
from collections import Counter
st = Counter(); ex = []
def key_of(b): return (frozenset(b.nest[0]['pos']), tuple(tiles_of(b)))
def sig(steps): return [(hexkey(s), s['typ']) for s in steps]
def pre(a, b):
    k = 0
    while k < min(len(a), len(b)) and a[k] == b[k]: k += 1
    return k
def analyse(kind, s1x, s2x, s1y, s2y):
    a1, a2, b1, b2 = sig(s1x), sig(s2x), sig(s1y), sig(s2y)
    if a1 == a2 and b1 == b2: st[kind+' identical journeys'] += 1; return
    st[kind+' route changes'] += 1
    yreg = set(r for r, _ in b1 + b2); xreg = set(r for r, _ in a1 + a2)
    dx = [r for r, _ in a1 + a2 if (r, _) not in set(a1) & set(a2)]
    dy = [r for r, _ in b1 + b2 if (r, _) not in set(b1) & set(b2)]
    st[kind+' x-diff regions all meet y-regions'] += all(any(r & q for q in yreg) for r in dx)
    st[kind+' x-diff regions share an edge (2 verts) w/ y-regions'] += all(any(len(r & q) >= 2 for q in yreg) for r in dx)
    st[kind+' x,y journeys share a full hexagon region'] += bool(set(r for r,_ in a1+a2) & yreg)
    st[(kind, 'prefix x', pre(a1, a2) > 0, 'suffix x', pre(a1[::-1], a2[::-1]) > 0)] += 1
for key, b0 in (census() if __name__ == "__main__" else iter(())):
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    stack=[b0]; seen=set(); ns=0
    while stack and ns < 60:
        b = stack.pop(); k = key_of(b)
        if k in seen: continue
        seen.add(k); ns += 1; cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind in ('BC', 'BB'):
                    if kind == 'BB' and x > y: continue
                    b1 = copy.deepcopy(b); _, l1x, s1x = mig_logged(b1,0,x,DV[kind[0]]); _, l1y, s1y = mig_logged(b1,0,y,DV[kind[1]])
                    b2 = copy.deepcopy(b); _, l2y, s2y = mig_logged(b2,0,y,DV[kind[1]]); _, l2x, s2x = mig_logged(b2,0,x,DV[kind[0]])
                    st[kind+' diamonds'] += 1
                    analyse(kind, s1x, s2x, s1y, s2y)
if __name__ == "__main__":
    for k in sorted(st, key=str): print(k, st[k])
