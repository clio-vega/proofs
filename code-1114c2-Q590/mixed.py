"""Mixed-direction diamonds.  R = (ii) removal order on alpha (top row first, E->W), D = (iii) order (depth 1 first, E->W).
Vertex (i,j): R[:i] sent A->B (dir N_A) and D[:j] sent A->C (dir E_A), defined when disjoint.
Diamond test: path (i,j)->(i+1,j)->(i+1,j+1) vs (i,j)->(i,j+1)->(i+1,j+1): same tiling? same landing boxes?
Positive control: diamond (0,0) with R[0]==D[0] is excluded; plus the brief's zero-bit guard: count diamonds whose two moves
are hexagon-adjacent (both journeys pass >=1 hexagon)."""
from common import *
from collections import Counter
st = Counter(); first = []
def orders(cells):
    conj = Counter(c for c, r in cells)
    R = sorted(cells, key=lambda x: (-x[1], -x[0]))
    D = sorted(cells, key=lambda x: (conj[x[0]] - x[1], -x[0]))
    return R, D
def apply(b0, seq, uA, wA):
    b = copy.deepcopy(b0); land = []
    for kind, x in seq:
        Y, ij, s = mig1(b, 0, x, wA if kind == 'R' else uA); land.append((kind, x, Y, ij))
    return b, land
for key, b0 in census():
    A = b0.nest[0]; uA, wA = A['u'], A['w']; cells = set(A['pos']); R, D = orders(cells); m = len(R)
    cellok = True
    for i in range(m):
        for j in range(m):
            if set(R[:i+1]) & set(D[:j+1]): continue
            base = [('R', x) for x in R[:i]] + [('C', x) for x in D[:j]]
            p1 = base + [('R', R[i]), ('C', D[j])]; p2 = base + [('C', D[j]), ('R', R[i])]
            try:
                b1, l1 = apply(b0, p1, uA, wA); b2, l2 = apply(b0, p2, uA, wA)
            except Exception as e:
                st['ERR '+str(e)[:30]] += 1; continue
            st['diamonds'] += 1
            same = tiles_of(b1) == tiles_of(b2); st['same tiling'] += same
            lb1 = {x: (Y, ij) for k, x, Y, ij in l1[-2:]}; lb2 = {x: (Y, ij) for k, x, Y, ij in l2[-2:]}
            st['same landing boxes'] += lb1 == lb2
            st['dest nests as intended'] += all(Y == (1 if k == 'R' else 2) for k, x, Y, ij in l1 + l2)
            if not same:
                cellok = False
                if len(first) < 2: first.append((key, i, j, R[i], D[j], l1[-2:], l2[-2:]))
    st['cells'] += 1; st['cells all diamonds commute'] += cellok
print(dict(st)); [print(f) for f in first]
