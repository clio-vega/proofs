"""Exchange test.  M = census mosaic, alpha = A-nest (cells (col,row)).  R = (ii) order (top row first, E->W),
D = (iii) order (depth 1 first, E->W within depth).  (iii) on M: send D in order to C (dir u_A).  M_1 = M after sending the
top row R[:alpha_l] to B (dir w_A).  (iii) on M_1: send D(alpha-hat) to C.
Exchange: first m-alpha_l landings of (iii) on M == landings of (iii) on M_1.
Also record: the cell moved at each step on M vs on M_1, and whether they are related by the one-row translation
(c, r) on M  <->  (c, r-1) on M_1 for c < alpha_l (the 'flock shift' phi), identity otherwise."""
from instr import *
from collections import Counter
def orders(cells):
    conj = Counter(c for c, r in cells)
    R = sorted(cells, key=lambda x: (-x[1], -x[0]))
    D = sorted(cells, key=lambda x: (conj[x[0]] - x[1], -x[0]))
    return R, D
def run(b0, seq):
    b = copy.deepcopy(b0); out = []
    for kind, x, dv in seq:
        Y, ij, s = mig1(b, 0, x, dv); out.append((kind, x, Y, ij))
    return b, out
st = Counter(); first = []
if __name__ != "__main__": census = lambda: iter(())
for key, b0 in census():
    A = b0.nest[0]; uA, wA = A['u'], A['w']; cells = set(A['pos']); R, D = orders(cells); m = len(R)
    ell = max(r for c, r in cells) + 1; top = [x for x in R if x[1] == ell-1]; al = len(top)
    _, iiiM = run(b0, [('C', x, uA) for x in D])
    bM1, iiM1 = run(b0, [('B', x, wA) for x in top])
    hat = cells - set(top); Rh, Dh = orders(hat)
    _, iiiM1 = run(bM1, [('C', x, uA) for x in Dh])
    st['cells'] += 1
    lM = [ij for _, _, _, ij in iiiM][:m-al]; lM1 = [ij for _, _, _, ij in iiiM1]
    ok = lM == lM1; st['Exchange holds'] += ok
    st['Exchange holds as SETS'] += set(lM) == set(lM1)
    # where does the cell moved differ, and how
    for k in range(m-al):
        xM, xM1 = iiiM[k][1], iiiM1[k][1]
        if xM == xM1: st['step: same cell'] += 1
        elif xM == (xM1[0], xM1[1]+1): st['step: M cell = M1 cell shifted up one row'] += 1
        else: st['step: other relation'] += 1
        st['step: same landing'] += iiiM[k][3] == iiiM1[k][3]
    if not ok and len(first) < 3: first.append((key, D, Dh, lM, lM1))
    if len(Dh) == 0: st['alpha one row (vacuous)'] += 1
for k in sorted(st): print(k, st[k])
for f in first: print(f)
