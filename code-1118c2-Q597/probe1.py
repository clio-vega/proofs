"""Where does D = M_P triangle Mx_P meet the middle?  For each clean record (<=5):
 Mx_P = M + (other move's journey) + prefix;  M_P = M + prefix.  Report |D|, D cap gone(a), D cap windows,
 and whether the other journey J_o contains a translate of a or of steps of C."""
from replay import *
from collections import Counter
st = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side, A, B, Jo = clean(r)
    mA, mB, p, q = middles(A, B)
    M0 = board0(r); Mo = apply(M0, Jo)
    MP = apply(M0, A[:p]); MxP = apply(Mo, B[:p])
    D = MP ^ MxP
    a, C = mA[0], mA[1:]
    st[(side, r['Us'][0][1], 'gone(a) in Mx_P', a[1] <= MxP)] += 1
    st[(side, 'gone(a) cap D', len(a[1] & D))] += 1
    st[(side, 'first Mx step gone cap D', len(mB[0][1] & D))] += 1
    # translates of a / C inside Jo
    sh = lambda s: tr(s, tuple(-v for v in cen(rh(s[1]))))
    shJo = [sh(s) for s in Jo]
    st[(side, 'a-shape in Jo', sh(a) in shJo)] += 1
    st[(side, 'all C-shapes in Jo', all(sh(s) in shJo for s in C))] += 1
    # Jo steps touching the middle region
    reg = set().union(*[s[1] | s[2] for s in mA + mB])
    st[(side, '#Jo steps meeting middle tiles', sum(1 for s in Jo if (s[1] | s[2]) & reg))] += 1
for k in sorted(st, key=str): print(st[k], k)
