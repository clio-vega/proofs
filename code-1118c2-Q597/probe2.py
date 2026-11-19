"""Half-wake translation: D_rem = M_P \ Mx_P (x's wake before), D_add = Mx_P \ M_P (after).
How much of D_rem, translated by -alpha, lies in D_add?  Which middle-step tiles (removed by C in M, C-alpha in Mx) are in D / common?"""
from replay import *
from collections import Counter
def T(S, v): return {(k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S}
st = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side, A, B, Jo = clean(r)
    mA, mB, p, q = middles(A, B)
    M0 = board0(r); MP = apply(M0, A[:p]); MxP = apply(apply(M0, Jo), B[:p])
    rem, add = MP - MxP, MxP - MP
    a, C = mA[0], mA[1:]
    alpha = tuple(y-x for x, y in zip(cen(rh(a[1])), cen(rh(a[2]))))
    ma = tuple(-v for v in alpha)
    tr_ = T(rem, ma) & add
    st[(side, 'frac of wake(rem) translated by -alpha into add', round(len(tr_)/len(rem), 2) if rem else None)] += 0
    # tiles of the board M_1 = MP.a that C reads (removes), and their -alpha translates in MxP
    B1 = apply(MP, [a]); need = set()
    Bc = set(B1)
    for s in C:
        need |= (s[1] & B1)   # original tiles of M_1 consumed by C
        Bc = apply(Bc, [s])
    ok = T(need, ma) <= MxP
    st[(side, 'C-consumed tiles of M_1, shifted -alpha, present in Mx_P', ok)] += 1
    kinds = Counter('D' if t in rem else ('a-new' if t in a[2] else 'common') for t in need)
    st[(side, 'origin of C-consumed tiles', tuple(sorted(kinds.items())))] += 1
    # rhombus edge vectors of x
    xr = rh(Jo[0][1]); vs = sorted(xr)
    # a-consumed tiles of MP: translated +gamma present in Mx_P.(C-alpha)?
    gamma = tuple(y-x for x, y in zip(cen(rh(C[0][1])), cen(rh(C[-1][2]))))
    Bx = apply(MxP, mB[:-1])
    st[(side, 'a-consumed tiles shifted +gamma present after C-alpha', T(a[1], gamma) <= Bx)] += 1
    kinds = Counter('D' if t in rem else 'common' for t in a[1])
    st[(side, 'origin of a-consumed', tuple(sorted(kinds.items())))] += 1
for k in sorted(st, key=str):
    if st[k]: print(st[k], k)
