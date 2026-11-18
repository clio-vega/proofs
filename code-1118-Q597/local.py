"""Local configurations of the clean journey's middle (one-star and all |U|).
K = (mid_M, mid_Mx) translated so that the first M-middle rhombus has centroid 0.  For each record:
 (P) do prefix/suffix steps of the clean journey touch I?   (I = observed interaction set, both journeys)
 (V) X on I equals (value first read in Mx-middle) - (value first read in M-middle)?  [needs: cells of I all read by both middles]
 (E) on I: X + Y^x = s  (q: s=(0,-1) for y clean;  p: clean journey is x, the statement is Y + X^y = s)
 Count distinct K (up to translation) per (side,|U|), and whether E is a function of K (it must hold for every K)."""
from collections import Counter, defaultdict
from q597lib import *
SV = {0: (0, -1), 1: (0, 1)}
def firstread(mid):
    out = {}
    for s in mid:
        for k, v in pure(s[1]).items(): out.setdefault(k, v)
    return out
def keys(steps):
    out = set()
    for s in steps: out |= set(pure(s[1]))
    return out
st = Counter(); Ks = defaultdict(set)
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side, nU = r['Us'][0]
    if side == 0: A, B, D_other, D1, D2 = r['Jy'], r['Jyx'], r['X'], r['Y'], r['Yx']
    else:         A, B, D_other, D1, D2 = r['Jx'], r['Jxy'], r['Y'], r['X'], r['Xy']
    mA, mB, p, q = middles(A, B)
    I = set(addd(D2, neg(D1)))
    st[(side, nU, 'P prefix/suffix touch I', bool(keys(A[:p] + A[len(A)-q:]) & I))] += 1
    fA, fB = firstread(mA), firstread(mB)
    covered = I <= set(fA) and I <= set(fB)
    st[(side, nU, 'I read by both middles', covered)] += 1
    if covered:
        st[(side, nU, 'V other-delta on I == firstread diff', all(tuple(b-a for a, b in zip(fA[k], fB[k])) == D_other.get(k, (0, 0)) for k in I))] += 1
    st[(side, nU, 'E on I', all(tuple(a+b for a, b in zip(D_other.get(k, (0, 0)), D2.get(k, (0, 0)))) == SV[side] for k in I))] += 1
    o = tuple(-x for x in cen(rh(mA[0][1])))
    K = (tuple(tr(s, o) for s in mA), tuple(tr(s, o) for s in mB))
    Ks[(side, nU)].add(K)
for k in sorted(st, key=str): print(st[k], k)
for k in sorted(Ks): print('distinct local configurations K', k, len(Ks[k]))
import pickle; pickle.dump(dict(Ks), open('Ks.pkl', 'wb'))
