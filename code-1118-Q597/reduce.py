"""Reduction check. Clean journey: q-side -> y (second journey in Mx vs M); p-side -> x (in My vs M).
(R1) prefix/suffix steps are literally equal, so  D2 - D1 = sum delta(middle2) - sum delta(middle1)   [additivity]
(R2) under the block swap  mid1 = a.C,  mid2 = (C - alpha).(a + gamma):
     D2 - D1 = sum_{s in C} (T_{-alpha} delta(s) - delta(s)) + T_gamma delta(a) - delta(a),  T_v shifts keys by v_side.
(R3) one-star: support of D2-D1 is star(u) (6 cells) and the per-cell values; record the step-shape catalogue of the middles
     (each step up to translation) to see how many local configurations occur."""
import pickle
from collections import Counter
from q597lib import *
def T(d, v):   # shift delta d by translation v (Z^4): side-0 keys by v[:2], side-1 keys by v[2:]
    out = {}
    for (side, key), val in d.items():
        w = v[:2] if side == 0 else v[2:]
        out[(side, frozenset((a+w[0], b+w[1]) for a, b in key))] = val
    return out
def neg(d): return {k: tuple(-x for x in v) for k, v in d.items()}
def norm_(v): return tuple(int(x) if x.denominator == 1 else x for x in v)
def shape(step):
    typ, g, n = step; o = cen(rh(g)); z = tuple(-x for x in o)
    return tr(step, z)
st = Counter(); shapes = Counter()
for r in R:
    if r['kind'] != 'BC' or not r['Us']: continue
    side = r['Us'][0][0]
    A, B, D1, D2, who = (r['Jy'], r['Jyx'], r['Y'], r['Yx'], 'y') if side == 0 else (r['Jx'], r['Jxy'], r['X'], r['Xy'], 'x')
    p = pre(A, B); q = suf(A, B, p); mA, mB = A[p:len(A)-q], B[p:len(B)-q]
    diff = addd(D2, neg(D1))
    s1 = {}; s2 = {}
    for s in mA: s1 = addd(s1, sdelta(s))
    for s in mB: s2 = addd(s2, sdelta(s))
    st[('R1 additivity', addd(s2, neg(s1)) == diff)] += 1
    a, C = mA[0], mA[1:]
    alpha = tuple(b-x for x, b in zip(cen(rh(a[1])), cen(rh(a[2]))))
    gamma = tuple(b-x for x, b in zip(cen(rh(C[0][1])), cen(rh(C[-1][2]))))
    pred = {}
    for s in C: pred = addd(addd(pred, T(sdelta(s), tuple(-x for x in alpha))), neg(sdelta(s)))
    pred = addd(addd(pred, T(sdelta(a), gamma)), neg(sdelta(a)))
    st[('R2 block-swap formula', pred == diff)] += 1
    st[('|U|', r['Us'][0][1], 'support size', len(diff))] += 1
    st[('alpha', side, norm_(alpha))] += 1
    for s in mA + mB: shapes[shape(s)] += 1
for k in sorted(st, key=str): print(st[k], k)
print('distinct step shapes (up to translation) in clean middles:', len(shapes))
