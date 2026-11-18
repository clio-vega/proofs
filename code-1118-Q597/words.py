"""Middle words: typ word + rhombus-centroid displacement word of y-middle (M vs Mx) and x-middle (M vs My).
Tests: (W1) same multiset of displacements; (W2) rhombus after the middle identical; (W3) word shapes."""
import pickle
from fractions import Fraction as Fr
from collections import Counter
exec(open('analyse_align.py').read().split('T = Counter()')[0])
def rh(tset):
    for k, vs in tset:
        if k == 'rh': return vs
def cen(vs): return tuple(sum(Fr(v[i]) for v in vs)/len(vs) for i in range(4))
def disp(step):
    _, g, n = step; a, b = cen(rh(g)), cen(rh(n)); return tuple(y-x for x, y in zip(a, b))
def name(d):
    return {(0, -1, 0, 0): 'A0-', (0, 0, 0, -1): 'A1-'}.get(tuple(int(x) if x.denominator == 1 else x for x in d), str(tuple(float(x) for x in d)))
st = Counter(); shapes = Counter(); names = {}
for r in R:
    c = cls(r)
    if c in ('BC disjoint', 'BB com'): continue
    for who, A, B in (('y', r['Jy'], r['Jyx']), ('x', r['Jx'], r['Jxy'])):
        p = pre(A, B); q = suf(A, B, p); mA, mB = A[p:len(A)-q], B[p:len(B)-q]
        dA, dB = [disp(s) for s in mA], [disp(s) for s in mB]
        st[(c, who, 'W1 same multiset', Counter(dA) == Counter(dB))] += 1
        st[(c, who, 'W2 same rhombus after middle', rh(mA[-1][2]) == rh(mB[-1][2]))] += 1
        for d in dA + dB: names.setdefault(d, 'v%d' % len(names))
        if c.startswith('BC'):
            shapes[(c, who, ''.join(names[d] for d in dA), ''.join(names[d] for d in dB))] += 1
for k in sorted(st, key=str): print(st[k], k)
print({v: tuple(str(x) for x in k) for k, v in names.items()})
for k in sorted(shapes, key=str): print(shapes[k], k)
