"""Is the Mx-middle a cyclic rotation of the M-middle, step by step a TRANSLATE (gone and new tile sets moved by one Z^4 vector)?
For each record report: k (rotation), translation-vector word, or 'none'."""
import pickle
from fractions import Fraction as Fr
from collections import Counter
exec(open('words.py').read().split('st = Counter()')[0])
def tr(step, v):
    typ, g, n = step
    f = lambda S: frozenset((k, frozenset(tuple(a+b for a, b in zip(P, v)) for P in vs)) for k, vs in S)
    return (typ, f(g), f(n))
def tvec(s1, s2):
    """vector v with s2 = s1 + v, or None"""
    if s1[0] != s2[0]: return None
    v = tuple(b-a for a, b in zip(cen(rh(s1[1])), cen(rh(s2[1]))))
    return v if tr(s1, v) == s2 else None
st = Counter(); ex = {}
for r in R:
    c = cls(r)
    if c in ('BC disjoint', 'BB com'): continue
    for who, A, B in (('y', r['Jy'], r['Jyx']), ('x', r['Jx'], r['Jxy'])):
        p = pre(A, B); q = suf(A, B, p); mA, mB = A[p:len(A)-q], B[p:len(B)-q]
        found = None
        if len(mA) == len(mB):
            n = len(mA)
            for k in range(n):
                vs = [tvec(mA[(i+k) % n], mB[i]) for i in range(n)]
                if all(v is not None for v in vs): found = (k, n); break
        tag = 'none' if found is None else ('rot k=%d of n=%d' % found if found[0] not in (1, found[1]-1) else ('front->back' if found[0] == 1 else 'back->front'))
        st[(c, who, tag)] += 1
        ex.setdefault((c, who, tag), r)
for k in sorted(st, key=str): print(st[k], k)
pickle.dump(ex, open('transl-ex.pkl', 'wb'))
