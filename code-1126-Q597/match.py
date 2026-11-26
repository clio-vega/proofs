"""Match E's (and L's) middle steps between journey and re-run as rigid Z^4 translates.
usage: python3 match.py  (reads /tmp/step0_5_*.pkl)"""
import pickle, glob, sys
from collections import Counter
def canon(s):
    """step -> (shape, anchor): shape = removed+created tile sets translated so that min point is 0"""
    pts = [p for t in s[1] for p in t[1]]
    m = min(pts)
    sh = lambda S: frozenset((t[0], frozenset(tuple(a-b for a, b in zip(p, m)) for p in t[1])) for t in S)
    return (s[0], sh(s[1]), sh(s[2])), m
def vsub(a, b): return tuple(x-y for x, y in zip(a, b))
def matchmid(A, Ar, p, la, lb):
    """for each re-run middle step, the index (rel. to p) of the journey middle step it translates, and the vector"""
    O = [canon(A[i]) for i in range(p, p+la)]
    R = [canon(Ar[i]) for i in range(p, p+lb)]
    out = []
    used = set()
    for j, (shr, mr) in enumerate(R):
        cands = [i for i, (sho, mo) in enumerate(O) if sho == shr and i not in used]
        if not cands: out.append((None, None)); continue
        i = cands[0]; used.add(i); out.append((i, vsub(mr, O[i][1])))
    return out
if __name__ == '__main__':
    st = Counter(); ex = {}
    for f in sorted(glob.glob('/tmp/step0_5_*.pkl')):
        for case, L, info in pickle.load(open(f, 'rb')):
            for r, (T, A, Ar, p, la, lb) in enumerate(info):
                role = 'L' if r == L else 'E'
                m = matchmid(A, Ar, p, la, lb)
                perm = tuple(i for i, v in m)
                vecs = tuple(v for i, v in m)
                full = None not in perm
                st[(case, role, 'all rerun mid steps are translates of journey mid steps', full)] += 1
                if full:
                    key = (case, role, T, 'perm', perm, 'vecs', vecs)
                    st[key] += 1
    for k in sorted(st, key=str): print(st[k], k)
