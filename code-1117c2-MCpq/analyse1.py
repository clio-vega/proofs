"""Step 1 report + step 2 candidate laws.  Vectors: net change d(k) = new - old in Z^2 for each pure cell k=(i,t).
Laws tested (each on BC and on BB-noncommuting control):
 L0  composites agree on (p,q)                         [= MC by PQ-RIGID]
 La  supp(X) & supp(Y) == {}  (outright disjoint)
 Lb  Yx == Y and Xy == X     (no interaction at all)
 Lc  additive exchange: (Yx - Y) == (Xy - X) as vector fields   [MC restated]
 Ld  Yx - Y == -(X restricted to supp Y) ... printed as structure stats
"""
import pickle, sys
from collections import Counter
recs = pickle.load(open('recs.pkl', 'rb'))
def vec(D):
    return {k: (n[0]-o[0], n[1]-o[1]) for k, (o, n) in D.items()}
def sub(A, B):
    out = {}
    for k in set(A) | set(B):
        a = A.get(k, (0, 0)); b = B.get(k, (0, 0)); d = (a[0]-b[0], a[1]-b[1])
        if d != (0, 0): out[k] = d
    return out
BC = [r for r in recs if r['kind'] == 'BC']; BB = [r for r in recs if r['kind'] == 'BB']
ov = lambda r: bool(set(r['X']) & set(r['Y']))
print('BC total', len(BC), 'overlapping', sum(map(ov, BC)), 'commute', sum(r['com'] for r in BC), 'pq-same', sum(r['pqsame'] for r in BC))
print('BB total', len(BB), 'overlapping', sum(map(ov, BB)), 'noncommuting', sum(not r['com'] for r in BB))
print('com == pqsame everywhere (PQ-RIGID sanity):', all(r['com'] == r['pqsame'] for r in recs))
BCo = [r for r in BC if ov(r)]; BBn = [r for r in BB if not r['com']]
pickle.dump((BCo, BBn), open('sets.pkl', 'wb'))
# class by |add|,|rem| on key sets, as in exchange.py
def cls(r):
    C, Cp = set(r['Y']), set(r['Yx']); return (len(Cp - C), len(C - Cp))
print('BC overlapping classes', Counter(cls(r) for r in BCo))
triv = [r for r in BCo if cls(r) == (0, 0)]
print('class (0,0): composite agree trivially (Yx==Y and Xy==X):', sum(r['Yx'] == r['Y'] and r['Xy'] == r['X'] for r in triv), '/', len(triv))
for name, S in (('BC-overlap', BCo), ('BB-noncomm', BBn)):
    st = Counter()
    for r in S:
        X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
        cY, cX = sub(Yx, Y), sub(Xy, X)
        st['Lc exchange cY==cX'] += cY == cX
        st['Lb no interaction'] += (cY == {} and cX == {})
        st['supp cY <= supp X'] += set(cY) <= set(X)
        st['supp cY <= supp X & supp Y'] += set(cY) <= set(X) & set(Y)
        st['types in supp X (q,p)'] += 0
        st[('X types', tuple(sorted(Counter(k[0] for k in X).items())), 'Y types', tuple(sorted(Counter(k[0] for k in Y).items())))] += 1
    print(name, len(S)); [print('   ', k, v) for k, v in sorted(st.items(), key=str)]
