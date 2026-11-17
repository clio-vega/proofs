import pickle
from collections import Counter
BCo, BBn = pickle.load(open('sets.pkl', 'rb'))
exec(open('analyse1.py').read().split('BC = [')[0].split("recs = pickle")[0])
def vec(D): return {k: (n[0]-o[0], n[1]-o[1]) for k, (o, n) in D.items()}
def sub(A, B):
    out = {}
    for k in set(A) | set(B):
        a = A.get(k, (0, 0)); b = B.get(k, (0, 0)); d = (a[0]-b[0], a[1]-b[1])
        if d != (0, 0): out[k] = d
    return out
st = Counter()
for name, S in (('BC', BCo), ('BB', BBn)):
  for r in S:
    X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
    for k, v in X.items(): st[(name, 'X vec', k[0], v)] += 1
    for k, v in Y.items(): st[(name, 'Y vec', k[0], v)] += 1
    cY = sub(Yx, Y)
    for k, v in cY.items(): st[(name, 'cY vec', k[0], v, 'inX', k in X, 'inY', k in Y)] += 1
for k, v in sorted(st.items(), key=str): print(v, k)
