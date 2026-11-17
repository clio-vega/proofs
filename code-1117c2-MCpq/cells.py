"""Per-cell pattern (X(k), Y(k), Yx(k), Xy(k)) over all cells touched by any of the four moves.  0 = no move."""
import pickle, sys
from collections import Counter
BCo, BBn = pickle.load(open('sets.pkl', 'rb'))
def vec(D): return {k: (n[0]-o[0], n[1]-o[1]) for k, (o, n) in D.items()}
Z = (0, 0)
for name, S in (('BC', BCo), ('BB', BBn)):
    pat = Counter(); per = Counter()
    for r in S:
        X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
        K = set(X) | set(Y) | set(Yx) | set(Xy); n = 0
        for k in K:
            t = (k[0], X.get(k, Z), Y.get(k, Z), Yx.get(k, Z), Xy.get(k, Z))
            if not (t[3] == t[2] and t[4] == t[1]): pat[t] += 1; n += 1
        per[n] += 1
    print(name, 'interacting cells per pair', sorted(per.items()))
    print('  side  X       Y       Yx      Xy')
    for k, v in sorted(pat.items(), key=lambda kv: -kv[1]): print('  ', v, k)
