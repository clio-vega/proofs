import pickle
from collections import Counter
exec(open('stars.py').read().split('st = Counter()')[0])
st = Counter(); shown = 0
for r in BCo:
    X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
    K = set(X) | set(Y) | set(Yx) | set(Xy)
    I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
    if not I: continue
    side = next(iter(I))[0]; tris = {k[1] for k in I}
    U = sorted(u for T in tris for u in T if star(u) <= tris); U = sorted(set(U))
    for u in U:
        vals = Counter(r['P'][(side, T)] for T in star(u))
        st[('per-star value multiset', tuple(sorted(vals.values())))] += 1
    vals = {r['P'][k] for k in I}
    if len(vals) > 1 and len(U) == 1 and shown < 3:
        shown += 1; u = U[0]; print('--- side', side, 'u', u)
        for T in sorted(star(u), key=sorted):
            k = (side, T); print('   ', sorted(T), 'val', r['P'][k], 'X', X.get(k, Z), 'Y', Y.get(k, Z), 'Yx', Yx.get(k, Z), 'Xy', Xy.get(k, Z))
for k, c in sorted(st.items(), key=str): print(c, k)
