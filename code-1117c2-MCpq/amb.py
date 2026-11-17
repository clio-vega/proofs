import pickle
exec(open('determ.py').read().split("for kind in ('BC', 'BB'):")[0])
target = (((0, 0), (0, 0)), ((0, -1), (0, -1)), ((1, -1), (0, -1)), ((0, 0), (0, 0)), ((1, -1), (0, 0)), ((0, -1), (0, -1)))
for r in recs:
    if r['kind'] != 'BC': continue
    X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
    K = set(X) | set(Y) | set(Yx) | set(Xy)
    I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
    tris = {k[1] for k in r['P'] if k[0] == 0}
    for u in {u for T in tris for u in T}:
        P6 = pos(u)
        if not all(T in tris for T in P6): continue
        sig = tuple((X.get((0, T), Z), Y.get((0, T), Z)) for T in P6)
        if sig == target and not all((0, T) in I for T in P6):
            print('record', r['key'], r['x'], r['y'], 'u', u, '|I|', len(I), 'overlap', bool(set(X) & set(Y)))
            for k in sorted(K, key=lambda k: (k[0], sorted(k[1]))):
                print('  ', k[0], sorted(k[1]), 'val', r['P'][k], 'X', X.get(k, Z), 'Y', Y.get(k, Z), 'Yx', Yx.get(k, Z), 'Xy', Xy.get(k, Z), 'I' if k in I else '')
