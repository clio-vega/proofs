import pickle
from collections import Counter
BCo, BBn = pickle.load(open('sets.pkl', 'rb'))
def vec(D): return {k: (n[0]-o[0], n[1]-o[1]) for k, (o, n) in D.items()}
Z = (0, 0)
S_DIR = {0: (0, -1), 1: (0, 1)}
stat = Counter(); shown = 0
for r in BCo:
    X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
    K = set(X) | set(Y) | set(Yx) | set(Xy)
    I = [k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))]
    if not I: continue
    Ssides = Counter(); Tside = Counter()
    for k in I:
        pt = (X.get(k, Z), Y.get(k, Z))
        s = S_DIR[k[0]]
        if pt in ((s, s), (Z, Z)): Ssides[(k[0], pt == (s, s))] += 1
        else: Tside[k[0]] += 1
    stat[('T', tuple(sorted(Tside.items())), 'S', tuple(sorted(Ssides.items())))] += 1
    # composite total on all cells vs X+Y
    D = {k: tuple(a+b for a, b in zip(X.get(k, Z), Yx.get(k, Z))) for k in K}
    stat[('all interacting cells composite = s', all(D[k] == S_DIR[k[0]] for k in I))] += 1
    stat[('noninteracting both-moved cells', sum(1 for k in K if k not in I and k in X and k in Y))] += 1
    if shown < 2 and len(I) == 6:
        shown += 1; print('--- record', r['key'], r['x'], r['y'])
        for k in sorted(I, key=str):
            print('  side', k[0], sorted(k[1]), 'val', r['P'][k], 'X', X.get(k, Z), 'Y', Y.get(k, Z), 'Yx', Yx.get(k, Z), 'Xy', Xy.get(k, Z))
for k, v in sorted(stat.items(), key=str): print(v, k)
