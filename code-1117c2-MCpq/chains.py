"""Chain-level determinacy of I from solo data.
Pass 1 (learn): for each BC record, U per side, ordered along d_side ((1,-1) on q, (1,0) on p); role of each u: single/first/inner/last;
   collect signature sets per (side, role).
Pass 2 (predict): from solo data only, U_pred = union of maximal runs u, u+d, ..., u+jd with sig(u) in first (or single if j=0),
   sig(inner) for middle, sig(last) for end.  Compare I_pred = union of stars with true I on all BC records (disjoint ones too),
   and report how BB records look under the same rule.  NOTE pass 2 uses signatures learned in pass 1 on the same data: this tests
   consistency (a function exists), not out-of-sample prediction.  Out-of-sample: learn on even-indexed records, test on odd."""
import pickle
from collections import Counter, defaultdict
exec(open('determ.py').read().split("for kind in ('BC', 'BB'):")[0])
D = {0: (1, -1), 1: (1, 0)}
def solo(r):
    X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
    K = set(X) | set(Y) | set(Yx) | set(Xy)
    I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
    return X, Y, I
def sigs(r, X, Y, side):
    tris = {k[1] for k in r['P'] if k[0] == side}; out = {}
    for u in {u for T in tris for u in T}:
        P6 = pos(u)
        if all(T in tris for T in P6): out[u] = tuple((X.get((side, T), Z), Y.get((side, T), Z)) for T in P6)
    return out
def add(u, d, j=1): return (u[0]+j*d[0], u[1]+j*d[1])
def learn(R):
    S = defaultdict(set)
    for r in R:
        X, Y, I = solo(r)
        for side in (0, 1):
            sg = sigs(r, X, Y, side)
            U = sorted(u for u in sg if all((side, T) in I for T in pos(u)))
            Us = set(U)
            for u in U:
                pre, nxt = add(u, D[side], -1) in Us, add(u, D[side]) in Us
                role = {(False, False): 'single', (False, True): 'first', (True, True): 'inner', (True, False): 'last'}[(pre, nxt)]
                S[(side, role)].add(sg[u])
    return S
def predict(r, S):
    X, Y, I = solo(r); Ipred = set()
    for side in (0, 1):
        sg = sigs(r, X, Y, side); d = D[side]
        for u, s in sg.items():
            if s in S[(side, 'single')]: Ipred |= {(side, T) for T in pos(u)}
            if s in S[(side, 'first')]:
                j = 1
                while sg.get(add(u, d, j)) in S[(side, 'inner')]: j += 1
                if sg.get(add(u, d, j)) in S[(side, 'last')]:
                    for i in range(j+1): Ipred |= {(side, T) for T in pos(add(u, d, i))}
    return Ipred, I
BC = [r for r in recs if r['kind'] == 'BC']; BB = [r for r in recs if r['kind'] == 'BB']
S = learn(BC)
for k in sorted(S): print('role', k, len(S[k]))
roles = list(S); print('role signature sets pairwise disjoint (within side):',
      all(not (S[a] & S[b]) for a in roles for b in roles if a < b and a[0] == b[0]))
res = Counter()
for r in BC:
    Ip, I = predict(r, S); res[('BC', Ip == I, 'I empty', not I)] += 1
for r in BB:
    Ip, I = predict(r, S); res[('BB', 'pred nonempty', bool(Ip))] += 1
print(dict(res))
S2 = learn(BC[0::2]); oos = Counter()
for r in BC[1::2]:
    Ip, I = predict(r, S2); oos[(Ip == I, 'I empty', not I)] += 1
print('out-of-sample (learn even, test odd):', dict(oos))
pickle.dump(dict(S), open('rolesigs.pkl', 'wb'))
