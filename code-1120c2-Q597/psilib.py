"""PROVE 11-20 c2: shared helpers (psi, rh, tr, edgevecs, ETA), copied from future.py."""
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
def rh(S): return next(t for t in S if t[0] == 'rh')
def psi(t, TT):
    eta = ETA[TT]
    return (min(-eta(P[2:]) for P in t[1]), min(eta(P[:2]) for P in t[1]))
def tr(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
def edgevecs(R):
    V = list(R[1]); out = set()
    for P in V:
        for Q in V:
            d = tuple(b-a for a, b in zip(P, Q))
            nz0, nz1 = any(d[:2]), any(d[2:])
            if (nz0 != nz1) and max(abs(x) for x in d) <= 1: out.add(d)   # unit vector in one lattice
    return out
