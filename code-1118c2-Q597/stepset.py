"""Board.step re-implemented on a bare tile set (kind, frozenset(vertices)), mirroring migrate.Board.step line by line,
but RETURNING the candidate list and tie-break data, so candidate-level covariance can be tested."""
import sys; sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import migrate
from migrate import G, prot, norm, Fr
from collections import Counter, defaultdict
from itertools import combinations
UNIT = {(1, 0), (0, 1), (-1, 0), (0, -1), (1, -1), (-1, 1)}
def is_edge(P, Q):
    d = tuple(b-a for a, b in zip(P, Q))
    return (d[2:] == (0, 0) and d[:2] in UNIT) or (d[:2] == (0, 0) and d[2:] in UNIT)
def tedges(t): return [frozenset(e) for e in combinations(t[1], 2) if is_edge(*e)]
def index(Bd):
    E = defaultdict(set)
    for t in Bd:
        for e in tedges(t): E[e].add(t)
    return E
def nbrs(t, E):
    out = set()
    for e in tedges(t): out |= E[e]
    out.discard(t); return out
def bcycle(S):
    cnt = Counter()
    for t in S:
        for e in tedges(t): cnt[e] += 1
    bd = [e for e, c in cnt.items() if c == 1]; adj = {}
    for e in bd:
        p, q = tuple(e); adj.setdefault(p, []).append(q); adj.setdefault(q, []).append(p)
    if any(len(v) != 2 for v in adj.values()): return None
    start = next(iter(adj)); cyc = [start]; prev = None; cur = start
    while True:
        nxt = [q for q in adj[cur] if q != prev][0] if prev is not None else adj[cur][0]
        if nxt == start: break
        cyc.append(nxt); prev, cur = cur, nxt
        if len(cyc) > 20: return None
    return cyc if len(cyc) == len(bd) else None
def candidates(Bd, R, d, E=None):
    E = E or index(Bd); rot = d.conjugate()
    x = lambda P: (G.phi(P)*rot).real
    xmin = min(x(P) for P in R[1]); out = []
    frontier = {frozenset([R])}
    for size in (2, 3, 4):
        new = set()
        for S in frontier:
            for t in S:
                for nb in nbrs(t, E):
                    if nb in S or nb[0] == 'rh': continue
                    new.add(S | {nb})
        frontier = new
    for S in frontier:
        if Counter(t[0] for t in S) != Counter(rh=1, tri=2, sq=1): continue
        cyc = bcycle(S)
        if cyc is None: continue
        zz = [G.phi(P) for P in cyc]; m_ = len(zz)
        cyc = [cyc[i] for i in range(m_) if abs(((zz[i]-zz[i-1]).conjugate()*(zz[(i+1) % m_]-zz[i])).imag) > 1e-9]
        if len(cyc) != 6: continue
        if min(x(P) for t in S for P in t[1]) < xmin - 1e-9: continue
        c = tuple(sum(Fr(P[k]) for P in cyc) / 6 for k in range(4)); V = set(tuple(Fr(p) for p in P) for P in cyc)
        rc = lambda P, k: tuple(a+b for a, b in zip(c, prot(tuple(Fr(p)-q for p, q in zip(P, c)), k)))
        if {rc(P, 3) for P in V} == V: typ = 'a'
        elif {rc(P, 2) for P in V} == V: typ = 'b'
        else: continue
        out.append((typ, frozenset(S), c))
    return out
def across(Bd, R, d, E=None):
    """tie-break tile: across the horizontal edge of R at its (unique) leftmost vertex; None if not applicable"""
    E = E or index(Bd); rot = d.conjugate(); vs = sorted(R[1])
    # need cyclic order of R's vertices: order by angle around centroid
    import cmath
    zc = sum(G.phi(P) for P in vs)/4; vs = sorted(vs, key=lambda P: cmath.phase(G.phi(P)-zc))
    z = [G.phi(P)*rot for P in vs]; xs = [w.real for w in z]; lm = min(range(4), key=lambda i: xs[i])
    if sum(1 for w in xs if abs(w-xs[lm]) < 1e-9) != 1: return None
    for nb_i in ((lm+1) % 4, (lm-1) % 4):
        if abs((z[nb_i]-z[lm]).imag) < 1e-9:
            return frozenset(E.get(frozenset((vs[lm], vs[nb_i])), set()) - {R})
    return None
def step(Bd, R, d):
    E = index(Bd); cands = candidates(Bd, R, d, E)
    if len(cands) > 1:
        ac = across(Bd, R, d, E)
        if ac is not None: cands = [z for z in cands if ac & z[1]]
    if len(cands) != 1: return None, cands
    typ, S, c = cands[0]
    rc = lambda P, k: norm(tuple(a+b for a, b in zip(c, prot(tuple(Fr(p)-q for p, q in zip(P, c)), k))))
    if typ == 'a': k = 3
    else:
        rot = d.conjugate(); ks = []
        for k in (2, 4):
            z = [G.phi(rc(P, k))*rot for P in R[1]]
            # adjacency in cyclic order needed: test all pairs that are edges after rotation
            Rv = [rc(P, k) for P in R[1]]
            ok = any(abs((G.phi(P)*rot - G.phi(Q)*rot).imag) < 1e-9 or abs((G.phi(P)*rot - G.phi(Q)*rot).real) < 1e-9
                     for P, Q in combinations(Rv, 2) if is_edge(P, Q))
            if ok: ks.append(k)
        if len(ks) != 1: return None, cands
        k = ks[0]
    new = frozenset((t[0], frozenset(rc(P, k) for P in t[1])) for t in S)
    return (typ, S, new), cands
