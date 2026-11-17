"""Structural probes on ALL subdivisions (all_a_b.pkl):
 (A) is every vertex of S a vertex of a pure cell or of the boundary?
 (B) is S determined by its vertex set V(S)?
 (C) is the wall gamma_e a shortest path in T1 between its endpoints (and is that shortest path unique)?"""
import sys, pickle, itertools
from collections import defaultdict, Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E
def dist(b):
    P = E.tri_pts(b); adj = defaultdict(set)
    for e in E.unit_edges(b):
        u, v = tuple(e); adj[u].add(v); adj[v].add(u)
    D = {}
    for s in P:
        d = {s: 0}; fr = [s]
        while fr:
            nf = []
            for u in fr:
                for v in adj[u]:
                    if v not in d: d[v] = d[u]+1; nf.append(v)
            fr = nf
        D[s] = d
    return D
for fn in sys.argv[1:]:
    a, b, M = pickle.load(open(fn, 'rb'))
    st = Counter(); gV = defaultdict(set); D1 = dist(b)
    for T in M:
        V = set(v for _, vs in T for v in vs)
        Vp = set(v for k, vs in T if k != 'par' for v in vs)
        st['vertices not on pure cell'] += len(V - Vp); st['vertices'] += len(V)
        gV[frozenset(V)].add(T)
        # walls gamma_e: for each 0-edge e, the set of f with e x f in S
        ge = defaultdict(list)
        for k, vs in T:
            if k == 'par':
                B0, B1 = E.parts(vs); ge[B0].append(B1)
        q = {}
        for k, vs in T:
            if k == 't0': B0, B1 = E.parts(vs); q[B0] = next(iter(B1))
        for e, fs in ge.items():
            # endpoints: vertices of odd degree in the multigraph of fs
            deg = Counter(x for f in fs for x in f); ends = [x for x, c in deg.items() if c % 2]
            if len(ends) == 2:
                st['interior-ish walls'] += 1
                if len(fs) == D1[ends[0]][ends[1]]: st['wall is geodesic'] += 1
    print(fn, dict(st), '| V(S)-classes', len(gV), 'of', len(M), 'collisions', sum(len(s) > 1 for s in gV.values()))
