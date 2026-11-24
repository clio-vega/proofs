"""Global jump rule.  For every M0 (all BC pairs), both k in {0,1}, every edge E_i(x) of Delta (both adjacent triangles
having a unique Lambda_k-lift) with no rhombus over it: jump h(t_{i,i+1}) - h(t_{i-1,i}) = K * n with K = #squares.
Collect n per (k, i).  If n_{i+3} = -n_i and n_0 + n_2 = n_1, then closure (Cl) is a corollary of (H1)+(H2).
Also: do the squares over E have Lambda_{1-k}-edges forming the unit chain from h(t_{i-1,i}) (H4 globally)?
usage: python3 jumpglob.py SIZE"""
import sys
from collections import Counter, defaultdict
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/bandrule.py').read().split('st = Counter(); shapes')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
st = Counter(); NRM = defaultdict(Counter); seen = set()
for P in PAIRS:
    key = frozenset(P['M0'])
    if key in seen: continue
    seen.add(key); M0 = set(P['M0'])
    for k in (0, 1):
        lifts = defaultdict(list)
        for u in M0:
            if zisL(u, k): lifts[zpk(u, k)].append(zfix(min(u[1]), k))
        SQ = defaultdict(list); RH = Counter()
        for u in M0:
            if u[0] == 'sq': SQ[zpk(u, k)].append(frozenset(zfix(Q, k) for Q in u[1]))
            if u[0] == 'rh': RH[zpk(u, k)] += 1
        st[('k', k, 'H1 unique lifts', all(len(v) == 1 for v in lifts.values()))] += 1
        verts = {x for tau in lifts for x in tau}
        for x in verts:
            for i, d in enumerate(DIRS):
                if i >= 3: continue  # each edge once, oriented by i in {0,1,2}
                E = frozenset((x, (x[0] + d[0], x[1] + d[1])))
                a, b = tri_at(x, i, -1), tri_at(x, i, 1)
                if a not in lifts or b not in lifts: continue
                if RH[E]: st[('k', k, 'edge with rhombus')] += 1; continue
                h1, h2 = lifts[a][0], lifts[b][0]; J = (h2[0] - h1[0], h2[1] - h1[1]); K = len(SQ[E])
                if K == 0:
                    st[('k', k, 'K=0 jump zero', J == (0, 0))] += 1; continue
                if J[0] % K or J[1] % K: st[('k', k, 'jump not multiple')] += 1; continue
                n = (J[0] // K, J[1] // K); NRM[(k, i)][n] += 1
                chain = {frozenset(((h1[0] + r * n[0], h1[1] + r * n[1]), (h1[0] + (r + 1) * n[0], h1[1] + (r + 1) * n[1]))) for r in range(K)}
                st[('k', k, 'H4 chain', set(SQ[E]) == chain and len(SQ[E]) == K)] += 1
print('distinct M0:', len(seen))
for kk in sorted(st, key=str): print(st[kk], kk)
for kk in sorted(NRM): print('n', kk, dict(NRM[kk]))
