"""Refine the falsifier: split by whether p is surjective onto Vert(P0) (no degenerate vertex)."""
import sys, time, pickle
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E, census_all as C
from collections import defaultdict, Counter
a, b = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time()
tiles, poly, R = E.enumerate_all(a, b)
M = [T for T in R if E.is_mixed(tiles, T)]
A0, A1 = set(E.tri_pts(a)), set(E.tri_pts(b))
gp, gq = defaultdict(set), defaultdict(set)
for T in M:
    p, q = C.pq(tiles, T); gp[p].add(T); gq[q].add(T)
def report(name, G, target):
    surj = {k: s for k, s in G.items() if set(v for _, v in k) == target}
    sizes = Counter(len(s) for s in G.values()); ssz = Counter(len(s) for s in surj.values())
    print(f'  {name}: classes={len(G)} size-hist={dict(sorted(sizes.items()))} | surjective classes={len(surj)} '
          f'tilings={sum(len(s) for s in surj.values())} size-hist={dict(sorted(ssz.items()))}')
print(f'(a,b)=({a},{b}) mixed={len(M)} (edge-to-edge tilings {len(R)}) time={time.time()-t0:.0f}s')
report('p', gp, A0); report('q', gq, A1)
pickle.dump((a, b, [frozenset(tiles[i] for i in T) for T in M]), open(f'all_{a}_{b}.pkl', 'wb'))
