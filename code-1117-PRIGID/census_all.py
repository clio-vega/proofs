"""Step 1 falsifier.  For each (a,b): enumerate edge-to-edge tilings, keep mixed subdivisions, cross-check by Cayley,
then test whether p alone (and q alone) determines the subdivision.  Pre-registered: any collision kills P-RIGID."""
import sys, importlib
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E
from collections import defaultdict
def pq(tiles, T):
    p, q = {}, {}
    for i in T:
        k, vs = tiles[i]; B0, B1 = E.parts(vs)
        if k == 't0': q[B0] = next(iter(B1))
        elif k == 't1': p[B1] = next(iter(B0))
    return frozenset(p.items()), frozenset(q.items())
def run(a, b):
    tiles, poly, R = E.enumerate_all(a, b)
    M = [T for T in R if E.is_mixed(tiles, T)]
    vols = set(); cay_ok = 0
    for T in R:
        V, why = E.cayley_ok(tiles, T, a, b)
        mixed = E.is_mixed(tiles, T)
        if V is not None: vols.add(V)
        if (V is not None) != mixed: print('  DISAGREE mixed=%s cayley=%s' % (mixed, why))
        cay_ok += V is not None
    gp, gq, gn = defaultdict(set), defaultdict(set), defaultdict(set)
    degen = 0
    A0 = set(E.tri_pts(a))
    for T in M:
        p, q = pq(tiles, T); gp[p].add(T); gq[q].add(T)
        if set(v for _, v in p) != A0: degen += 1
    cp = sum(len(s) > 1 for s in gp.values()); cq = sum(len(s) > 1 for s in gq.values())
    print(f'(a,b)=({a},{b}) SIGN={E.SIGN} tilings={len(R)} mixed={len(M)} cayley_ok={cay_ok} volumes={len(vols)} '
          f'with_degenerate_vertex={degen} p-classes={len(gp)} p-collisions={cp} q-classes={len(gq)} q-collisions={cq}')
    return tiles, M, gp, gq
if __name__ == '__main__':
    E.SIGN = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    run(int(sys.argv[1]), int(sys.argv[2]))
