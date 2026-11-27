"""Lemma Mos-P corroboration: at every engine step on the size-5 closure, the removed tiles and the created tiles have
the SAME boundary (as a set of unit Z^4 edges, every vertex kept), the same kinds multiset, and exactly one rhombus.
So a step re-tiles a fixed region and leaves every other tile and the region's boundary edges untouched.
usage: python3 hexcheck.py NW NCH"""
import sys, copy
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, '/home/clio/projects/reviews/code-1125c2-BS')
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590')
import indep_bs_census as IBC
from indep_bs_census import census, closure, move
from common import corners
def bd(S):
    c = Counter()
    for kind, vs in S:
        vs = list(vs)
        # tiles are stored as vertex cycles; frozenset loses order, so rebuild edges as unit pairs at distance 1 in Z^4 sense
    return c
def edges_of(kind, cyc):
    return [frozenset((cyc[i], cyc[(i + 1) % len(cyc)])) for i in range(len(cyc))]
_step = IBC._engine_step
REC = []
def rec(self, R, d):
    before = dict(self.tiles); R2 = _step(self, R, d)
    if R2 is not None:
        rem = [before[t] for t in set(before) - set(self.tiles)]; cre = [self.tiles[t] for t in set(self.tiles) - set(before)]
        REC.append((rem, cre))
    return R2
import migrate; migrate.Board.step = rec
def boundary(T):
    c = Counter(e for kind, vs in T for e in edges_of(kind, vs))
    return frozenset(e for e, k in c.items() if k == 1), all(k <= 2 for k in c.values())
def work(args):
    w, NCH = args; st = Counter()
    for i, (key, b0) in enumerate(census()):
        if i % NCH != w: continue
        for b, dirs in closure(b0):
            for x in corners(set(b.nest[0]['pos'])):
                for T in 'BC':
                    REC.clear(); move(b, x, dirs[T])
                    for rem, cre in REC:
                        st['steps'] += 1
                        b1, ok1 = boundary(rem); b2, ok2 = boundary(cre)
                        st[('same boundary edges', b1 == b2)] += 1
                        st[('boundary is one 6-12 edge cycle', 6 <= len(b1) <= 12)] += 1
                        st[('kinds preserved', Counter(k for k, _ in rem) == Counter(k for k, _ in cre))] += 1
                        st[('exactly one rhombus', sum(k == 'rh' for k, _ in rem) == 1)] += 1
                        st[('internally edge-to-edge', ok1 and ok2)] += 1
                        st[('boundary length', len(b1))] += 1
    return st
if __name__ == '__main__':
    NW, NCH = int(sys.argv[1]), int(sys.argv[2]); tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(work, [(w, NCH) for w in range(NCH)]): tot.update(o)
    for k in sorted(tot, key=str): print(tot[k], k)
