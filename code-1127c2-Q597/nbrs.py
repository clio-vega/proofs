"""Gap-sharpening probe for (NS) (not a proof): at every pre-step state of every census journey (size-5 closure), what
tiles share an edge with the moving rhombus R?  In particular: is R ever edge-adjacent to a nest rhombus?
usage: python3 nbrs.py NW NCH"""
import sys, copy
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, '/home/clio/projects/reviews/code-1125c2-BS')
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590')
import indep_bs_census as IBC
from indep_bs_census import census, closure
from common import corners
from flock import unit_dir
import migrate
_step = IBC._engine_step
ST = Counter()
def rec(self, R, d):
    ks = Counter(self.tiles[t][0] for t in self.neighbours(R))
    ST[('neighbour kinds of R', tuple(sorted(ks.items())))] += 1
    return _step(self, R, d)
migrate.Board.step = rec
def work(args):
    w, NCH = args; ST.clear()
    for i, (key, b0) in enumerate(census()):
        if i % NCH != w: continue
        for b, dirs in closure(b0):
            for x in corners(set(b.nest[0]['pos'])):
                for T in 'BC':
                    nb = copy.deepcopy(b); t = nb.nest[0]['pos'].pop(x); nb.journey(t, 0, unit_dir(dirs[T]))
    return Counter(ST)
if __name__ == '__main__':
    NW, NCH = int(sys.argv[1]), int(sys.argv[2]); tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(work, [(w, NCH) for w in range(NCH)]): tot.update(o)
    for k in sorted(tot, key=str): print(tot[k], k)
