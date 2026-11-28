"""PROVE 2026-11-28 lesions for the (NS) proof, on a 1-in-5 sample of size-5 closure states (census order).
 CTRL : unlesioned journeys (all corners, B and C): stuck / raised count (must be 0).
 L1   : start at a NON-corner nest-A cell (breaks 'nest A minus x is an ideal', used in Case II at v_in):
        k=0 #candidates, and whether a nest rhombus lies across a FRONT edge of R; outcome.
 L2   : delete one tri/sq tile t of the board (breaks (ii) 'tiles cover H'); journeys of every corner: outcome
        (terminates / stuck / raised), and for stuck ones whether R's obtuse-corner fans touch the hole.
 L4   : x-rule off: #accepted fan configurations per pre-step state (proof predicts 2 at gash-interior states)."""
import sys, copy, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1128-Q597')
sys.path.insert(0, '/home/clio/projects/reviews/code-1125c2-BS')
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590')
import nslib as NL
import indep_bs_census as IBC
from indep_bs_census import census, closure
from common import corners
from flock import unit_dir
import migrate
from collections import Counter
from multiprocessing import Pool
_step = IBC._engine_step
ST = Counter(); MODE = {}
def rec(self, R, d):
    if MODE.get('l4'):
        n = len(NL._cands_generic(self, R)); ST[('L4 x-rule off: #accepted configs', n)] += 1
    return _step(self, R, d)
migrate.Board.step = rec
def outcome(nb, t, d):
    try: nb.journey(t, 0, d); return 'terminates'
    except RuntimeError as ex: return 'RuntimeError: ' + str(ex)[:30]
    except AssertionError as ex: return 'AssertionError'
def work(args):
    w, NCH = args; ST.clear(); rng = random.Random(1128 + w)
    for i, (key, b0) in enumerate(census()):
        if i % NCH != w: continue
        MODE['l4'] = False; cl = list(closure(b0))
        for si, (b, dirs) in enumerate(cl):
            if si % 5: continue
            pos = set(b.nest[0]['pos']); cs = corners(pos)
            for T in 'BC':
                d = unit_dir(dirs[T])
                for x in cs:
                    MODE['l4'] = True
                    nb = copy.deepcopy(b); t = nb.nest[0]['pos'].pop(x); ST[('CTRL', outcome(nb, t, d))] += 1
                    MODE['l4'] = False
                for x in pos - set(cs):
                    nb = copy.deepcopy(b); t = nb.nest[0]['pos'].pop(x)
                    nk = NL.cands(nb, t, d)
                    fr_rh = 'rh' in [nb.tiles[s][0] for s in nb.neighbours(t)]
                    ST[('L1 non-corner: k=0 #cands, rh neighbour', len(nk), fr_rh)] += 1
                    ST[('L1 outcome', outcome(nb, t, d))] += 1
                # L2: one random tri/sq hole per (state, T)
                nb = copy.deepcopy(b)
                hole = rng.choice([s for s, (k, vs) in nb.tiles.items() if k != 'rh'])
                hv = set(nb.tiles[hole][1]); nb.remove(hole)
                for x in cs:
                    nb2 = copy.deepcopy(nb); t = nb2.nest[0]['pos'].pop(x)
                    ST[('L2 hole: outcome', outcome(nb2, t, d))] += 1
    return Counter(ST)
if __name__ == '__main__':
    NW, NCH = int(sys.argv[1]), int(sys.argv[2]); tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(work, [(w, NCH) for w in range(NCH)]): tot.update(o)
    for k in sorted(tot, key=str): print(tot[k], k)
