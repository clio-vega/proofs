"""PROVE 2026-11-28 Step 0: (NS) on the uncapped size-5 closure (all states, all nest-A corners, both targets).
At every call of Board.step (= every non-packed pre-step state):
  N2 #good candidates, N3 #bad (+ is it the undo of the previous step; k = step index), N4 nest-rhombus edge-neighbours,
  G1 tri/sq psi-images tile T_n (multiset of unit triangles == that of the census start board),
  G2 labels agree across interior puzzle edges except the gash; record the gash side labels.
usage: python3 ns_census.py NW NCH"""
import sys, copy
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, '/home/clio/projects/proofs/code-1128-Q597')
sys.path.insert(0, '/home/clio/projects/reviews/code-1125c2-BS')
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590')
import nslib as NL
import indep_bs_census as IBC
from indep_bs_census import census, closure
from common import corners
from flock import unit_dir
import migrate
_step = IBC._engine_step
ST = Counter(); CTX = {}
def tn_tris(b):
    out = Counter()
    for t, (kind, vs) in b.tiles.items():
        if kind == 'rh': continue
        for u in NL.unit_tris(kind, vs)[0]: out[u] += 1
    return out
def puzzle_checks(self, R):
    ok, P, M, Q, u, base = NL.gash(self, R)
    ST[('G0 R collapses to a length-2 segment', ok)] += 1
    if not ok: return
    tri = Counter(); edgelab = {}
    for t, (kind, vs) in self.tiles.items():
        if kind == 'rh': continue
        ts, es = NL.unit_tris(kind, vs)
        for x in ts: tri[x] += 1
        for e, l in es.items(): edgelab.setdefault(e, []).append(l)
    ST[('G1 tri/sq images tile T_n', tri == CTX['Tn'])] += 1
    gas = {frozenset((P, M)), frozenset((M, Q))}
    bad = [e for e, ls in edgelab.items() if len(ls) == 2 and ls[0] != ls[1] and e not in gas]
    ST[('G2 labels agree off the gash', not bad)] += 1
    for e in gas:
        ST[('G2 gash half-edge: #pieces, labels', len(edgelab.get(e, [])), tuple(sorted(map(str, edgelab.get(e, [])))))] += 1
def rec(self, R, d):
    if not CTX.get('on'): return _step(self, R, d)
    k = CTX['k']; CTX['k'] += 1
    e, res = NL.classify(self, R, d)
    ng = sum(1 for r in res if r[0] == 'good'); nb = sum(1 for r in res if r[0] == 'bad')
    ST[('N2 #good', ng)] += 1
    ST[('N3 #bad, k==0', nb, k == 0)] += 1
    if nb:
        und = all(r[2] == CTX['last_created'] for r in res if r[0] == 'bad')
        ST[('N3 bad == undo (tile ids of previous step)', und)] += 1
    nk = Counter(self.tiles[t][0] for t in self.neighbours(R))
    ST[('N4 R has a nest-rhombus edge-neighbour; #good', nk['rh'] > 0, ng)] += 1
    ST[('N5 across edge defined; #good, #bad', e is not None, ng, nb)] += 1
    puzzle_checks(self, R)
    before = set(self.tiles)
    R2 = _step(self, R, d)
    CTX['last_created'] = frozenset(set(self.tiles) - before)
    return R2
migrate.Board.step = rec
def work(args):
    w, NCH = args; ST.clear()
    for i, (key, b0) in enumerate(census()):
        if i % NCH != w: continue
        CTX['on'] = False; cl = list(closure(b0)); CTX['on'] = True
        for b, dirs in cl:
            CTX['Tn'] = tn_tris(b); ST['states'] += 1
            for x in corners(set(b.nest[0]['pos'])):
                for T in 'BC':
                    nb = copy.deepcopy(b); t = nb.nest[0]['pos'].pop(x); CTX['k'] = 0; CTX['last_created'] = None
                    nb.journey(t, 0, unit_dir(dirs[T])); ST[('journeys', T)] += 1
        CTX['on'] = False   # census() runs migrate_flock (Board.step) for the NEXT key: keep it out
    return Counter(ST)
if __name__ == '__main__':
    NW, NCH = int(sys.argv[1]), int(sys.argv[2]); tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(work, [(w, NCH) for w in range(NCH)]): tot.update(o)
    for k in sorted(tot, key=str): print(tot[k], k)
