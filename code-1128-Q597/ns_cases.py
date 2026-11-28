"""PROVE 2026-11-28: the case split of Theorem (NS!) on the census, and engine identification.
Per pre-step state: Case I (M interior to T_n: all 6 unit triangles at M are in T_n) or Case II;
 Case II: is R admissibly packed in nest A (engine packed_at + ideal)?  #candidates (S1)?  is the gash on side A?
 nest-rhombus EDGE-neighbour of R  =>  Case II ?
 engine identification: the tile set the engine actually removes == the unique good candidate of nslib.classify.
usage: python3 ns_cases.py NW NCH"""
import sys, copy
sys.path.insert(0, '/home/clio/projects/proofs/code-1128-Q597')
exec(open('/home/clio/projects/proofs/code-1128-Q597/ns_census.py').read().split("if __name__")[0])
def around(p):
    out = []
    for i in range(6):
        a, b = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)][i], [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)][(i+1) % 6]
        out.append(frozenset((p, (p[0]+a[0], p[1]+a[1]), (p[0]+b[0], p[1]+b[1]))))
    return out
def rec(self, R, d):
    if not CTX.get('on'): return _step(self, R, d)
    ok, P, M, Q, u, base = NL.gash(self, R)
    case = 'I' if all(t in CTX['Tn'] for t in around(M)) else 'II'
    e, res = NL.classify(self, R, d)
    nk = Counter(self.tiles[t][0] for t in self.neighbours(R))
    ST[('edge-nest-neighbour => case II', nk['rh'] == 0 or case == 'II')] += 1
    if case == 'II':
        A = self.nest[0]; ij = self.packed_at(A, self.tiles[R][1])
        adm = ij is not None and ij not in A['pos'] and self.is_ideal(set(A['pos']) | {ij})
        ST[('case II: admissibly packed in A; #cands; k==0', adm, len(res), CTX['k'] == 0)] += 1
    else:
        ST[('case I: #cands (good,bad)', sum(r[0] == 'good' for r in res), sum(r[0] == 'bad' for r in res))] += 1
    good = [r[2] for r in res if r[0] == 'good']
    before = dict(self.tiles); CTX['k'] += 1
    R2 = _step(self, R, d)
    removed = frozenset(t for t in before if t not in self.tiles)
    ST[('engine removed == the good candidate', len(good) == 1 and removed == good[0])] += 1
    return R2
migrate.Board.step = rec
if __name__ == '__main__':
    NW, NCH = int(sys.argv[1]), int(sys.argv[2]); tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(work, [(w, NCH) for w in range(NCH)]): tot.update(o)
    for k in sorted(tot, key=str): print(tot[k], k)
