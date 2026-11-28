"""Gash paths (puzzle coords) of K,H,G,Gp for (SW) instances N0..N1.  A gash = psi-image of the moving rhombus,
printed as (midpoint, direction) with midpoint in lattice coords (a,b) meaning a + b*omega6."""
import sys
exec(open('/home/clio/projects/proofs/code-1128c2-Q585-SW/swfull.py').read().split("if __name__")[0])
POS = []
_s3 = migrate.Board.step
def rec(self, R, d):
    R2 = _s3(self, R, d)
    if R2 is not None: POS.append(self.tiles[R2][1])
    return R2
migrate.Board.step = rec
def gash(vs):
    ps = sorted(set((P[0]+P[2], P[1]+P[3]) for P in vs))
    return ps  # 3 points: ends + middle
def path(b, ij, dv):
    POS.clear(); IBC.STEPLOG.clear(); LOG.clear()
    t = b.nest[0]['pos'].pop(ij); v0 = b.tiles[t][1]; R, Y, ij2, s = b.journey(t, 0, unit_dir(dv))
    return [gash(v0)] + [gash(v) for v in POS], ij2
N0, N1 = int(sys.argv[1]), int(sys.argv[2]); cnt = 0
for key, b0 in census():
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            y = (x[0], x[1]-1)
            if y not in S or y not in corners(S - {x}): continue
            if cnt < N0: cnt += 1; continue
            if cnt >= N1: sys.exit()
            cnt += 1
            b1 = copy.deepcopy(b); K, lb = path(b1, x, dirs['B']); G, lc = path(b1, y, dirs['C'])
            b2 = copy.deepcopy(b); H, lc2 = path(b2, x, dirs['C']); Gp, lb2 = path(b2, y, dirs['B'])
            print('#', cnt-1, key[2], key[3], 'x', x, 'y', y, 'boxes', lb, lc, lc2, lb2, 'same', tiles_of(b1) == tiles_of(b2))
            for nm, P in (('K x>B', K), ('H x>C', H), ('G y>C', G), ('Gp y>B', Gp)):
                print('  %-7s' % nm, ' '.join(''.join(str(c) for c in (p[0][0], p[0][1])) + '-' + ''.join(str(c) for c in (p[2][0], p[2][1])) for p in P))
