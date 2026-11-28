"""PROVE 2026-11-28: for every good x-passing fan configuration of fanenum.py (allowed classes), run the ENGINE's own
Board.step on the 4-tile board: it must not raise (in particular the type-b rotation assert has exactly one k) and the
new rhombus must again have an allowed class (W1)."""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1128-Q597')
from fanenum import *
cnt = Counter()
for T in 'BC':
    for u in UNIT:
        if eta(T, u) == 0: continue
        U0, U1 = (u[0], u[1], 0, 0), (0, 0, u[0], u[1]); O = (0, 0, 0, 0)
        Rvs = (O, U0, add(U0, U1), U1)
        z = [G.phi(P) for P in Rvs]
        if sum((z[i].conjugate()*z[(i+1) % 4]).imag for i in range(4)) < 0: Rvs = tuple(reversed(Rvs))
        for v in (U0, U1):
            i = Rvs.index(v); ef, et = NL.vsub(Rvs[i-1], v), NL.vsub(Rvs[(i+1) % 4], v)
            for F in fillings(v, ef, et):
                tiles = [('rh', Rvs)] + [(k, tile_vs(v, k, c, n)) for k, c, n in F]
                b, ids = board(tiles)
                e, res = NL.classify(b, ids[0], D[T])
                if not res or res[0][0] != 'good': continue
                try:
                    R2 = migrate.Board.step(b, ids[0], D[T]); vs = b.tiles[R2][1]
                    ok, P, M, Q, u2, base = NL.gash(b, R2)
                    cnt[(T, 'engine step ok; new class allowed', ok and eta(T, u2) != 0)] += 1
                except Exception as ex:
                    cnt[(T, 'engine step RAISED', str(ex)[:40])] += 1
for k in sorted(cnt, key=str): print(cnt[k], k)
