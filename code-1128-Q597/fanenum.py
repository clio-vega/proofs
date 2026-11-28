"""PROVE 2026-11-28, Lemma Fan (finite, exhaustive, board-free).
For T in {B,C} and each of the 6 unit vectors u: R = [0,U0]x[0,U1], U0=(u,0), U1=(0,u).  At each obtuse corner v of R,
enumerate EVERY edge-to-edge filling of the complementary angle at v by tiles of the non-rhombus forms
(0-triangle, 1-triangle, 90-degree square e x f), by DFS around v.  For each filling F, build the 4-tile board R+F and run
the ENGINE's candidate test (nslib.cands = verbatim copy of Board.step's S1 part) and the S2 across edge.
Prints per (T, u, corner): #fillings, and per filling (types, accepted?, hexagon type, x-rule, good/bad)."""
import sys, cmath, math
sys.path.insert(0, '/home/clio/projects/proofs/code-1128-Q597')
import nslib as NL, migrate
from collections import Counter
G = migrate.G
UNIT = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
VEC = [(a, b, 0, 0) for a, b in UNIT] + [(0, 0, a, b) for a, b in UNIT]
def add(P, Q): return tuple(p + q for p, q in zip(P, Q))
def neg(P): return tuple(-p for p in P)
def ang(v): return math.degrees(cmath.phase(G.phi(v))) % 360
def ccw(a, b):   # ccw angle from direction a to direction b, in (0,360)
    t = (ang(b) - ang(a)) % 360; return t
def tile_kind(c, n):
    """tile at v with edges c, n (n reached from c by a ccw turn of the tile's angle) or None."""
    l = lambda e: 0 if not any(e[2:]) else 1
    t = ccw(c, n)
    if l(c) == l(n) and abs(t - 60) < 1e-6: return 'tri'
    if l(c) != l(n) and abs(t - 90) < 1e-6: return 'sq'
    return None
def fillings(v, e_from, e_to):
    """all sequences of tiles filling the ccw sector at v from edge e_from to edge e_to."""
    total = ccw(e_from, e_to); out = []
    def dfs(cur, acc, seq):
        if abs(acc - total) < 1e-6:
            if cur == e_to: out.append(list(seq))
            return
        for n in VEC:
            k = tile_kind(cur, n)
            if k is None: continue
            a = acc + ccw(cur, n)
            if a > total + 1e-6: continue
            seq.append((k, cur, n)); dfs(n, a, seq); seq.pop()
    dfs(e_from, 0.0, []); return out
def tile_vs(v, k, c, n):
    if k == 'tri': return (v, add(v, c), add(v, n))
    return (v, add(v, c), add(add(v, c), n), add(v, n))
def board(tiles):
    b = migrate.Board.__new__(migrate.Board); b.tiles = {}; b.edges = {}; b.nid = 0
    ids = [b.add(k, vs) for k, vs in tiles]; return b, ids
D = {'B': G.phi((0, 0, 0, 1)), 'C': -G.phi((0, 1, 0, 0))}
for T in D: D[T] /= abs(D[T])
def eta(T, u): return u[1] if T == 'B' else u[0] + u[1]
if __name__ == '__main__':
    summary = Counter()
    for T in 'BC':
        for u in UNIT:
            U0, U1 = (u[0], u[1], 0, 0), (0, 0, u[0], u[1]); O = (0, 0, 0, 0)
            Rvs = (O, U0, add(U0, U1), U1)
            allowed = eta(T, u) != 0
            # orient R ccw in the plane
            z = [G.phi(P) for P in Rvs]; area = sum((z[i].conjugate()*z[(i+1) % 4]).imag for i in range(4))
            if area < 0: Rvs = tuple(reversed(Rvs))
            n_ok_corners = []
            for v in (U0, U1):
                i = Rvs.index(v); prev, nxt = Rvs[i-1], Rvs[(i+1) % 4]
                # R occupies the ccw sector at v from (nxt - v) to (prev - v); the fan is from (prev-v) ccw to (nxt-v)
                ef, et = NL.vsub(prev, v), NL.vsub(nxt, v)
                Fs = fillings(v, ef, et)
                assert abs(ccw(ef, et) - 210) < 1e-6
                rows = []
                for F in Fs:
                    tiles = [('rh', Rvs)] + [(k, tile_vs(v, k, c, n)) for k, c, n in F]
                    b, ids = board(tiles)
                    e, res = NL.classify(b, ids[0], D[T])
                    full = [r for r in res if r[2] == frozenset(ids)]
                    # x-rule-free acceptance: rerun with a direction-free test (cands with x-rule disabled)
                    acc = NL.cands_nox(b, ids[0]) if hasattr(NL, 'cands_nox') else None
                    word = ''.join('t' if k == 'tri' else 's' for k, _, _ in F)
                    verdict = full[0][0] if full else 'x-fail'
                    typ = full[0][1] if full else '-'
                    acc_ok = any(S == frozenset(ids) for _, S, _ in acc) if acc is not None else None
                    rows.append((word, acc_ok, typ, verdict, len(res)))
                    summary[(T, allowed, 'accepted-hexagon (no x-rule)', acc_ok)] += 1
                n_ok_corners.append(rows)
                print(T, u, 'allowed' if allowed else 'FORBIDDEN', 'corner', 'U0' if v == U0 else 'U1',
                      '#fillings', len(Fs), rows, 'across', 'defined' if NL.across_edge(board([('rh', Rvs)])[0], 0, D[T]) else 'undef')
            if allowed:
                vs = [[r[3] for r in rows] for rows in n_ok_corners]
                goodc = [all(x == 'good' for x in c) for c in vs]
                otherc = [all(x in ('bad', 'x-fail') for x in c) for c in vs]
                ok = (goodc[0] and otherc[1]) or (goodc[1] and otherc[0])
                summary[(T, 'allowed', 'one corner all good, other all bad/x-fail', ok)] += 1
                summary[(T, 'allowed', '#fillings per corner', tuple(len(r) for r in n_ok_corners))] += 1
    print()
    for k in sorted(summary, key=str): print(summary[k], k)
