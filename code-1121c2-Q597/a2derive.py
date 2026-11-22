"""PROVE 11-21 c2: (A2) from (CC), the finite part.  In each (CC) configuration, list the non-rhombus tiles of A's
hexagons H_p and H_{p+1} that are J_o tiles, with J_o's step (c0 / c0+1), its sheet_pre, and whether the tile is
A-SHAPED (-alpha = -delta_T parallel to an edge; delta_B = w_A = (0,0;0,1), delta_C = u_A = (0,-1;0,0)).
(A2) can only fail at an A-shaped tile read at a step >= p+1 with sheet_pre = 1, i.e. the c0 tile in H_{p+1} when c0 is (1,b).
Control: the same tiles tested with the OTHER frame vector (delta_T'), which must make some of them A-shaped.
usage: python3 a2derive.py"""
import sys
exec(open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split("for cc in CCS:")[0])
DELTA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
def edges(t):
    V = list(t[1]); out = set()
    for P in V:
        for Q in V:
            d = tuple(Fr(b) - Fr(a) for a, b in zip(P, Q))
            if any(d): out.add(d)
    return out
def parallel(t, v):
    for d in edges(t):
        ks = {Fr(a) / b for a, b in zip(d, v) if b != 0}
        if len(ks) == 1 and all(a == 0 for a, b in zip(d, v) if b == 0) and next(iter(ks)) != 0: return True
    return False
for T, t0, t1, tA, off in CCS:
    Tp = OTHER[T]
    RJ0 = canon(Tp, t0[0]); H0, RJ1, C0 = step(Tp, RJ0, t0[1]); H1, RJ2, C1 = step(Tp, RJ1, t1[1])
    RAc = canon(T, tA[0]); sh = tuple(b + o - a for a, b, o in zip(base(RAc), base(RJ0), off)); RAp = trs(RAc, sh)
    Hp, RAp1, Cp = step(T, RAp, tA[1]); Hb, Rb, Cb = step(T, RAp1, 'b')
    mal = tuple(-x for x in DELTA[T]); malc = tuple(-x for x in DELTA[Tp])
    sp = {0: 1 if t0 == (1, 'b') else 0, 1: 0}
    for name, H, k in (('H_p', Hp, 'p'), ('H_p+1', Hb, 'p+1')):
        for u in H:
            if u[0] == 'rh': continue
            src = 'c0' if u in H0 else ('c0+1' if u in H1 else 'not J')
            if src == 'not J': continue
            s = sp[0] if src == 'c0' else sp[1]
            print(T, t0, '| tile of', name, u[0], 'from J step', src, 'sheet_pre', s, '| A-shaped (-delta_T):', parallel(u, mal),
                  '| control (-delta_T\'):', parallel(u, malc), '| (A2) threat' if (k == 'p+1' and s == 1 and parallel(u, mal)) else '')
