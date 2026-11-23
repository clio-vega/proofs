"""PROVE 11-23 c2, brief step 3 lesion: drop the survivor clause (|H_p & W| = 3, J_o's c0..c0+2 read).
Placements: (a) the |S| = 2 placements (c0+2 type not a survivor), (b) the contact offset shifted by each of 12 unit vectors.
On each placeable lesioned contact evaluate the three finite ingredients of (Stale_a):
  DC  the D-composition claim (re-run-only tiles in Jcre or C'_p);
  FR  g-separation of J's steps >= c0 (T=B) / >= c0+2 (T=C) from A's steps >= p+2, plus the one-sided window on leftovers;
  P   the past window from (R^A_{p+2}, R^J_{c0}) (box 4): number of feasible tiles.
Also (Ea)'s ingredient, which is survivor-free by design: printed for completeness.   usage: python3 lesionA.py"""
import sys
from collections import Counter
sys.argv = ['wakeA.py', '4']
exec(open('wakeA.py').read().split('\nnr = lambda')[0])
nr = lambda S: {u for u in S if u[0] != 'rh'}
def winA(T, t, R):
    oA = offs(t, R, T)
    return [q for q in QS[T] if oA[0] <= q[0] and oA[1] <= q[1] and oA[2] >= q[2] and oA[3] >= q[3]]
def evaluate(T, t0, t1, t2, off):
    Tp = OTHER[T]
    RJ0 = canon(Tp, t0[0]); H0, RJ1, C0 = step(Tp, RJ0, t0[1]); H1, RJ2, C1 = step(Tp, RJ1, t1[1]); H2, RJ3, C2 = step(Tp, RJ2, t2)
    RAc = canon(T, 1); sh = tuple(b + o - a for a, b, o in zip(base(RAc), base(RJ0), off)); RAp = trs(RAc, sh)
    Hp, RAp1, Cp = step(T, RAp, STR[T]); Hp1, RAp2, Cp1 = step(T, RAp1, STR[Tp])
    Hq, RBp1, Cq = step(T, RAp, STR[Tp]); Hq1, RBp2, Cq1 = step(T, RBp1, STR[T])
    Jrem = nr(H0) | nr(H1) | nr(H2); Jcre = nr(C0) | nr(C1) | nr(C2)
    S = len(nr(Hp) & Jrem)
    rest = (nr(Hq) | nr(Hq1)) - nr(Hp) - nr(Hp1)
    DC = rest <= (Jcre | nr(Cq))
    gs = [g(R) for R in (RJ0, RJ1, RJ2, RJ3)]; j0 = 0 if T == 'B' else 2
    sepg = (gs[j0] < g(RAp2)) if T == 'B' else (gs[j0] > g(RAp2))
    left = [u for u in (Jrem if T == 'B' else (nr(H0) | nr(H1))) - nr(Hp) - nr(Hp1)]
    FR = sepg and all(not winA(T, u, RAp2) for u in left)
    _, _, nf, _ = window(T, Tp, RAp2, RJ0)
    return S, DC, FR, nf
st = Counter()
for (T, t0, t1, off) in CASEA:
    for t2 in ('a0', 'a1', 'b'):
        for tag, offs_ in (('true offset', [off]), ('shifted offset', [tuple(a + b for a, b in zip(off, v)) for v in W12])):
            for o in offs_:
                try: S, DC, FR, nf = evaluate(T, t0, t1, t2, o)
                except (KeyError, AssertionError): st[(T, tag, 'not placeable')] += 1; continue
                st[(T, tag, '|S|', S, 'DC', DC, 'FR', FR, 'P feasible==0', nf == 0, 'all three', DC and FR and nf == 0)] += 1
for k in sorted(st, key=str): print(st[k], k)
