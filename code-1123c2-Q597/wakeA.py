"""PROVE 11-23 c2: (W_a) in case a: A (toward T) removes, at steps k >= p+2, no tile that J_o (toward T') removes.
Placement of the contact exactly as caseA.py (Theorem local survivors, |S| = 3): J_o's steps c0, c0+1, c0+2 from
R^J_{c0}; A at R_p takes a_T, then a_T' (the square), so R^A_{p+2} = R_p + alpha + alpha'.
Part FR (J_o steps j >= c0, future + read steps): Lemma phi.  phi = a0 + b0, g(R) = phi(min R).
  T = B: A's tile at k: phi(t) = g(R^A_k) >= g(R^A_{p+2});  J (to C): phi(t) <= g(R^J_j) <= g(R^J_{c0}).
         so it suffices that g(R^J_{c0}) < g(R^A_{p+2}).
  T = C: A's tile: phi(t) <= g(R^A_{p+2});  J (to B): phi(t) = g(R^J_j) >= g(R^J_{c0}).  suffices g(R^J_{c0}) > g(R^A_{p+2}).
Part P (J_o steps j < c0, the past): window + words, as shpast.py, with reference rhombi R^A_{p+2} (A forward) and
  R^J_{c0} (J backward).  Feasible tiles in a box; D_A, D_J; then all T-words of length <= 2D_A+2 from R^A_{p+2} and all
  T'-words of length <= 2D_J+2 ending at R^J_{c0} (class of the c0 step); the removed-tile unions must be disjoint.
Controls: (i) the g-gap itself is printed; LESION: drop A's step p+1 (reference R^A_{p+1}): the FR gap must close.
(ii) liveness: number of the 12 unit translates of S_A meeting S_J.   usage: python3 wakeA.py [box]"""
import os, sys, itertools
from fractions import Fraction as Fr
BOX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
STR = {'B': 'a1', 'C': 'a0'}
CASEA = [('B', (0, 'a0'), (0, 'b'), (0, 0, 0, 0)), ('B', (1, 'b'), (0, 'b'), (0, 0, 1, -1)),
         ('C', (0, 'a1'), (0, 'b'), (0, 1, 0, 1)), ('C', (1, 'b'), (0, 'b'), (1, 0, 0, 0))]
def g(R): return min(R[1])[0] + min(R[1])[1]
def phi(t): return min(t[1])[0] + min(t[1])[1]
def Q(T):
    return {offs(t, Rc, T) for (ck, typ), (Rc, rem, Rn, cre) in CF[T].items() for t in rem}
QS = {T: Q(T) for T in 'BC'}
shapes = set()
for T in 'BC':
    for (ck, typ), (Rc, rem, Rn, cre) in CF[T].items():
        for t in rem | cre: shapes.add(trs(t, tuple(-x for x in base(t))))
def l1(a, b): return sum(abs(x - y) for x, y in zip(a, b))
def place(T, t0, off):
    Tp = OTHER[T]
    RJ0 = canon(Tp, t0[0])
    RAc = canon(T, 1); sh = tuple(b + o - a for a, b, o in zip(base(RAc), base(RJ0), off)); RAp = trs(RAc, sh)
    Hp, RAp1, Cp = step(T, RAp, STR[T]); Hq, RAp2, Cq = step(T, RAp1, STR[Tp])
    return RJ0, RAp, RAp1, RAp2, Hp | Hq
def window(T, Tp, RA, RJ):
    DA = DJ = 0; nf = ob = 0
    for s in shapes:
        for v in itertools.product(range(-BOX, BOX + 1), repeat=4):
            t = trs(s, tuple(a + b for a, b in zip(base(RA), v)))
            oA, oJ = offs(t, RA, T), offs(t, RJ, Tp)
            qa = [q for q in QS[T] if oA[0] <= q[0] and oA[1] <= q[1] and oA[2] >= q[2] and oA[3] >= q[3]]
            qj = [q for q in QS[Tp] if oJ[0] >= q[0] and oJ[1] >= q[1] and oJ[2] <= q[2] and oJ[3] <= q[3]]
            if not qa or not qj: continue
            nf += 1; ob += max(abs(x) for x in v) == BOX
            DA = max(DA, max(l1(oA, q) for q in qa)); DJ = max(DJ, max(l1(oJ, q) for q in qj))
    return DA, DJ, nf, ob
def run_word(T, R, word):
    out = []
    for typ in word:
        Rc, rem, Rn, cre = CF[T][(classkey(R), typ)]
        sh = base(R); assert trs(Rc, sh) == R
        out.append(frozenset(trs(t, sh) for t in rem)); R = trs(Rn, sh)
    return out, R
def SA(T, RA, LA):
    S = set()
    for n in range(1, LA + 1):
        for w in itertools.product(('a0', 'a1', 'b'), repeat=n):
            try: steps, _ = run_word(T, RA, w)
            except KeyError: continue
            for r in steps: S |= r
    return S
def SJ(Tp, RJ, LJ):
    S = set(); nw = 0
    for c in {k[0] for k in CF[Tp]}:
        R0 = CF[Tp][(c, 'a0')][0]
        for n in range(1, LJ + 1):
            for w in itertools.product(('a0', 'a1', 'b'), repeat=n):
                try: steps, Rend = run_word(Tp, R0, w)
                except KeyError: continue
                if classkey(Rend) != classkey(RJ): continue
                sh = tuple(a - b for a, b in zip(base(RJ), base(Rend))); nw += 1
                for r in steps: S |= {trs(t, sh) for t in r}
    return S, nw
nr = lambda S: {u for u in S if u[0] != 'rh'}
for (T, t0, t1, off) in CASEA:
    Tp = OTHER[T]
    RJ0, RAp, RAp1, RAp2, Hmid = place(T, t0, off)
    gJ, gA2, gA1, gAp = g(RJ0), g(RAp2), g(RAp1), g(RAp)
    okFR = gJ < gA2 if T == 'B' else gJ > gA2
    lesFR = gJ < gA1 if T == 'B' else gJ > gA1
    print(T, t0, t1, 'g(R^J_c0)', gJ, 'g(R_p)', gAp, 'g(R_p+1)', gA1, 'g(R_p+2)', gA2, '| FR separates:', okFR,
          '| LESION ref R_p+1 separates:', lesFR)
    DA, DJ, nf, ob = window(T, Tp, RAp2, RJ0)
    LA, LJ = 2 * DA + 2, 2 * DJ + 2
    Sa = nr(SA(T, RAp2, LA)); Sj, nw = SJ(Tp, RJ0, LJ); Sj = nr(Sj)
    live = sum(1 for v in W12 if {trs(t, v) for t in Sa} & Sj)
    print('   P: feasible', nf, 'on box bd', ob, 'D_A', DA, 'D_J', DJ, 'LA', LA, 'LJ', LJ, '|S_A|', len(Sa), 'J words', nw,
          '|S_J|', len(Sj), 'S_A & S_J', len(Sa & Sj), '| liveness translates meeting', live, '/12',
          '| CONTROL S_J & mid hexagons', len(Sj & nr(Hmid)))
