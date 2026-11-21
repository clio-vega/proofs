"""PROVE 11-21 c2: (Sh3-past) as a finite check.  Claim: after the split, A (toward T) never removes a tile that J_o
(toward T') removed at a step j < e (before the split).
Step 1 (window): for such a tile t, with R* the split rhombus (class 0, shared), monotonicity of the 4-offset
f(i) = om(t, R_i) along each journey (Lemma window, G1c tex) gives
  A:  om_T(t,R*) psi <= q_psi, lam >= q_lam  for some q in Q_T;     k - split <= 2 max_q |om_T(t,R*) - q|_1 + 1
  J:  om_T'(t,R*) psi >= q'_psi, lam <= q'_lam for some q' in Q_T';  e - j  <= 2 max_q' |om_T'(t,R*) - q'|_1 + 1
We enumerate every tile (all shapes occurring in the configuration table) at every translate in a box, keep the
feasible ones, and take the max of the two bounds (box adequacy: no feasible tile on the box boundary).
Step 2 (words): S_J = tiles removed at steps < e by every T'-word of length <= LJ whose last step is the class-0 b at R*;
S_A = tiles removed after the split by every T-word of length <= LA continuing from the T-image of R*.
Words are independent in the model, so a realisable violation needs S_J & S_A nonempty.   usage: python3 shpast.py [box]"""
import os, sys, itertools
from collections import Counter
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '.')
exec(open('words.py').read().split('\ndef run(T, L, st):')[0])
BOX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
k0 = frozenset({(0, 1), (0, -1)})
CF = {T: configs(T) for T in 'BC'}
def Q(T):
    out = set()
    for (ck, typ), (Rc, rem, Rn, cre) in CF[T].items():
        for t in rem: out.add(offs(t, Rc, T))
    return out
QS = {T: Q(T) for T in 'BC'}
print('Q_B', sorted(QS['B'])); print('Q_C', sorted(QS['C']))
Rstar = CF['B'][(k0, 'b')][0]; assert Rstar == CF['C'][(k0, 'b')][0]
shapes = set()
for T in 'BC':
    for (ck, typ), (Rc, rem, Rn, cre) in CF[T].items():
        for t in rem | cre:
            sh = tuple(-x for x in base(t)); shapes.add(trs(t, sh))
print('tile shapes', len(shapes))
def l1(a, b): return sum(abs(x - y) for x, y in zip(a, b))
def window(T, Tp):
    DA = DJ = 0; nfeas = 0; onbd = 0
    for s in shapes:
        for v in itertools.product(range(-BOX, BOX + 1), repeat=4):
            t = trs(s, v)
            oA, oJ = offs(t, Rstar, T), offs(t, Rstar, Tp)
            qa = [q for q in QS[T] if oA[0] <= q[0] and oA[1] <= q[1] and oA[2] >= q[2] and oA[3] >= q[3]]
            qj = [q for q in QS[Tp] if oJ[0] >= q[0] and oJ[1] >= q[1] and oJ[2] <= q[2] and oJ[3] <= q[3]]
            if not qa or not qj: continue
            nfeas += 1
            if max(abs(x) for x in v) == BOX: onbd += 1
            DA = max(DA, max(l1(oA, q) for q in qa)); DJ = max(DJ, max(l1(oJ, q) for q in qj))
    return DA, DJ, nfeas, onbd
def run_word(T, R, word):
    """yield (step index, R before, removed tiles) for word from R"""
    out = []
    for i, typ in enumerate(word):
        Rc, rem, Rn, cre = CF[T][(classkey(R), typ)]
        sh = base(R); assert trs(Rc, sh) == R
        out.append((i, R, frozenset(trs(t, sh) for t in rem))); R = trs(Rn, sh)
    return out, R
def SJ(Tp, LJ):
    """tiles removed before the last step by T'-words of length <= LJ ending in the class-0 b at R* (translated)"""
    S = set(); nwords = 0
    for c in {k[0] for k in CF[Tp]}:
        R0 = CF[Tp][(c, 'a0')][0]
        for n in range(0, LJ):
            for w in itertools.product(('a0', 'a1', 'b'), repeat=n):
                steps, Rend = run_word(Tp, R0, w)
                if classkey(Rend) != k0: continue
                sh = tuple(a - b for a, b in zip(base(Rstar), base(Rend)))
                nwords += 1
                for i, Rb, rem in steps: S |= {trs(t, sh) for t in rem}
    return S, nwords
def SA(T, LA):
    Rc, rem, Rn, cre = CF[T][(k0, 'b')]
    sh = base(Rstar); R1 = trs(Rn, sh); split = frozenset(trs(t, sh) for t in rem)
    S = set(); nwords = 0
    for n in range(1, LA + 1):
        for w in itertools.product(('a0', 'a1', 'b'), repeat=n):
            steps, _ = run_word(T, R1, w); nwords += 1
            for i, Rb, r in steps: S |= r
    return S, nwords, split
if __name__ == '__main__':
    for T, Tp in (('B', 'C'), ('C', 'B')):
        DA, DJ, nf, ob = window(T, Tp)
        LA, LJ = 2 * DA + 1, 2 * DJ + 2
        print(T, 'feasible tiles', nf, 'on box boundary', ob, 'D_A', DA, 'D_J', DJ, '-> LA', LA, 'LJ (letters incl. split)', LJ)
        Sj, nj = SJ(Tp, LJ); Sa, na, split = SA(T, LA)
        inter = Sj & Sa
        print(T, 'J words', nj, '|S_J|', len(Sj), 'A words', na, '|S_A|', len(Sa), 'S_J & S_A', len(inter),
              'CONTROL split hexagon & S_A', len(split & Sa), 'CONTROL S_J & split', len(Sj & split))
        for t in list(inter)[:5]: print('   violation tile', t)
