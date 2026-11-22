"""PROVE 11-22 c2: case a (L's re-run takes (1,a_T')), finite derivation on the configuration table (placement as hsderive.py).
Contacts (Theorem local survivors with |S| = 3): A at (1,a_T), re-run (1,a_T'), J_o steps c0, c0+1, c0+2.
A side: certified present on X_{p+1} = J's removed tiles (W; on X_p because the prefix reads no W) minus H_p, plus A's created
tiles at p.  Which class-1 T-configurations at R_{p+1} have all 3 non-rh tiles certified?  (expect: only a_T').
Re-run side: certified on Y_{p+1} = J's created tiles (never removed by J: no-revisit; in Y_p by Lemma dinv) minus the re-run's
step-p hexagon, plus the re-run's created tiles at p.  Which class-1 configurations at R'_{p+1} are certified?  (expect: only a_T).
Also: A's step p+1 (a_T' at R_{p+1}) == re-run step p + alpha, and re-run step p+1 == A step p + alpha' (the square closes).
Lesion: offsets shifted by the 12 unit vectors.   usage: python3 caseA.py"""
import os, sys
from fractions import Fraction as Fr
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
STR = {'B': 'a1', 'C': 'a0'}
CASEA = [('B', (0, 'a0'), (0, 'b'), (0, 0, 0, 0)), ('B', (1, 'b'), (0, 'b'), (0, 0, 1, -1)),
         ('C', (0, 'a1'), (0, 'b'), (0, 1, 0, 1)), ('C', (1, 'b'), (0, 'b'), (1, 0, 0, 0))]
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
W12 = [u + (0, 0) for u in UN] + [(0, 0) + u for u in UN]
nr = lambda S: {u for u in S if u[0] != 'rh'}
def run(T, t0, t1, t2typ, off):
    Tp = OTHER[T]; aT, aTp = STR[T], STR[Tp]
    RJ0 = canon(Tp, t0[0]); H0, RJ1, C0 = step(Tp, RJ0, t0[1]); H1, RJ2, C1 = step(Tp, RJ1, t1[1]); H2, RJ3, C2 = step(Tp, RJ2, t2typ)
    RAc = canon(T, 1); sh = tuple(b + o - a for a, b, o in zip(base(RAc), base(RJ0), off)); RAp = trs(RAc, sh)
    Hp, RAp1, Cp = step(T, RAp, aT)          # A at p: straight step of T
    Hq, RBp1, Cq = step(T, RAp, aTp)         # re-run at p: a_T'
    Jrem = nr(H0) | nr(H1) | nr(H2); Jcre = nr(C0) | nr(C1) | nr(C2)
    S = nr(Hp) & Jrem
    out = {'|S|': len(S), 'G_t\' meets Jcre': bool(nr(Hq) & Jcre)}
    Xknown = (Jrem - nr(Hp)) | nr(Cp)
    Yknown = (Jcre - nr(Hq)) | nr(Cq)
    out['A certified'] = tuple(t for t in ('a0', 'a1', 'b') if nr(step(T, RAp1, t)[0]) <= Xknown)
    out['rerun certified'] = tuple(t for t in ('a0', 'a1', 'b') if nr(step(T, RBp1, t)[0]) <= Yknown)
    al = tuple(b - a for a, b in zip(base(RAp), base(RAp1))); alp = tuple(b - a for a, b in zip(base(RAp), base(RBp1)))
    A1 = step(T, RAp1, aTp); B1 = step(T, RBp1, aT)
    out['A[p+1] == rerun[p] + alpha'] = (frozenset(trs(t, al) for t in Hq), trs(RBp1, al)) == (A1[0], A1[1])
    out['rerun[p+1] == A[p] + alpha\''] = (frozenset(trs(t, alp) for t in Hp), trs(RAp1, alp)) == (B1[0], B1[1])
    out['square closes'] = A1[1] == B1[1]
    return out
for (T, t0, t1, off) in CASEA:
    for t2 in ('a0', 'a1', 'b'):
        try: r = run(T, t0, t1, t2, off)
        except (KeyError, AssertionError) as ex: r = ('not placeable', type(ex).__name__)
        print(T, t0, t1, 'c0+2', t2, r)
print('LESION: shifted offsets with |S| = 3 and re-run meeting Jcre, and A certified == (a_T\',)')
for (T, t0, t1, off) in CASEA:
    for t2 in ('a0', 'a1', 'b'):
        n = 0; m = 0
        for v in W12:
            try:
                r = run(T, t0, t1, t2, tuple(a + b for a, b in zip(off, v)))
                if r['|S|'] == 3: n += 1; m += r['A certified'] == (STR[OTHER[T]],)
            except (KeyError, AssertionError): pass
        print(' ', T, t0, t1, t2, 'shifted offsets with |S|=3:', n, ' of which A certified only a_T\':', m)
