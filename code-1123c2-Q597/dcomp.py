"""PROVE 11-23 c2, D-composition in case a (finite, on the placement of caseA.py).
Claim: every non-rh tile of the re-run's hexagons H'_p (a_T' at R_p) and H'_{p+1} (a_T at R'_{p+1}) that is not in
H_p or H_{p+1} (A's hexagons) lies in Jcre (J_o's created tiles at c0..c0+2) or in C'_p (the re-run's created tiles at p).
Then D = X_{p+2} - Y_{p+2} consists of wake tiles and of non-rh tiles created by A at p, p+1 only.
LESION: replace Jcre by Jrem (J_o's removed tiles): the claim must fail.     usage: python3 dcomp.py"""
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
STR = {'B': 'a1', 'C': 'a0'}
CASEA = [('B', (0, 'a0'), (0, 'b'), (0, 0, 0, 0)), ('B', (1, 'b'), (0, 'b'), (0, 0, 1, -1)),
         ('C', (0, 'a1'), (0, 'b'), (0, 1, 0, 1)), ('C', (1, 'b'), (0, 'b'), (1, 0, 0, 0))]
nr = lambda S: {u for u in S if u[0] != 'rh'}
for (T, t0, t1, off) in CASEA:
    Tp = OTHER[T]
    for t2 in ('a0', 'a1', 'b'):
        RJ0 = canon(Tp, t0[0]); H0, RJ1, C0 = step(Tp, RJ0, t0[1]); H1, RJ2, C1 = step(Tp, RJ1, t1[1])
        try: H2, RJ3, C2 = step(Tp, RJ2, t2)
        except KeyError: continue
        RAc = canon(T, 1); sh = tuple(b + o - a for a, b, o in zip(base(RAc), base(RJ0), off)); RAp = trs(RAc, sh)
        Hp, RAp1, Cp = step(T, RAp, STR[T]); Hp1, RAp2, Cp1 = step(T, RAp1, STR[Tp])
        Hq, RBp1, Cq = step(T, RAp, STR[Tp]); Hq1, RBp2, Cq1 = step(T, RBp1, STR[T])
        Jrem = nr(H0) | nr(H1) | nr(H2); Jcre = nr(C0) | nr(C1) | nr(C2)
        if len(nr(Hp) & Jrem) != 3: continue
        rest = (nr(Hq) | nr(Hq1)) - nr(Hp) - nr(Hp1)
        print(T, t0, t1, 'c0+2', t2, '| re-run tiles outside A hexagons:', len(rest), 'in Jcre:', len(rest & Jcre),
              'in C\'_p:', len(rest & nr(Cq)), 'CLAIM', rest <= (Jcre | nr(Cq)), '| H\'_p+1 & C\'_p (no-revisit, must be 0):', len(nr(Hq1) & nr(Cq)),
              '| R_p+2 == R\'_p+2', RAp2 == RBp2, '| LESION Jrem instead:', rest <= (Jrem | nr(Cq)))
