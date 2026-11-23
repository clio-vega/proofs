"""Retreat (case b): (Comp)(i) as a finite check at the split.  R* = class-0 rhombus where A (toward T) and J_o (toward T')
take their (0,b) step on the same hexagon (Lemma shared: identical hexagons).  R# = A's image (class 1 of T).
Candidate = the (1,a_T) good configuration of T at R# - alpha (alpha = displacement of T's (1,a_T)): this is H_p + gamma,
since gamma := R_{p+m+1} - R_{p+1} and R_{p+m+1} = R#.  Claim: its non-rh tiles are created by J_o at e (its (0,b) at R*)
or at e+1 (any placeable type at J_o's image R^o_{e+1}).  Then they are not in M0 (Lemma nr), so not in M0 & Mx.
LESION: alpha replaced by alpha' (the (1,a_T') displacement), and by 0.      usage: python3 compi.py"""
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
STR = {'B': 'a1', 'C': 'a0'}
nr = lambda S: {u for u in S if u[0] != 'rh'}
def disp_of(T, typ, cls=1):
    R = canon(T, cls); H, Rn, C = step(T, R, typ); return tuple(b - a for a, b in zip(base(R), base(Rn)))
for T in 'BC':
    Tp = OTHER[T]
    Rs = CF[T][(k0, 'b')][0]; assert Rs == CF[Tp][(k0, 'b')][0]
    HA, Rsharp, CA = step(T, Rs, 'b'); HJ, RJe1, CJ = step(Tp, Rs, 'b')
    assert HA == HJ
    for lab, al in (('alpha', disp_of(T, STR[T])), ('LESION alpha\'', disp_of(T, STR[Tp])), ('LESION 0', (0, 0, 0, 0))):
        Rc = trs(Rsharp, tuple(-x for x in al))
        try: Hc, Rcn, Cc = step(T, Rc, STR[T])
        except (KeyError, AssertionError): print(T, lab, 'not placeable'); continue
        for t1 in ('a0', 'a1', 'b'):
            try: HJ1, RJe2, CJ1 = step(Tp, RJe1, t1)
            except (KeyError, AssertionError): print(T, lab, 'J e+1', t1, 'not placeable'); continue
            cand = nr(Hc)
            print(T, lab, 'J e+1 type', t1, '| cand tiles', sorted(u[0] for u in cand), '| in Cre_J(e):', len(cand & nr(CJ)),
                  'in Cre_J(e+1):', len(cand & nr(CJ1)), 'CLAIM', cand <= nr(CJ) | nr(CJ1),
                  '| image R_c + alpha == R#:', Rcn == Rsharp, '| J e+1 hexagon meets cand:', len(nr(HJ1) & cand))
