"""No stop at R^E_{k+3}, by contradiction through L^E (Theorem bsall: L^E_{k+1} = L_{k+2}-alpha, L^E_{k+2} = a+gamma =
T_L's (1,a_{T_L}) at R^L_{k+3}-alpha).  If J_E stopped at R^E_{k+3}, that rhombus would be a tile of M0.J_E that L^E never
removes; so no tile removed by L^E_{k+1} or L^E_{k+2} may overlap it.  Print the overlaps.
Also: is a+gamma really T_L's (1,a_TL) at R^L_{k+3}-alpha?  (check via the census in leprov.py / shape)."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
aT = {'B': 'a1', 'C': 'a0'}
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0)
    _, RE3, creE = step(TE, R, 'b')
    remLb_a, RL3a, creLb_a = step(TL, tr(R, mal), 'b')          # L^E_{k+1} = L_{k+2} - alpha
    assert RL3a == tr(step(TL, R, 'b')[1], mal)
    remA, RA, creA = step(TL, RL3a, aT[TL])                     # L^E_{k+2}
    for name, rem in (('L^E_{k+1}', remLb_a), ('L^E_{k+2}', remA)):
        ov = [g[0] for g in nonrh(rem) if overlap(g, RE3)]
        eq = [g[0] for g in rem if g == RE3]
        print('T_L', TL, name, 'removed nonrh tiles overlapping R^E_{k+3}:', ov, '| equal:', eq,
              '| overlapping nonrh C(E_{k+2}):', [g[0] for g in nonrh(rem) if any(overlap(g, h) and g != h for h in nonrh(creE))])
    # control: R^E_{k+3} vs the rhombi of L^E (must be distinct)
    print('   distinct rhombi: R^E_{k+3} != R_{k+2}-alpha', RE3 != tr(R, mal), ', != R^L_{k+3}-alpha', RE3 != RL3a)
