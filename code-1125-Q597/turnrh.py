"""PROVE 11-25: where is the re-run's last rhombus R_p+gamma relative to J_o at the turn?  Word model.
R* class 0, A (toward T) takes (0,b) to R_end; R_p+gamma = R_end - alpha.  J_o (toward T') takes (0,b) at R* to R^o_{e+1},
then type t at e+1 to R^o_{e+2}.  Report equalities, and the vector R_p+gamma - R^o_{e+1}.  Also: rem/cre disjoint in all 12 configs."""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597')
import words as W
from wmodel import *
CL0 = frozenset({(0, 1), (0, -1)}); CL1 = {'B': frozenset({(1, -1), (-1, 1)}), 'C': frozenset({(1, 0), (-1, 0)})}
AT = {'B': 'a1', 'C': 'a0'}; OTHER = {'B': 'C', 'C': 'B'}; CF = {T: W.configs(T) for T in 'BC'}
for T in 'BC':
    for k, (Rc, rem, Rn, cre) in CF[T].items(): assert not (rem & cre), 'rem/cre overlap'
print('rem and cre disjoint in all 12 configurations')
for T in 'BC':
    Tp = OTHER[T]
    Rstar = next(v[0] for k, v in CF[T].items() if k[0] == CL0)
    _, _, Rend = step(T, Rstar, 'b')
    Rp = W.canon if False else None
    _, _, R1 = step(T, next(v[0] for k, v in CF[T].items() if k[0] == CL1[T]), AT[T])
    alpha = vsub(W.base(R1), W.base(next(v[0] for k, v in CF[T].items() if k[0] == CL1[T])))
    Rpg = W.trs(Rend, tuple(-x for x in alpha))
    _, creE, Re1 = step(Tp, Rstar, 'b')
    print(T, 'alpha', alpha, 'R_p+gamma == R^o_{e+1}:', Rpg == Re1, ' R_p+gamma in cre(e):', Rpg in creE,
          ' vector R_p+gamma - R^o_{e+1}:', vsub(W.base(Rpg), W.base(Re1)))
    for t in ('a0', 'a1', 'b'):
        _, creE1, Re2 = step(Tp, Re1, t)
        print('   e+1 type', t, 'straight' if t == AT[Tp] else '', 'R_p+gamma == R^o_{e+2}:', Rpg == Re2)
