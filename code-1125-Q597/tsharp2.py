"""(T#) in the word model: for T, linear functionals that are constant on nest-T frame cells:
   C: frame u=(1,0;0,0), w=(0,0;-1,0): l0 = b0, l1 = b1 (interior side of Delta_k: l >= c_k; frame cells: l = c_k)
   B: frame u=(-1,1;0,0), w=(0,0;1,-1): l0 = -(a0+b0), l1 = -(a1+b1) (interior side l >= c_k)
 Witnesses: the non-rh tiles of A's (0,b) hexagon at R* (M0 wake tiles in a record).  Report min l over their vertices
 minus l(R_p+gamma) (constant on the rhombus).  Need: < 0 for some k  => R_p+gamma is not a frame cell.
 Lesion: the same test for R_{p+m+1} = R_p+gamma+alpha (A's rhombus, which IS a packed cell in 180 C records): must NOT certify."""
from wmodel import *
L = {'C': (lambda P: P[1], lambda P: P[3]), 'B': (lambda P: -(P[0] + P[1]), lambda P: -(P[2] + P[3]))}
for T in 'BC':
    R_p = canon(T); aT = AT[T]
    _, _, Rp1 = step(T, R_p, aT); alpha = vsub(W.base(Rp1), W.base(R_p))
    Rstar = next(v[0] for k, v in CF[T].items() if k[0] == CL0)
    remE, _, Rend = step(T, Rstar, 'b')
    Rpg = W.trs(Rend, tuple(-x for x in alpha))
    for name, R in (('R_p+gamma', Rpg), ('LESION R_{p+m+1}', Rend), ('R*', Rstar)):
        out = []
        for k in (0, 1):
            vals = {L[T][k](P) for P in R[1]}
            const = len(vals) == 1
            mw = min(L[T][k](P) for t in remE for P in t[1])
            out.append(('k', k, 'l const on R', const, 'min_w l - l(R)', mw - min(vals)))
        print(T, name, out)
