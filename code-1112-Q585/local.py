"""Local matrices of the KT puzzle model at the Schur (nonequivariant) point.
Labels 0,1,2 (2 := '10').  Pieces from /home/clio/git/puzzles/solve_puzzle.py (Sage's set):
 Delta (bottom, left/, right\\): (0,0,0),(1,0,2),(2,1,0),(1,1,1),(0,2,1)
 Nabla (top, left\\, right/):    see below, derived as 180-degree rotation of Delta.
Diamond = Delta (above) + Nabla (below) sharing the horizontal edge x.
 inputs: UL = Delta.left (/), UR = Delta.right (\\);  outputs: LL = Nabla.left (\\), LR = Nabla.right (/).
"""
from itertools import product
L=(0,1,2)
DELTA={(0,0,0),(1,0,2),(2,1,0),(1,1,1),(0,2,1)}
# Nabla triples from solve_puzzle.py docstring (north-, SE\\, SW/): (0,0,0),(0,10,1),(10,1,0),(1,0,10),(1,1,1)
# In a down triangle the left side is \\ (SW... ) -- take Sage convention (north, SE, SW) and map to (top, right=SE?, left=SW?)
NAB_SAGE={(0,0,0),(0,2,1),(2,1,0),(1,0,2),(1,1,1)}
# down-triangle: top edge horizontal; its two other sides meet at the bottom vertex.
# SE side = right side (going down-left from top-right corner, a '\\'?? ) -- record both readings and test
def diamond(NAB):
    D={}
    for ul,ur,ll,lr in product(L,repeat=4):
        D[(ul,ur,ll,lr)]=sum((x,ul,ur) in DELTA and NAB(x,ll,lr) for x in L)
    return D
readings={'top,left=SW,right=SE':lambda x,ll,lr:(x,lr,ll) in NAB_SAGE,
          'top,left=SE,right=SW':lambda x,ll,lr:(x,ll,lr) in NAB_SAGE}
for name,f in readings.items():
    D=diamond(f)
    nz={k:v for k,v in D.items() if v}
    print(name, 'nonzero entries',len(nz),'max',max(nz.values()))
    for (ul,ur) in product(L,repeat=2):
        outs=[(ll,lr,v) for (a,b,ll,lr),v in nz.items() if (a,b)==(ul,ur)]
        print('  in',(ul,ur),'->',outs)

import sympy as sp
idx=list(product(L,repeat=2))
for name,f in readings.items():
    D=diamond(f)
    M=sp.Matrix(9,9,lambda o,i: D[idx[i]+idx[o]])   # column = input (ul,ur), row = output (ll,lr)
    Psw=sp.Matrix(9,9,lambda o,i: 1 if idx[o]==idx[i][::-1] else 0)  # tensor flip
    print(name,'rank',M.rank(),' (P M)^2 == P M ?',(Psw*M)**2==Psw*M,' (M P)^2==M P?',(M*Psw)**2==M*Psw)
    print('  eigenvalues of P*M:',(Psw*M).eigenvals())
