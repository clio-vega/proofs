import sys; sys.path.insert(0,'/home/clio/projects/scratch/20261021c2')
from g3 import *

print('=== POSITIVE CONTROL, named in the pre-registration: N=3, m=2, I={2,3}, type B, eps=1')
print('    (r=|I|=2>=1 so NOT the r=0 vacuum; I != [N] so the T-blocks are non-empty;')
print('     all four blocks of Inv u T are populated.)')
N=3; I=(2,3); m=2; w=wvars(N)
for lab,f in diag_factors(N,I,w,'B'):
    print(f'    {lab:16s}  {str(f):26s}  ord_S = {ord_S(f,w,m,1)}')
A=Atilde(N,w,I,I,X='B')
print(f'    SUM over factors                           ord_S = {sum(ord_S(f,w,m,1) for _,f in diag_factors(N,I,w,"B"))}')
print(f'    ord_S of Atilde^B_(I,I) from SYMMETRISATION ord_S = {ord_S(A,w,m,1)}')
print(f'    predicted [m in I]                               = {1 if m in I else 0}')
print()
print('=== the char-sensitivity of "order EXACTLY one", as an identity ===')
v=sp.symbols("v")
for h in (1,2,3):
    co=sp.simplify(sp.cancel((1-v**(-h))/(v-1)))
    print(f'    (1 - v^-{h})/(v-1) = {sp.factor(co)};  value at v=1 is {sp.simplify(co.subs(v,1))}')
print('    -> the cofactor at v=1 is exactly  h * (-1) ... i.e. |h|; nonzero iff h != 0 in the base field.')
print('       type B: h=1, char-free.   type C: h=2, needs char != 2.')
