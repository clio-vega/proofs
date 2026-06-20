import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from math import factorial
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5f=sp.lambdify((a,b,j),sp.sympify(d['H5sym']),'math')
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def vfact(r): return r-s2(r)
def vff(x,r):
    s=0
    for t in range(r): s+=v2(x-t)
    return s
def vpsi(A,B,J): return vff(A+2,J-4)+vff(B+1,J-4)+v2(H5v(A,B,J))
def Delta(A,B,J): return J-2*s2(J)+2*(vpsi(A,B,J)-2*vfact(J))

# 1) Among VALID a-even shapes (a even,b odd, a>=b>=10), find min Delta(j) for j=6,7,8,9, wide range
for J in [6,7,8,9]:
    for parity,(ap,bp) in [('aeven',(0,1)),('aodd',(1,0))]:
        mn_=10**9; arg=None
        for B in range(11 if bp else 10, 1200):
            if B%2!=bp: continue
            for A in range(B, 4000):
                if A%2!=ap: continue
                Dj=Delta(A,B,J)
                if Dj<mn_: mn_=Dj; arg=(A,B)
        print(f"j={J} {parity}: min Delta over valid (a>=b>=10) wide = {mn_} at {arg}")
