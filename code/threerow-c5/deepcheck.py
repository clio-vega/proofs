import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from math import comb
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5sym=sp.sympify(d['H5sym']); Q5sym=sp.sympify(d['Q5sym'])
H5f=sp.lambdify((a,b,j),H5sym,'math'); Q5f=sp.lambdify((a,b,j),Q5sym,'math')
def Q5v(A,B,J): return int(round(Q5f(A,B,J)))
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def v2C(n,k):
    if k<0 or k>n: return 10**9
    return v2(comb(n,k))
def Delta(A,B,J):  # Prop2, requires J<=b
    return J-2*s2(J)+2*(v2C(A+6,J)+v2C(B+5,J)+v2(Q5v(A,B,J))-v2(Q5v(A,B,0)))

# For each exceptional j, sweep a,b widely (b>=j, box-interior b>=10, a>=b, a+b odd)
exc=[10,11,16,17,18,19]
AMAX=2200
for J in exc:
    mins={'aeven':10**9,'aodd':10**9}
    argmin={'aeven':None,'aodd':None}
    for B in range(max(10,J), 600):
        for A in range(B, AMAX):
            if (A+B)%2==0: continue
            if B>A: continue
            Dj=Delta(A,B,J)
            par='aeven' if A%2==0 else 'aodd'
            if Dj<mins[par]:
                mins[par]=Dj; argmin[par]=(A,B)
    print(f"j={J}: min Delta a-even={mins['aeven']} at {argmin['aeven']}; a-odd={mins['aodd']} at {argmin['aodd']}")
