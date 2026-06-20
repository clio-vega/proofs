import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb
import sympy as sp, pickle, random

a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5sym=sp.sympify(d['H5sym']); Q5sym=sp.sympify(d['Q5sym'])
H5f=sp.lambdify((a,b,j),H5sym,'math'); Q5f=sp.lambdify((a,b,j),Q5sym,'math')
def Q5v(A,B,J): return int(round(Q5f(A,B,J)))
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def v2(n):
    n=abs(int(n)); 
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def vfact(r): return r-s2(r) if r>=0 else 0
def vff(x,r):  # v2 of falling factorial (x)_r = x(x-1)..(x-r+1)
    s=0
    for t in range(r): s+=v2(x-t)
    return s
def v2C(n,k):
    if k<0 or k>n: return 10**9
    return v2(comb(n,k))
def Tfun(A,B,J): return v2C(A+6,J)+v2C(B+5,J)+v2(Q5v(A,B,J))-v2(Q5v(A,B,0))

# Verify peeling identity for j>=4:
#  T(j) = v2((a+2)_{j-4}) + v2((b+1)_{j-4}) - 2 vfact(j) + v2 Q5(j) - v2(a-3) - v2(b-4)
bad=0;tested=0
for _ in range(300):
    B=random.randrange(6,40);A=random.randrange(B,B+50)
    if (A+B)%2==0:A+=1
    if A<B: continue
    for J in range(4,B+1):
        lhs=Tfun(A,B,J)
        rhs=vff(A+2,J-4)+vff(B+1,J-4)-2*vfact(J)+v2(Q5v(A,B,J))-v2(A-3)-v2(B-4)
        tested+=1
        if lhs!=rhs:
            bad+=1
            if bad<5: print("PEEL MISMATCH",A,B,J,lhs,rhs)
print(f"Peeling identity (j>=4): tested={tested} mismatches={bad}")

# g(j) for deep indices (Case B floor): g(j)=j-2s2(j)-4*v2(j(j-1)(j-2)(j-3))+6
def g(J): return J-2*s2(J)-4*(vfact(J)-vfact(J-4))+6
print("\ng(j) for j=10..40:")
exc=[]
for J in range(10,41):
    gv=g(J)
    if gv<=0: exc.append(J)
    print(f"  j={J}: g={gv}{'  <=0!' if gv<=0 else ''}")
print("exceptions (g<=0) in 10..40:",exc)
# analytic: for large j, v2(4 consec)<=2+floor(log2 j); s2<=floor(log2 j)+1
import math
print("\ncheck g(j)>0 for all j>=? via bound j-6log2(j)... :")
for J in [40,48,56,64,80,100,128]:
    print(f"  j={J}: g={g(J)}")
