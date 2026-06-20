import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb, factorial
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5sym=sp.sympify(d['H5sym']); H5f=sp.lambdify((a,b,j),H5sym,'math')
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
# Delta(j) for 1<=j<=9 via clean form (tip=0): j-2s2(j)+2[v2C(a+6,j)+v2C(b+5,j)+v2H5(j)-v2H5(0)]
def Delta_low(A,B,J):
    return J-2*s2(J)+2*(v2C(A+6,J)+v2C(B+5,J)+v2(H5v(A,B,J))-v2(H5v(A,B,0)))

# sanity vs MN-based Delta
def v2val(A,B,J):
    m=(A+B+5)//2; M=Mj((A,B,5),J,m)
    if M==0: return None
    return J+2*v2(comb(m,J)*M)
import random
bad=0
for _ in range(150):
    B=random.randrange(9,30);A=random.randrange(B,B+30)
    if (A+B)%2==0:A+=1
    if A<B: continue
    v0=v2val(A,B,0)
    for J in range(1,min(9,B)+1):
        vj=v2val(A,B,J)
        if vj is None: continue
        if vj-v0 != Delta_low(A,B,J): bad+=1
print("Delta_low vs MN check, bad=",bad)

# min over residues for j=1..9, both branches, plus tie class enumeration mod 4
for J in range(1,10):
    for parity in ['aeven','aodd']:
        ap=0 if parity=='aeven' else 1; bp=1-ap
        mn_=10**9; tieclasses=set()
        for A in range(max(J, 8)+ (0), max(J,8)+200):
            if A%2!=ap: continue
            for B in range(max(J,9), A+1):
                if B%2!=bp: continue
                Dj=Delta_low(A,B,J)
                if Dj<mn_: mn_=Dj
        print(f"  j={J} {parity}: minDelta={mn_}")
