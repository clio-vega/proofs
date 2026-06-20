import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
hcoeffs=[sp.sympify(s) for s in d['hcoeffs']]
def binom_poly(jj,k):
    r=sp.Integer(1)
    for t in range(k): r*=(jj-t)
    return r/sp.factorial(k)
H5=lambda J: sp.expand(sum(hcoeffs[k]*binom_poly(J,k) for k in range(9)))
# R-polynomials
R1=sp.expand((a+6)*(b+5)-20)
H1=H5(1); H0=H5(0); H2=H5(2); H3=H5(3); H5_5=H5(5)
# verify H5(1)=(a+3)(a+4)(a+5)(b+2)(b+3)(b+4)*R1
chk1=sp.expand(H1-(a+3)*(a+4)*(a+5)*(b+2)*(b+3)*(b+4)*R1)
print("H5(1) = prefac*R1 :", chk1==0, " R1=(a+6)(b+5)-20")
R2=sp.symbols('R2'); R2poly=sp.expand(sp.simplify(H2/((a+3)*(a+4)*(b+2)*(b+3))))
print("H5(2)/((a+3)(a+4)(b+2)(b+3)) integer poly R2 =", sp.factor(R2poly))
R3poly=sp.expand(sp.simplify(H3/((a+3)*(b+2))))
print("H5(3)/((a+3)(b+2)) = R3 (deg) ok:", sp.simplify(H3-(a+3)*(b+2)*R3poly)==0)

# Numeric verification of reductions: Delta(1)=-1+2 v2(R1); Delta(2)=2 v2(R2)-4; Delta(3)=2 v2(R3)-5
def v2(n):
    n=abs(int(n));
    if n==0: return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def v2C(n,k):
    if k<0 or k>n: return 10**9
    return v2(comb(n,k))
def H5v(A,B,J): return int(H5(J).subs({a:A,b:B}))
def Delta(A,B,J):
    return J-2*s2(J)+2*(v2C(A+6,J)+v2C(B+5,J)+v2(H5v(A,B,J))-v2(H5v(A,B,0)))
R1f=sp.lambdify((a,b),R1,'math'); R2f=sp.lambdify((a,b),R2poly,'math'); R3f=sp.lambdify((a,b),R3poly,'math')
import random
bad=0
for _ in range(300):
    B=random.randrange(9,40);A=random.randrange(B,B+40)
    if (A+B)%2==0:A+=1
    if A<B: continue
    if Delta(A,B,1)!=-1+2*v2(round(R1f(A,B))): bad+=1
    if Delta(A,B,2)!=2*v2(round(R2f(A,B)))-4: bad+=1
    if Delta(A,B,3)!=2*v2(round(R3f(A,B)))-5: bad+=1
print("reductions Delta(1,2,3) numeric check, bad=",bad)

# Residue checks for R-valuation floors (stable two moduli), both branches:
#  a-even: v2(R1)>=2, v2(R2)>=2, v2(R3)>=4 ; a-odd: v2(R1)=0, v2(R2)=1, v2(R3)=1 (exact -> check <= and >=)
def floorcheck(Rf, par_a, par_b, need, K):
    M=2**K; off=M; fails=0
    for A in range(off+par_a, off+M,2):
        for B in range(off+par_b, off+M,2):
            v=v2(round(Rf(A,B)))
            if v<need: fails+=1
    return fails
for K in [6,8]:
    fe1=floorcheck(R1f,0,1,2,K); fe2=floorcheck(R2f,0,1,2,K); fe3=floorcheck(R3f,0,1,4,K)
    print(f"K={K} a-even: v2(R1)>=2 fails={fe1}, v2(R2)>=2 fails={fe2}, v2(R3)>=4 fails={fe3}")
    # a-odd exact: v2(R1)=0 (odd), v2(R2)=1, v2(R3)=1
    fo1=floorcheck(R1f,1,0,1,K)  # >=1 should FAIL fully? we want =0 i.e. odd -> v2>=1 has many; check max instead
    # check v2 exactly
    def exact(Rf,pa,pb,val,K):
        M=2**K; off=M; bad=0
        for A in range(off+pa,off+M,2):
            for B in range(off+pb,off+M,2):
                if v2(round(Rf(A,B)))!=val: bad+=1
        return bad
    print(f"   a-odd exact: v2(R1)=0 bad={exact(R1f,1,0,0,K)}, v2(R2)=1 bad={exact(R2f,1,0,1,K)}, v2(R3)=1 bad={exact(R3f,1,0,1,K)}")

# Tie classes determined mod 4 (stable mod 16): j=2 aeven Delta=0 <=> (0,1); j=5 aodd Delta=-3 <=> (3,0)
def tie_by_class(J,pa,pb,mod):
    from collections import defaultdict
    cls=defaultdict(set)
    base=64
    for A in range(base, base+ mod*8):
        if A%2!=pa: continue
        for B in range(base, base+mod*8):
            if B%2!=pb: continue
            cls[(A%mod,B%mod)].add(Delta(A,B,J))
    return cls
for mod in [4,8]:
    c2=tie_by_class(2,0,1,mod); c5=tie_by_class(5,1,0,mod)
    t2=[k for k,v in c2.items() if 0 in v]; t5=[k for k,v in c5.items() if -3 in v]
    # also check each class is single-valued at the min? print min per class for the tie classes
    print(f"mod {mod}: j=2 aeven Delta=0 classes={t2}; j=5 aodd Delta=-3 classes={t5}")
