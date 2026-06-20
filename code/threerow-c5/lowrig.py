import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from math import comb, factorial
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5sym=sp.sympify(d['H5sym']); H5f=sp.lambdify((a,b,j),H5sym,'math')
hcoeffs=[sp.sympify(s) for s in d['hcoeffs']]
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def vfact(r): return r-s2(r)
def v2C(n,k):
    if k<0 or k>n: return 10**9
    return v2(comb(n,k))
def binom(x,r):
    num=1
    for t in range(r): num*=(x-t)
    return num//factorial(r)
def Delta_low(A,B,J):
    return J-2*s2(J)+2*(v2C(A+6,J)+v2C(B+5,J)+v2(H5v(A,B,J))-v2(H5v(A,B,0)))

# (a) constancy of Delta(1),(2),(3) on each branch
from collections import Counter
for parity in ['aeven','aodd']:
    ap=0 if parity=='aeven' else 1; bp=1-ap
    cnt={1:Counter(),2:Counter(),3:Counter()}
    for A in range(8,260):
        if A%2!=ap: continue
        for B in range(9,A+1):
            if B%2!=bp: continue
            for J in [1,2,3]:
                cnt[J][Delta_low(A,B,J)]+=1
    print(f"{parity}: Delta(1) dist={dict(cnt[1])}, Delta(2)={dict(cnt[2])}, Delta(3)={dict(cnt[3])}")

# (b) tie class for j=2 (a-even) and j=5 (a-odd): which (a%4,b%4) give the minimal Delta
print("\n--- j=2 a-even: Delta(2) by (a%4,b%4) ---")
seen={}
for A in range(8,400):
    if A%2!=0: continue
    for B in range(9,A+1):
        if B%2!=1: continue
        key=(A%4,B%4); seen.setdefault(key,set()).add(Delta_low(A,B,2))
for k in sorted(seen): print("  ",k, sorted(seen[k]))
print("\n--- j=5 a-odd: Delta(5) by (a%4,b%4) ---")
seen={}
for A in range(8,400):
    if A%2!=1: continue
    for B in range(9,A+1):
        if B%2!=0: continue
        key=(A%4,B%4); seen.setdefault(key,set()).add(Delta_low(A,B,5))
for k in sorted(seen): print("  ",k, sorted(seen[k]))

# (c) min v2(Psi_j) for j=4..9 per parity, Psi_j=(a+2)_{j-4}(b+1)_{j-4}H5(j)
print("\n--- j=4..9: min v2(Psi_j) and implied Delta floor ---")
def vff(x,r):
    s=0
    for t in range(r): s+=v2(x-t)
    return s
for J in range(4,10):
    r=J-4
    for parity in ['aeven','aodd']:
        ap=0 if parity=='aeven' else 1; bp=1-ap
        mn_=10**9
        for A in range(max(8,J),300):
            if A%2!=ap: continue
            for B in range(max(9,J),A+1):
                if B%2!=bp: continue
                vpsi=vff(A+2,r)+vff(B+1,r)+v2(H5v(A,B,J))
                if vpsi<mn_: mn_=vpsi
        Dfloor=J-2*s2(J)+2*(mn_-2*vfact(J))
        print(f"  j={J} {parity}: min v2(Psi)={mn_}, implied min Delta={Dfloor}")
