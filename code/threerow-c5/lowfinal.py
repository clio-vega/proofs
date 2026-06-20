import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from math import factorial
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5sym=sp.sympify(d['H5sym']); H5f=sp.lambdify((a,b,j),H5sym,'math')
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def s2(n): return bin(n).count('1')
def vfact(r): return r-s2(r)
def binom(x,r):
    num=1
    for t in range(r): num*=(x-t)
    return num//factorial(r)
def div_ok(n,w): return (abs(int(n))%(2**w))==0

# required Delta floor per (j,parity)
floor_even={4:2,5:1,6:2,7:1,8:2,9:1}   # a-even: need Delta>0
floor_odd ={4:-2,5:-3,6:-2,7:-1,8:-2,9:-1}  # a-odd: no new tie (except j=5)
def kappa(J,f): 
    return 2*vfact(J)+(f-J+2*s2(J))//2

# Psi_j = (a+2)_{j-4} (b+1)_{j-4} H5(j). Verify 2^{kappa}|Psi_j over full residues mod 2^K, K and K+2.
def run(K):
    out={}
    for J in range(4,10):
        r=J-4
        for parity,floors in [('aeven',floor_even),('aodd',floor_odd)]:
            ap=0 if parity=='aeven' else 1; bp=1-ap
            kap=kappa(J,floors[J])
            M=2**K; off=M
            fails=0
            for A in range(off+ap, off+M, 2):
                Ca=binom(A+2,r)
                for B in range(off+bp, off+M, 2):
                    if not div_ok(Ca*binom(B+1,r)*H5v(A,B,J), kap): fails+=1
            out[(J,parity,kap)]=fails
    return out
import time
for K in [8,10]:
    t=time.time(); o=run(K); bad=sum(o.values())
    print(f"K={K}: total failures(j=4..9,both parities)={bad}  (kappa used: { {k[0:2]:k[2] for k in o} })  time={time.time()-t:.1f}s")
    if bad: 
        for k,v in o.items():
            if v: print("  FAIL",k,v)

# tie/offset class determination: check Delta(2)(aeven) and Delta(5)(aodd) depend only on (a,b) mod 4
# by verifying mod 8 gives consistent value per mod-4 class
import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def v2C(n,k):
    from math import comb
    if k<0 or k>n: return 10**9
    return v2(comb(n,k))
def Delta_low(A,B,J):
    return J-2*s2(J)+2*(v2C(A+6,J)+v2C(B+5,J)+v2(H5v(A,B,J))-v2(H5v(A,B,0)))
print("\nTie determined-mod-4 check:")
for (J,ap,bp,lbl) in [(2,0,1,'j=2 aeven'),(5,1,0,'j=5 aodd')]:
    from collections import defaultdict
    byclass=defaultdict(set)
    for A in range(16,16+400):
        if A%2!=ap: continue
        for B in range(16,A+1):
            if B%2!=bp: continue
            byclass[(A%4,B%4)].add(Delta_low(A,B,J))
    # report classes whose Delta hits the tie value
    tieval = 0 if J==2 else -3
    tieclasses=[k for k,v in byclass.items() if tieval in v]
    nontie_min={k:min(v) for k,v in byclass.items()}
    print(f"  {lbl}: classes achieving tie value {tieval}: {tieclasses}; min per class: {dict(nontie_min)}")
