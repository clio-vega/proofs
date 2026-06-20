import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from math import factorial
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5sym=sp.sympify(d['H5sym']); H5f=sp.lambdify((a,b,j),H5sym,'math')
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def s2(n): return bin(n).count('1')
def vfact(r): return r-s2(r) if r>=0 else 0
def binom(x,r):
    num=1
    for t in range(r): num*=(x-t)
    return num//factorial(r)
def g(J): return J-2*s2(J)-4*(vfact(J)-vfact(J-4))+6
def wj(J):
    target=2 if J%2==0 else 1
    return (target-g(J)+6)//2
def div_ok(n,w): return (abs(int(n)) % (2**w))==0  # 2^w | n

exc=[10,11,16,17,18,19]
# RIGOROUS: 2^{w_j} | C(a+2,r)*C(b+1,r)*H5(a,b,j) for all a,b of given parity.
# Check over full residue window mod 2^K, for K and K+2, confirm identical (0 failures).
def run(K):
    res={}
    for J in exc:
        r=J-4; w=wj(J); M=2**K
        for parity in ['aeven','aodd']:
            a_par=0 if parity=='aeven' else 1; b_par=1-a_par
            fails=0
            # a,b range over one full period mod M, shifted up so a+2>=r
            aoff=M  # a in [M, 2M)
            for A in range(aoff+a_par, aoff+M, 2):
                Cab=binom(A+2,r)
                for B in range(aoff+b_par, aoff+M, 2):
                    prod=Cab*binom(B+1,r)*H5v(A,B,J)
                    if not div_ok(prod,w): fails+=1
            res[(J,parity)]=fails
    return res
import time
for K in [9,11]:
    t=time.time()
    r=run(K)
    bad=sum(v for v in r.values())
    print(f"K={K}: total failures across all (j,parity) = {bad}   (time {time.time()-t:.1f}s)")
    if bad>0:
        for k,v in r.items():
            if v: print("   FAIL",k,v)
print("\nw_j values:", {J:wj(J) for J in exc})
