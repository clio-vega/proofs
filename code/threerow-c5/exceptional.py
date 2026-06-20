import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from math import comb, factorial
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5sym=sp.sympify(d['H5sym']); H5f=sp.lambdify((a,b,j),H5sym,'math')
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def v2cap(n,K):
    n=abs(int(n))
    if n==0: return K
    k=0
    while n%2==0 and k<K: n//=2;k+=1
    return k
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
exc=[10,11,16,17,18,19]
allok=True
for J in exc:
    r=J-4; w=wj(J); K=w+5; M=2**K; base=128*M
    for parity in ['aeven','aodd']:
        a_par=0 if parity=='aeven' else 1; b_par=1-a_par
        # collect a-residues with va<w
        ares=[]
        for ra in range(M):
            if ra%2!=a_par: continue
            A=ra+base
            va=v2cap(binom(A+2,r),K)
            if va<w: ares.append((ra,A,va))
        bres=[]
        for rb in range(M):
            if rb%2!=b_par: continue
            B=rb+base
            vb=v2cap(binom(B+1,r),K)
            if vb<w: bres.append((rb,B,vb))
        minW=999; arg=None
        for ra,A,va in ares:
            for rb,B,vb in bres:
                if va+vb>=w: 
                    W=va+vb
                else:
                    vh=v2cap(H5v(A,B,J),K)
                    W=va+vb+vh
                if W<minW: minW=W; arg=(ra,rb,va,vb,W-va-vb if va+vb<w else 0)
        ok = minW>=w
        allok &= ok
        print(f"j={J} ({parity}): need W>={w}; |ares|={len(ares)} |bres|={len(bres)}; minW={minW} arg(ra,rb,va,vb,vh)={arg}  {'OK' if ok else 'FAIL'}")
print("\nALL EXCEPTIONAL CHECKS:", "PASS" if allok else "FAIL")
