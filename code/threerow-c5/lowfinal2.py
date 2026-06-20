import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5f=sp.lambdify((a,b,j),sp.sympify(d['H5sym']),'math')
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def v2cap(n,K):
    n=abs(int(n))
    if n==0: return K
    k=0
    while n%2==0 and k<K: n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def vfact(r): return r-s2(r)
def vffcap(x,r,K):
    s=0
    for t in range(r): s+=v2cap(x-t,K)
    return s
floor_even={4:2,5:1,6:2,7:1,8:2,9:1}
floor_odd ={4:-2,5:-3,6:-2,7:-1,8:-2,9:-1}
def kappa(J,f): return 2*vfact(J)+(f-J+2*s2(J))//2
def run(K):
    out={}
    for J in range(4,10):
        r=J-4
        for parity,floors in [('aeven',floor_even),('aodd',floor_odd)]:
            ap=0 if parity=='aeven' else 1; bp=1-ap
            kap=kappa(J,floors[J]); M=2**K; off=M; fails=0; mnv=999
            for A in range(off+ap, off+M, 2):
                vA=vffcap(A+2,r,K)
                for B in range(off+bp, off+M, 2):
                    vp=vA+vffcap(B+1,r,K)+v2cap(H5v(A,B,J),K)
                    if vp<mnv: mnv=vp
                    if vp<kap: fails+=1
            out[(J,parity)]=(kap,mnv,fails)
    return out
import time
for K in [8,10]:
    t=time.time(); o=run(K); bad=sum(v[2] for v in o.values())
    print(f"K={K}: total failures={bad}  time={time.time()-t:.1f}s")
    for k,v in o.items(): 
        if v[2]: print("  FAIL",k,"kappa,minv,fails=",v)
    print("  (kappa,minv) per (j,par):", {k:(v[0],v[1]) for k,v in o.items()})
