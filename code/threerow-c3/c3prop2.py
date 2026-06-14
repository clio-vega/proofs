import sympy as sp
from math import comb
from mn import Mj

a,b,j = sp.symbols('a b j', integer=True)
Q3 = ((a-1)*(a+3)*(a+4)*(b-2)*(b+2)*(b+3)
      + (-j**6+15*j**5 + j**4*(3*a*b-6*a-3*b-79) + j**3*(-12*a*b+24*a+12*b+201)
         + j**2*(-3*a**2*b**2+3*a**2*b+6*a**2-3*a*b**2+24*a*b-36*a+6*b**2-27*b-244)
         + j*(-3*a**2*b**2-3*a**2*b+18*a**2-9*a*b**2-15*a*b+66*a+12*b**2+18*b+36)))

def v2(n):
    n=int(n)
    if n==0: return 10**9
    c=0
    while n%2==0: n//=2; c+=1
    return c
def s2(n):
    return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N: return 10**9
    return s2(k)+s2(N-k)-s2(N)

def Delta_pred(av,bv,jv):
    q0=int(Q3.subs({a:av,b:bv,j:0}))
    qj=int(Q3.subs({a:av,b:bv,j:jv}))
    return jv-2*s2(jv)+2*v2C(av+4,jv)+2*v2C(bv+3,jv)+2*(v2(qj)-v2(q0))

# direct val(j)-val(0)
def val(av,bv,jv,mv):
    M=Mj((av,bv,3),jv,mv)
    if M==0: return None
    return jv+2*(v2(comb(mv,jv))+v2(M))

bad=0; tested=0
for mv in range(3,22):
    n=2*mv
    for av in range(1,n+1):
        for bv in range(3,av+1):
            if n-av-bv!=3: continue
            v0=val(av,bv,0,mv)
            for jv in range(1,bv+1):
                vj=val(av,bv,jv,mv)
                if vj is None: continue
                tested+=1
                pred=Delta_pred(av,bv,jv)
                if vj-v0!=pred:
                    bad+=1
                    if bad<20: print("MISMATCH",(av,bv,3),"m",mv,"j",jv,"direct",vj-v0,"pred",pred)
print(f"Prop-2 check: tested={tested} mismatches={bad}")
