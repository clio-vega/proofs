import sympy as sp
a,b,j=sp.symbols('a b j')
E=(a*b+a+2*b)*sp.binomial(j,2)-6*sp.binomial(j,3)-12*sp.binomial(j,4)
Phi=a*b+a+2*b-(j-1)*(j-2)
print("E - C(j,2)*Phi =", sp.simplify(sp.expand(E-sp.binomial(j,2)*Phi)))

from math import comb
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N:return 10**9
    return s2(k)+s2(N-k)-s2(N)

# Pi2 bound: 1+v2C(a+4,j)+v2C(b+3,j)+v2C(j,2) >= a3p+b2p+g(j)  [a3p=v2(a+3),b2p=v2(b+2)]
# = 2x: 2+2v2C(a+4,j)+2v2C(b+3,j)+2v2C(j,2) >= 2a3p+2b2p+2s2(j)-(j+3)
# Candidate sublemmas (drop terms / single absorption):
c=dict()
for name in ["Pi2-full","i:vCb>=b2p+g (drop vC(a+2,j-2))","ii:vCa>=a3p+g (drop vC(b+1,j-2))",
             "iii a-side:vC(a+2,j-2)+vCb>=b2p+g","iv b-side:vCa+vC(b+1,j-2)>=a3p+g"]:
    c[name]=0
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            a3p,b2p=v2(av+3),v2(bv+2)
            for jv in range(4,bv+1):
                G2=2*s2(jv)-(jv+3)  # =2g(j)
                vCa,vCb=v2C(av+4,jv),v2C(bv+3,jv)
                if 2+2*vCa+2*vCb+2*v2C(jv,2) < 2*a3p+2*b2p+G2: c["Pi2-full"]+=1
                if 2*vCb < 2*b2p+G2: c["i:vCb>=b2p+g (drop vC(a+2,j-2))"]+=1
                if 2*vCa < 2*a3p+G2: c["ii:vCa>=a3p+g (drop vC(b+1,j-2))"]+=1
                # a-side: absorb C(j,2) into a: v2 Pi2 = a3p + v2C(a+2,j-2)+vCb ; need >= a3p+b2p+g
                if 2*(v2C(av+2,jv-2)+vCb) < 2*b2p+G2: c["iii a-side:vC(a+2,j-2)+vCb>=b2p+g"]+=1
                if 2*(vCa+v2C(bv+1,jv-2)) < 2*a3p+G2: c["iv b-side:vCa+vC(b+1,j-2)>=a3p+g"]+=1
print("E=C(j,2)Phi confirmed above.")
for k_,v_ in c.items(): print(f"{v_:6d}  {k_}")
