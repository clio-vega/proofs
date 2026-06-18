import sympy as sp
a_,b_,j_ = sp.symbols('a b j', integer=True)
H = (a_**3*b_**3 + 9*a_**3*b_**2 + 26*a_**3*b_ + 24*a_**3 + 12*a_**2*b_**3
 - 4*a_**2*b_**2*j_**2 - 8*a_**2*b_**2*j_ + 108*a_**2*b_**2 - 12*a_**2*b_*j_**2
 - 48*a_**2*b_*j_ + 312*a_**2*b_ - 8*a_**2*j_**2 - 64*a_**2*j_ + 288*a_**2
 + 47*a_*b_**3 - 20*a_*b_**2*j_**2 - 64*a_*b_**2*j_ + 423*a_*b_**2 + 6*a_*b_*j_**4
 - 12*a_*b_*j_**3 - 30*a_*b_*j_**2 - 384*a_*b_*j_ + 1222*a_*b_ + 24*a_*j_**3
 - 40*a_*j_**2 - 488*a_*j_ + 1128*a_ + 60*b_**3 - 24*b_**2*j_**2 - 120*b_**2*j_
 + 540*b_**2 + 6*b_*j_**4 + 12*b_*j_**3 - 42*b_*j_**2 - 696*b_*j_ + 1560*b_
 - 4*j_**6 + 48*j_**5 - 244*j_**4 + 648*j_**3 - 664*j_**2 - 648*j_ + 1440)
P8 = sp.prod([j_-t for t in range(8)])
Q4 = sp.expand((a_-2)*(b_-3)*H + P8)
Q4f = sp.lambdify((a_,b_,j_),Q4,'math')
Hf  = sp.lambdify((a_,b_,j_),H,'math')

def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def vC(n,k):
    if k<0 or k>n: return 10**9
    return s2(k)+s2(n-k)-s2(n)
def vfall(top,r):
    return 0 if r<=0 else sum(v2(top-t) for t in range(r))
def vfact(r): return sum(v2(t) for t in range(1,r+1))

def Tpeel(av,bv,jv):
    return (vfall(av+2,jv-3)+vfall(bv+1,jv-3)-2*vfact(jv)+v2(Q4f(av,bv,jv))-v2(av-2)-v2(bv-3))
def Delta(av,bv,jv):
    return jv-2*s2(jv)+2*Tpeel(av,bv,jv)

# Lemma checks j in 4,6 : the sum bounds
print("=== j=4: min of v2(a+2)+v2(b+1)+v2H(4)  (need >=6) ===")
mn=10**9; arg=None
for av in range(4,400):
    for bv in range(4,av+1):
        if (av+bv)%2: continue
        s=v2(av+2)+v2(bv+1)+v2(Hf(av,bv,4))
        if s<mn: mn=s; arg=(av,bv)
print("min=",mn,"at",arg)

print("=== j=6: min of v2((a+2)(a+1)a)+v2((b+1)b(b-1))+v2H(6) (need >=8) ===")
mn=10**9; arg=None
for av in range(6,400):
    for bv in range(6,av+1):
        if (av+bv)%2: continue
        s=vfall(av+2,3)+vfall(bv+1,3)+v2(Hf(av,bv,6))
        if s<mn: mn=s; arg=(av,bv)
print("min=",mn,"at",arg)

# odd j=3,5,7: Delta
print("\n=== odd small j: min Delta ===")
for jj in [3,5,7]:
    mn=10**9; arg=None
    for av in range(jj,400):
        for bv in range(jj,av+1):
            if (av+bv)%2: continue
            d=Delta(av,bv,jj)
            if d<mn: mn=d; arg=(av,bv)
    print(f"  j={jj}: min Delta={mn} at {arg}")

# j>=8 regime: min Delta per j, and the minimizing shapes
print("\n=== j>=8: min Delta(j) (interior b>=j) over wide a,b ===")
for jj in range(8,21):
    mn=10**9; arg=None
    for av in range(jj,260):
        for bv in range(jj,av+1):
            if (av+bv)%2: continue
            d=Delta(av,bv,jj)
            if d<mn: mn=d; arg=(av,bv)
    print(f"  j={jj}: minDelta={mn} at {arg}")
