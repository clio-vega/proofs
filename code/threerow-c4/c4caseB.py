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
Hf=sp.lambdify((a_,b_,j_),H,'math')
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n):return bin(n).count('1')
def vfall(top,r):return 0 if r<=0 else sum(v2(top-t) for t in range(r))
def vfact(r):return sum(v2(t) for t in range(1,r+1))
def Dhat(av,bv,jv):  # heavy-free quantity
    return jv-2*s2(jv)+2*(vfall(av+2,jv-3)+vfall(bv+1,jv-3)-2*vfact(jv)+v2(Hf(av,bv,jv)))

# Is Dhat(j)>0 for ALL a==b mod 2 (b>=j), j=8..12 ?
print("min Dhat(j) over ALL shapes a>=b>=j, a==b mod2:")
for jj in range(8,13):
    mn=10**9;arg=None
    for av in range(jj,200):
        for bv in range(jj,av+1):
            if (av+bv)%2:continue
            d=Dhat(av,bv,jj)
            if d<mn:mn=d;arg=(av,bv)
    print(f"  j={jj}: min Dhat={mn} at {arg}")
