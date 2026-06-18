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
Q4f = sp.lambdify((a_,b_,j_), Q4, 'math')

def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2; k+=1
    return k
def s2(n): return bin(n).count('1')
def vC(n,k):
    if k<0 or k>n: return 10**9
    return s2(k)+s2(n-k)-s2(n)

def T(av,bv,jv):
    return vC(av+5,jv)+vC(bv+4,jv)+v2(Q4f(av,bv,jv))-v2(Q4f(av,bv,0))

# 1. Check Compensation Lemma  T(j) >= 1 - v2(j)  i.e. T(j)+v2(j) >= 1
print("=== Compensation Lemma T(j)+v2(j) >= 1 ? ===")
minTV=10**9; argm=None
fails=0
for av in range(4,220):
    for bv in range(4,av+1):
        if (av+bv)%2!=0: continue
        for jv in range(1,min(bv,50)+1):
            tv=T(av,bv,jv)+v2(jv)
            if tv<minTV: minTV=tv; argm=(av,bv,jv)
            if tv<1:
                fails+=1
                if fails<8: print("  FAIL",(av,bv,jv),"T+v2=",tv)
print(f"min T(j)+v2(j) = {minTV} at {argm}; fails(<1)={fails}")

# 2. At j=4 specifically (tip-free, P8(4)=0): what is min Delta(4)?
print("\n=== j=4 detail: Delta(4)=4-2*1+2T(4)=2+2T(4) ===")
minT4=10**9; arg4=None
for av in range(4,300):
    for bv in range(4,av+1):
        if (av+bv)%2!=0: continue
        t4=T(av,bv,4)
        if t4<minT4: minT4=t4; arg4=(av,bv)
print(f"min T(4)={minT4} at {arg4}; so min Delta(4)=2+2*minT4={2+2*minT4}")

# 3. j=2 detail
print("\n=== j=2 detail: Delta(2)=2T(2) ===")
minT2=10**9
ties2=0; tot2=0
for av in range(4,200):
    for bv in range(4,av+1):
        if (av+bv)%2!=0: continue
        t2=T(av,bv,2)
        minT2=min(minT2,t2)
        tot2+=1
        if t2==0: ties2+=1
print(f"min T(2)={minT2}; tie(T2=0) count={ties2}/{tot2}")
