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
Q4f=sp.lambdify((a_,b_,j_),Q4,'math'); Hf=sp.lambdify((a_,b_,j_),H,'math')
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n):return bin(n).count('1')
def vC(n,k):
    if k<0 or k>n:return 10**9
    return s2(k)+s2(n-k)-s2(n)
def Delta(av,bv,jv):
    return jv-2*s2(jv)+2*(vC(av+5,jv)+vC(bv+4,jv)+v2(Q4f(av,bv,jv))-v2(Q4f(av,bv,0)))

# 1) Confirm case split & case A bound: when v2(heavy*H(j)) >= v2(P8(j)) [case A],
#    claimed Delta(j) >= j-4-2 s2(j-8) > 0.
print("=== j>=8: verify case-A bound and classify ===")
caseA=caseB=0
fail=[]
worstB=10**9; worstBat=None
for jv in range(8,40):
    for av in range(jv,170):
        for bv in range(jv,av+1):
            if (av+bv)%2: continue
            heavy=v2(av-2)+v2(bv-3)
            vP8=v2(P8.subs(j_,jv))  # constant
            vheavyH=heavy+v2(Hf(av,bv,jv))
            d=Delta(av,bv,jv)
            if vheavyH>=vP8:
                caseA+=1
                bound=jv-4-2*s2(jv-8)
                if d<bound: fail.append(('A',av,bv,jv,d,bound))
            else:
                caseB+=1
                if d<worstB: worstB=d; worstBat=(av,bv,jv)
            if d<=0: fail.append(('VIOL',av,bv,jv,d))
vP8 = v2(int(sp.prod([8-t for t in range(8)])))
print(f"caseA={caseA} caseB={caseB}; caseA-bound failures: {[f for f in fail if f[0]=='A'][:5]}")
print(f"any Delta<=0: {[f for f in fail if f[0]=='VIOL'][:5]}")
print(f"worst Delta in case B = {worstB} at {worstBat}  (still >0)")

# 2) Wide certification of WHOLE theorem via closed form, larger range
print("\n=== Wide closed-form certification: Delta(j)>0 off j=2, Delta(2)>=0 ===")
viol=0
for av in range(4,300):
    for bv in range(4,av+1):
        if (av+bv)%2: continue
        for jv in range(1,bv+1):
            d=Delta(av,bv,jv)
            if jv==2 and d<0: viol+=1
            if jv!=2 and d<=0: viol+=1
print("violations (a<300):",viol)
