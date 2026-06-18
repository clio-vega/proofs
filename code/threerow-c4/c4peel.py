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
def vfall(top,r):  # v2 of (top)_r = top(top-1)...(top-r+1)
    if r<=0: return 0
    return sum(v2(top-t) for t in range(r))
def vfact(r): return sum(v2(t) for t in range(1,r+1))

def T(av,bv,jv):
    return vC(av+5,jv)+vC(bv+4,jv)+v2(Q4f(av,bv,jv))-v2(Q4f(av,bv,0))

def T_peel(av,bv,jv):
    # claimed = vfall(a+2,j-3)+vfall(b+1,j-3) - 2 vfact(j) + v2Q4(j) - v2(a-2) - v2(b-3)
    return (vfall(av+2,jv-3)+vfall(bv+1,jv-3) - 2*vfact(jv)
            + v2(Q4f(av,bv,jv)) - v2(av-2) - v2(bv-3))

# Verify peeling identity for j>=3
print("=== Verify peeling identity (j>=3) ===")
bad=0
for av in range(4,120):
    for bv in range(4,av+1):
        if (av+bv)%2!=0: continue
        for jv in range(3,min(bv,30)+1):
            if T(av,bv,jv)!=T_peel(av,bv,jv):
                bad+=1
                if bad<6: print("  MISMATCH",(av,bv,jv),T(av,bv,jv),T_peel(av,bv,jv))
print("mismatches:",bad)

# H(j) factored for small j
print("\n=== H(j) factored ===")
for jj in [0,1,2,3,4,5,6,7]:
    print(f"H({jj}) =", sp.factor(sp.expand(H.subs(j_,jj))))
