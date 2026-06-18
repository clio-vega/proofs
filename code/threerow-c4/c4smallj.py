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
def T(av,bv,jv):
    return vC(av+5,jv)+vC(bv+4,jv)+v2(Q4f(av,bv,jv))-v2(Q4f(av,bv,0))
def Delta(av,bv,jv):
    return jv-2*s2(jv)+2*T(av,bv,jv)

# small j, with VALID constraint a>=b>=4
print("min Delta(j) over a>=b>=4 (a,b same parity), j<=b:")
for jj in range(1,8):
    mn=10**9;arg=None
    for av in range(4,400):
        for bv in range(4,av+1):
            if bv<jj: continue
            if (av+bv)%2: continue
            d=Delta(av,bv,jj)
            if d<mn:mn=d;arg=(av,bv)
    print(f"  j={jj}: minDelta={mn} at {arg}")

# Double-check: is there ANY interior j (1<=j<=b) off {0,2} with Delta<=0, a>=b>=4?
print("\nGlobal interior check a<=200:")
viol=[]
for av in range(4,200):
    for bv in range(4,av+1):
        if (av+bv)%2: continue
        for jj in range(1,bv+1):
            d=Delta(av,bv,jj)
            if d<=0 and jj!=2:
                viol.append((av,bv,jj,d))
            if d<0 and jj==2:
                viol.append(('NEG2',av,bv,jj,d))
print("violations (Delta<=0 off j=2, or Delta(2)<0):",len(viol),viol[:10])
