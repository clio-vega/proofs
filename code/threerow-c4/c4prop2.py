import sys
sys.path.insert(0,'/home/clio/projects/proofs/code/threerow-c3')
from mn import Mj
from math import comb
import sympy as sp

a_,b_,j_ = sp.symbols('a b j', integer=True)
Q4 = (a_-2)*(b_-3)*0  # placeholder; build via decomposition
# Use decomposition Q4 = (a-2)(b-3)H + P8 with H from c4verify output
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

def v2(n):
    n=abs(n)
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2; k+=1
    return k

def Qv(av,bv,jv): return int(Q4.subs({a_:av,b_:bv,j_:jv}))
def s2(n):
    return bin(n).count('1')

def vC(n,k):
    # v2 C(n,k) via Kummer = number of carries; use s2
    if k<0 or k>n: return 10**9
    return (s2(k)+s2(n-k)-s2(n))

def Delta_closed(av,bv,jv):
    # Prop 2 (c=4): Delta(j)= j - 2 s2(j) + 2 vC(a+5,j)+2 vC(b+4,j) + 2[v2Q(j)-v2Q(0)]
    return (jv - 2*s2(jv) + 2*vC(av+5,jv) + 2*vC(bv+4,jv)
            + 2*(v2(Qv(av,bv,jv)) - v2(Qv(av,bv,0))))

# 1. Verify Prop2 vs direct val(j)-val(0) via MN
print("=== Verify Prop 2 vs direct ===")
bad=0;tested=0
for m in range(4,16):
    n=2*m
    for av in range(4,n):
        bv=n-4-av
        if bv<4 or bv>av: continue
        M0=Mj((av,bv,4),0,m)
        if M0==0: continue
        val0 = 0 + 2*v2(comb(m,0)*M0)
        for jv in range(1,bv+1):
            M=Mj((av,bv,4),jv,m)
            if M==0: continue
            valj = jv + 2*v2(comb(m,jv)*M)
            d_direct = valj-val0
            d_closed = Delta_closed(av,bv,jv)
            tested+=1
            if d_direct!=d_closed:
                bad+=1
                if bad<8: print("  MISMATCH",(av,bv),jv,d_direct,d_closed)
print(f"tested={tested} mismatch={bad}")

# 2. Wide scan: is Delta(j)>0 for all interior j not in {0,2}, and Delta(2)>=0?
print("\n=== Wide scan Delta over closed form (a,b same parity, b>=4) ===")
viol_off=[]   # j not in {0,2} with Delta<=0
neg2=[]       # Delta(2)<0
mn4=10**9     # min Delta(4)
mn6=10**9; mn8=10**9
tie_at=set()
M=160
for av in range(4, M):
    for bv in range(4, av+1):
        if (av+bv)%2!=0: continue   # a,b same parity
        for jv in range(1,min(bv, 40)+1):
            d=Delta_closed(av,bv,jv)
            if jv==2:
                if d<0: neg2.append((av,bv,d))
                if d==0: tie_at.add(2)
            else:
                if d<=0:
                    if d==0: tie_at.add(jv)
                    if jv not in (2,) and len(viol_off)<30:
                        viol_off.append((av,bv,jv,d))
            if jv==4: mn4=min(mn4,d)
            if jv==6: mn6=min(mn6,d)
            if jv==8: mn8=min(mn8,d)
print("tie j-values seen:", sorted(tie_at))
print("Delta(2)<0 cases:", neg2[:10], "count",len(neg2))
print("Delta<=0 off {0,2} cases:", viol_off[:20], " total", len(viol_off))
print("min Delta(4)=",mn4," min Delta(6)=",mn6," min Delta(8)=",mn8)
