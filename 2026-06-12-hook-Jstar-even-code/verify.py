"""
Verification for: 2026-06-12-hook-Jstar-even.md
Theorem: for every hook lambda=(a,1^b) |- 2m, |J*| in {1,2} (hence even when >=2).
Checks, end to end:
  (1) M_j = C(2m-1-j, a-1)  via brute-force symmetric-function Pieri chains.
  (2) Prop 2 identity: val(j)-val(0) = j + 2[v2C(m,j)+v2C(b,j)-v2C(2m-1,j)].
  (3) g(j)>=0 all j, g(j)>0 for j not in {0,2}  ==> J* subset {0,2}.
  (4) Lemma K and Lemma K' bounds.
"""
from sympy import binomial
from collections import defaultdict

def v2(n):
    if n==0: return None
    c=0
    while n%2==0: n//=2;c+=1
    return c
def C(n,k):
    if k<0 or k>n or n<0: return 0
    return int(binomial(n,k))

# ---- (1) brute force M_j via Pieri ----
def add_single(mu):
    res=[];mu=list(mu);L=len(mu)
    for i in range(L):
        if i==0 or mu[i-1]>mu[i]:
            nu=mu[:];nu[i]+=1;res.append(tuple(nu))
    res.append(tuple(mu+[1]));return res
def add_vstrip2(mu):
    mu=list(mu);ext=mu+[0,0];n=len(ext);res=[]
    for i in range(n):
        for k in range(i+1,n):
            nu=ext[:];nu[i]+=1;nu[k]+=1
            if all(nu[t]>=nu[t+1] for t in range(n-1)):
                while nu and nu[-1]==0: nu.pop()
                res.append(tuple(nu))
    return res
def mult(state,op):
    new=defaultdict(int)
    for mu,c in state.items():
        for nu in op(mu): new[nu]+=c
    return new
def Mj_brute(lam,j,m):
    state={():1}
    for _ in range(j): state=mult(state,add_vstrip2)
    for _ in range(2*(m-j)): state=mult(state,add_single)
    return state.get(tuple(lam),0)

ok1=True
for m in range(1,7):
    for a in range(1,2*m+1):
        b=2*m-a; lam=tuple([a]+[1]*b) if b>0 else (a,)
        for j in range(m+1):
            if Mj_brute(lam,j,m)!=C(2*m-1-j,a-1): ok1=False
print("(1) M_j = C(2m-1-j,a-1)  [m<=6]:", ok1)

# ---- (2),(3) Prop 2 + g(j) structure, and conclude J* subset {0,2} ----
def carries(x,y):
    c=0;cy=0
    while x>0 or y>0 or cy>0:
        s=(x&1)+(y&1)+cy; cy=s>>1; c+=cy; x>>=1;y>>=1
    return c
def vC(n,k):
    if k<0 or k>n: return None
    return carries(k,n-k)

ok2=ok3=True; maxJ=0
for m in range(2,300):
    for a in range(1,2*m+1):
        b=2*m-a; r=a-1
        def val(j):
            Mj=C(2*m-1-j,r); cj=C(m,j)
            if Mj==0 or cj==0: return None
            return j+2*(v2(cj)+v2(Mj))
        v0=val(0); J=[0]
        for j in range(1,b+1):
            vj=val(j)
            if vj is None: continue
            g=vj-v0
            if g!=j+2*(vC(m,j)+vC(b,j)-vC(2*m-1,j)): ok2=False   # Prop 2
            if g<0: ok3=False                                    # 0 is min
            if g==0 and j!=0 and j!=2: ok3=False                 # only {0,2} tie
            if g>0 and j==0: ok3=False
            if g==0: J.append(j)
        maxJ=max(maxJ,len(J))
        if any(jj not in (0,2) for jj in J): ok3=False
print("(2) Prop 2 identity [m<300]:", ok2)
print("(3) g>=0, J* subset {0,2}, max|J*| =", maxJ, "  [m<300]:", ok3 and maxJ<=2)

# ---- (4) Lemma K / K' ----
okK=okKp=True
for m in range(2,2000):
    for t in range(1,m//2+1):
        if vC(m-1,t)-vC(m,2*t) > t: okK=False
    for s in range(0,(m-1)//2+1):
        if 2*s+1<=m and vC(m-1,s)-vC(m,2*s+1) > s: okKp=False
print("(4) Lemma K (<=t) [m<2000]:", okK, " | Lemma K' (<=s) [m<2000]:", okKp)
