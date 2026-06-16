from mn import Mj
from math import comb, factorial
def v2(n):
    n=abs(n)
    if n==0: return 99
    k=0
    while n%2==0: n//=2; k+=1
    return k
def N3(a,b): return a*(b+4)-(b**2+3*b+8)
def N2(a,b): return a**2*b**2+7*a**2*b+12*a**2-2*a*b**3-9*a*b**2-21*a*b-20*a+b**4+2*b**3+13*b**2+16*b-8
def N1(a,b): return a**2*b**3+9*a**2*b**2+26*a**2*b+24*a**2-2*a*b**4-7*a*b**3-13*a*b**2-14*a*b+b**5-2*b**4+17*b**3-4*b**2-84*b-96

# ---- (1) verify content bounds hold for ALL shapes in range ----
bad_content=[]
for m in range(8,200):
    for b in range(6,m):
        a=2*m-4-b
        if a<b: continue
        if (a+b)%2!=0: continue  # same parity guaranteed since a+b=2m-4
        n2=v2(N2(a,b)); n1=v2(N1(a,b)); n3=v2(N3(a,b))
        if n2<2: bad_content.append(('N2',a,b,n2))
        if b%2==0:
            if n3<1: bad_content.append(('N3even',a,b,n3))
            if n1<3: bad_content.append(('N1even',a,b,n1))
        else:
            if n1<2: bad_content.append(('N1odd',a,b,n1))
print("content-bound violations (m<200):", len(bad_content), bad_content[:8])

# ---- (2) CERTIFIED lower bound on Delta using ONLY proven facts ----
# proven product facts:
#  v2(prod of Q consecutive ints) >= floor(Q/2)  (count of evens)
#  removing 1 factor: >= floor(Q/2)-1
def evens_lb(Q): return Q//2
def Delta_certified_lb(a,b,m,i):
    P=m-b-4
    if b%2==1: # b odd
        L={1:(b-3)//2+1,2:(b-3)//2+2,3:(b-3)//2+3,4:(b-3)//2+4}[i]
        # Delta = (i-3) - 2c + 2 v2(N_i) + 2 v2(Pi_i)
        c={4:0,3:0,2:1,1:1}[i]
        nlb={4:0,3:0,2:2,1:2}[i]   # v2(N): N4=0;N3 odd=0;N2>=2;N1>=2
        prodlb=evens_lb(L)         # >=1 even when L>=2
        return (i-3)-2*c+2*nlb+2*prodlb
    else: # b even, absorb factor P+b/2+1 -> Pi' has L-1 terms
        L={1:b//2+0,2:b//2+1,3:b//2+2,4:b//2+3}[i]  # L_i=b/2+i-1
        c={4:0,3:0,2:1,1:1}[i]
        nlb={4:0,3:1,2:2,1:3}[i]   # N4=0;N3 even>=1;N2>=2;N1 even>=3
        prodlb=max(0,evens_lb(L)-1)  # remove one factor
        return (i-4)-2*c+2*nlb+2*prodlb

# verify certified lb <= true Delta AND certified lb >= 1
fail=[]
for m in range(8,120):
    for b in range(6,m):
        a=2*m-4-b
        if a<b: continue
        P=m-b-4
        lam=(a,b,4)
        M0=Mj(lam,0,m)
        if M0==0: continue
        v0=2*v2(comb(m,0)*M0)
        for i in [1,2,3,4]:
            if b+i>m: continue
            Mj_=Mj(lam,b+i,m)
            if Mj_==0: continue
            true_D=(b+i)+2*v2(comb(m,b+i)*Mj_)-v0
            lb=Delta_certified_lb(a,b,m,i)
            if lb>true_D: fail.append(('LB>true',a,b,m,i,lb,true_D))
            if lb<1: fail.append(('LB<1',a,b,m,i,lb,true_D))
print("certification failures (m<120):", len(fail), fail[:8])
