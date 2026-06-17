"""
generalc_certify.py -- prove session 2026-06-17

Test that the SAME hand-bound mechanism closes the boundary lemma for c=4..8.
Uniform inputs:
  Delta(b+i) = (eta-k) - 2 v2(k!) + 2 v2(N_i) - 2 v2(a-c+2)
               - 2 [k<=c-2] v2(a-b+1) + 2 v2(Pi_i),   k=c-i.
Hand lower bound: even deficits cost the "+1" (-2); their /2 factor is either
  IN Pi_i range (absorbed by Lemma P: remaining v2 >= floor(L/2)-d_in) or
  OUT of range (compensated by v2(N_i): comp = v2(N_i) - sum out v2(f) >= 0).
v2(N_i) computed EXACTLY from MN value of M_{b+i}.
theta: 0 unless (c odd and a odd) -> 3.
"""
from math import comb, factorial
from mn import Mj

def v2(n):
    n=abs(int(n))
    if n==0: return 999
    k=0
    while n%2==0: n//=2;k+=1
    return k

def v2_Ni(P,bb,c,i):
    aa=2*P+bb+c; m=P+bb+c; k=c-i
    Mbi=Mj((aa,bb,c),bb+i,m)
    if Mbi==0: return 999
    base=v2(Mbi)+v2(factorial(c))+v2(factorial(k))-v2(bb-c+1)
    if i==1: return base-v2(aa-bb+1)
    return base-sum(v2(bb+t) for t in range(2,i+1))

def true_delta(P,bb,c,i):
    aa=2*P+bb+c; m=P+bb+c; lam=(aa,bb,c)
    M0=Mj(lam,0,m)
    if M0==0: return None
    j=bb+i
    if j>m: return None
    try: Mji=Mj(lam,j,m)
    except AssertionError: return None
    if Mji==0: return None
    return (j+2*v2(comb(m,j)*Mji)) - 2*v2(M0)

def hand_bound(P,bb,c,i):
    aa=2*P+bb+c; k=c-i; w=bb+c
    eta=1 if w%2==1 else 2
    ell=-(-(w+3)//2); L=ell-k-1
    base=(eta-k)-2*v2(factorial(k))
    lo,hi=P+k+1, P+ell-1
    d_in=0; out_facs=[]
    if bb%2==0:                      # a-c+2 = 2(P+b/2+1)
        f1=P+bb//2+1; base-=2
        if lo<=f1<=hi: d_in+=1
        else: out_facs.append(f1)
    if i>=2 and c%2==1:              # a-b+1 = 2(P+(c+1)/2)
        f2=P+(c+1)//2; base-=2
        if lo<=f2<=hi: d_in+=1
        else: out_facs.append(f2)
    elif i>=2 and c%2==0:
        pass                         # a-b+1 odd, no deficit
    piterm=2*(L//2-d_in)
    comp=v2_Ni(P,bb,c,i)-sum(v2(f) for f in out_facs)
    return base+piterm+2*comp, comp

for c in [4,5,6,7,8]:
    b0=2*c
    checked=0; bad_valid=0; bad_suff=0; minmarg=999
    res={}
    for P in range(0,45):
        for bb in range(b0, 80):
            aa=2*P+bb+c; m=P+bb+c
            if aa<bb or m>110: continue
            theta=3 if (c%2==1 and aa%2==1) else 0
            for i in range(1,c+1):
                td=true_delta(P,bb,c,i)
                if td is None: continue
                hb,comp=hand_bound(P,bb,c,i)
                checked+=1
                if hb>td: bad_valid+=1
                if not (hb>-theta): bad_suff+=1
                minmarg=min(minmarg, hb+theta)
                key=(aa%2,i); r=res.get(key,[999,999])
                r[0]=min(r[0],td+theta); r[1]=min(r[1],hb+theta); res[key]=r
    print("c=%d (b>=%d, m<=110): checked=%d  invalid=%d  insufficient=%d  min hand margin=%d"
          %(c,b0,checked,bad_valid,bad_suff,minmarg))
