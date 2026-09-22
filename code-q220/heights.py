import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *

def alpha_fun(S,T,n):
    """alpha on the lift Z, base period computed from h(0)=0."""
    a={}; h=0
    for i in range(n):
        if i in T: h-=1
        if i in S: a[i]=h; h+=1
    d=len(S)-len(T)
    def A(j):
        q,r=divmod(j,n)
        return a[r]+q*d
    return A,d

def last_minimiser(S,T,x,n,A):
    """last S-position in the window (x, x+n) minimising alpha."""
    pts=[j for j in range(x+1,x+n) if j%n in S]
    mu=min(A(j) for j in pts)
    return max(j for j in pts if A(j)==mu)%n

# --- verify the bracket lemma: b_x = last minimiser of alpha over the window
bad=0; checked=0
for n in range(3,8):
    for S,T,_ in additive_pairs(n):
        A,d=alpha_fun(S,T,n)
        for x in X_set(S,T,n):
            r=ms_etilde(S,T,x,n)
            if r is None: continue
            checked+=1
            if r[0]!=last_minimiser(S,T,x,n,A): bad+=1
print(f"LEMMA (b_x = last alpha-minimiser in window): {checked} firings checked, {bad} failures")

# --- also: alpha_{b_x} == -|R_1|
bad2=0
for n in range(3,8):
    for S,T,_ in additive_pairs(n):
        A,d=alpha_fun(S,T,n)
        for x in X_set(S,T,n):
            L1,R1,_=ms_pairing(S,T,x,n)
            if not L1: continue
            pos={r:i for i,r in enumerate(ms_order(x,n))}
            b=min(L1,key=lambda r:pos[r])
            bl=[j for j in range(x+1,x+n) if j%n==b][0]
            # height at window start is 0
            if A(bl)!=-len(R1): bad2+=1
print(f"COROLLARY (alpha_b = -|R_1|): {bad2} failures")
