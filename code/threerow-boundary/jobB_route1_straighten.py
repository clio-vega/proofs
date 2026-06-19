"""
jobB_route1_straighten.py -- 2026-06-19 CODE, Job B route-1 (continuation).

The dimension-hook det D_q(a-j,b-j,c-j) VANISHES at the boundary (two negative
parts) while M_{b+i} != 0 (jobB_route1_qdet.py).  But the bialternant STRAIGHTENS:
vec = (a-j+2, b-j+1, c-j); sort descending with sign eps; nu = sorted - (2,1,0).
If nu is a partition with distinct vec entries, the f-continuation = eps * f^nu.

TEST: is |M_{b+i}| == f^nu for the straightened nu?  If yes, M is +-a genuine
DIMENSION, so its q-hook (a ratio of q-integers, leak=0) has tower == v2 EXACTLY,
and the full deficit tower reproduces g[c][k] -- route 1 CLOSES.  If no, report.
"""
import sympy as sp
from math import factorial
from functools import lru_cache

q = sp.symbols('q')

def v2int(n):
    n = abs(int(n))
    if n == 0:
        return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

# ---- alternant engine (self-contained, same as qdet) ----
def _f(n): return factorial(n)
def Ecoef(j, A, B, C):
    if A < 0 or B < 0 or C < 0 or (A + B + C) != 2*j: return 0
    s = A + B - C
    if s & 1: return 0
    al = s//2; be = B-al; ga = A-al
    if al < 0 or be < 0 or ga < 0 or al+be+ga != j: return 0
    return _f(j)//(_f(al)*_f(be)*_f(ga))
def _m3(n,i,j,k):
    if i<0 or j<0 or k<0 or i+j+k!=n: return 0
    return _f(n)//(_f(i)*_f(j)*_f(k))
def Tfast(p,qq,r,j,s):
    if p<0 or qq<0 or r<0 or p+qq+r!=2*j+s: return 0
    tot=0
    for P in range(s+1):
        for Q in range(s-P+1):
            R=s-P-Q; e=Ecoef(j,p-P,qq-Q,r-R)
            if e: tot+=_m3(s,P,Q,R)*e
    return tot
_S3=[((2,1,0),1),((1,2,0),-1),((2,0,1),-1),((0,2,1),1),((1,0,2),1),((0,1,2),-1)]
def Malt(a,b,c,j,m):
    s=2*m-2*j; base=(a+2,b+1,c); tot=0
    for (perm,sgn) in _S3:
        tot+=sgn*Tfast(base[0]-perm[0],base[1]-perm[1],base[2]-perm[2],j,s)
    return tot

def Ni_data(c,i,Pv,Bv,even):
    k=c-i; bv=2*Bv if even else 2*Bv+1; av=2*Pv+bv+c; m=Pv+bv+c; j=bv+i
    if bv<c+1 or j>m: return None
    M=Malt(av,bv,c,j,m); den=(bv-c+1)
    for t in range(2,i+1): den*=(bv+t)
    num=M*factorial(c)*factorial(k)
    if den==0 or num%den!=0: return None
    return num//den,M,den,av,bv,m,j

# ---- dimension of an arbitrary 3-vector via straightening ----
def straighten(vec):
    """vec = mu+delta (delta=(2,1,0)).  Return (eps, nu) with f-cont = eps*f^nu,
    or (0,None) if a repeat (=> alternant 0)."""
    import itertools
    pairs = list(enumerate(vec))
    # sign of permutation that sorts vec descending
    s = sorted(vec, reverse=True)
    if len(set(vec)) < len(vec):
        return 0, None
    # parity: count inversions needed
    perm = sorted(range(len(vec)), key=lambda i: -vec[i])
    inv = sum(1 for a in range(len(perm)) for b in range(a+1,len(perm)) if perm[a]>perm[b])
    eps = (-1)**inv
    delta = [2,1,0]
    nu = [s[i]-delta[i] for i in range(len(vec))]
    if any(x < 0 for x in nu) or any(nu[i] < nu[i+1] for i in range(len(nu)-1)):
        return 0, None   # not a partition
    return eps, tuple(x for x in nu)

@lru_cache(maxsize=None)
def fdim(mu):
    mu = tuple(x for x in mu if x > 0)
    n = sum(mu)
    if n == 0: return 1
    L = len(mu)
    detr = sp.Matrix([[sp.Rational(1, factorial(mu[i]-(i+1)+(jc+1)))
                       if mu[i]-(i+1)+(jc+1) >= 0 else 0
                       for jc in range(L)] for i in range(L)]).det()
    return int(factorial(n)*detr)

# ---- q-hook of a partition nu: f^nu_q = [n]_q!/prod[hook]_q (ratio of q-ints) ----
def qint(n):
    n=int(n); return sp.Integer(0) if n<=0 else sp.expand(sum(q**i for i in range(n)))
def qfact(n):
    r=sp.Integer(1)
    for i in range(1,int(n)+1): r*=qint(i)
    return sp.expand(r)
def hooks(nu):
    nu=[x for x in nu if x>0]; conj=[sum(1 for r in nu if r>cidx) for cidx in range(nu[0])] if nu else []
    H=[]
    for i,row in enumerate(nu):
        for jc in range(row):
            arm=row-jc-1; leg=conj[jc]-i-1; H.append(arm+leg+1)
    return H
_PHI={r:sp.Poly(sp.cyclotomic_poly(2**r,q),q) for r in range(1,16)}
def tower_qhook(nu):
    """Sum_r mult_{Phi_2^r} of f^nu_q = [n]!_q / prod[h]_q, via additive valuations
    on q-integers (each [n]_q = prod_{d|n,d>1} Phi_d) -> EXACT, no polynomial div."""
    from sympy import divisors
    def qint_tower(n):  # Sum_r mult_{Phi_2^r}([n]_q) = #{r>=1: 2^r | n}
        return sum(1 for r in range(1,16) if n % (2**r)==0)
    n=sum(nu)
    tot=sum(qint_tower(t) for t in range(1,n+1))       # [n]!_q
    for h in hooks(nu):
        tot-=qint_tower(h)
    return tot

def v2_qfacttower(n):  # tower of [n]!_q
    return sum(sum(1 for r in range(1,16) if t%(2**r)==0) for t in range(1,n+1))

print("="*80)
print("Route-1 straightening test: is |M_{b+i}| a genuine (straightened) dimension f^nu?")
print("  If yes -> q-hook of nu is leak-free -> tower(N^q)=g exactly -> route 1 CLOSES.")
print("="*80)
print(f"  {'c':>2} {'k':>2} {'i':>2} {'sl':>4} | {'g':>3} {'twN':>4} | "
      f"{'|M|=f^nu?':>9} {'nu':>14} | match")

def content_min(c,i,even,Bmax=40,Pmax=40):
    g=10**9; best=None
    for Bv in range((c+1)//2+1,Bmax+1):
        for Pv in range(0,Pmax+1):
            res=Ni_data(c,i,Pv,Bv,even)
            if res is None: continue
            vv=v2int(res[0])
            if vv<g: g=vv; best=res
    return g,best

rows=[]
for c in range(4,10):
    for k in (4,5,6):
        i=c-k
        if i<1: continue
        for even in (True,False):
            g,res=content_min(c,i,even)
            if res is None: continue
            Nval,M,den,av,bv,m,j=res
            vec=(av-j+2, bv-j+1, c-j)
            eps,nu=straighten(vec)
            sl='even' if even else 'odd'
            if nu is None:
                isdim=False; twN=None
            else:
                fnu=fdim(nu)
                isdim=(abs(M)==fnu)
                # full deficit q-tower: tower(f^nu_q) + tower([c]![k]!) - tower(den_q)
                den_tower=sum(1 for t in [bv-c+1]+[bv+tt for tt in range(2,i+1)]
                              for r in range(1,16) if t%(2**r)==0)
                tw_M=tower_qhook(nu)
                tw_ck=v2_qfacttower(c)+v2_qfacttower(k)
                twN=tw_M+tw_ck-den_tower
            match=(twN==g)
            rows.append((c,k,i,sl,g,twN,isdim,nu,match,M))
            print(f"  {c:>2} {k:>2} {i:>2} {sl:>4} | {g:>3} {str(twN):>4} | "
                  f"{str(isdim):>9} {str(nu):>14} | {'YES' if match else 'no'}")

ndim=sum(1 for r in rows if r[6])
nmatch=sum(1 for r in rows if r[8])
print(f"\n  |M| == f^nu (M is a straightened dimension): {ndim}/{len(rows)}")
print(f"  tower(N^q) == g[c][k]: {nmatch}/{len(rows)}")
if ndim==len(rows) and nmatch==len(rows):
    print("\n  VERDICT: YES -- M_{b+i} is +-a straightened dimension f^nu; its q-hook is a")
    print("  ratio of q-integers (leak=0), so the cyclotomic tower of the deficit equals")
    print("  the exact content g[c][k].  Route 1 closes the Content Lemma (single q-hook).")
elif ndim==len(rows):
    print("\n  VERDICT: PARTIAL -- |M|=f^nu but tower!=g; den/c!k! bookkeeping mismatch (investigate).")
else:
    print("\n  VERDICT: NO/PARTIAL -- M is NOT always a straightened dimension; alternant")
    print("  carries content beyond a single q-hook on the failing rows.")
