"""INDEPENDENT verification of the m=2 failures.
Instrument 2 never touches the y-box or the l1 formula: it builds f_nu from
gen.f_nu (interval convolutions) and reads the half-width off as (len(co)-1)/2,
exactly as mprofile.py did on 10-02.  Then it forms M(r) directly.
Calibration: must reproduce a KNOWN value first (brief control 2)."""
import gen, tent
from collections import Counter
from layers import is_concave_pospart

def M_direct(mu,lam,n,m,b):
    """M(r) from raw f_nu supports -- no y-box, no l1 formula."""
    hs=[]
    for nu in gen.slice_nus(mu,lam,n,m,b):
        off,co=gen.f_nu(nu,lam,n,m)
        if co: hs.append((len(co)-1)/2)
    if not hs: return None,None
    lo,hi=min(hs),max(hs)
    return lo,[sum(1 for h in hs if abs(h-(lo+j))<1e-9) for j in range(int(round(hi-lo))+1)]

# ---- CALIBRATION on values already held (paper sec:verif) ----
print("CALIBRATION")
off,co=gen.slice_sum((0,3),(4,7),8,2,2)
print("  k(.,2) at n=8,mu=(0,3),lam=(4,7):", co, " expected (1,3,6,7,6,3,1)")
off,co=gen.slice_sum((0,5),(5,12),10,2,4)
print("  G at n=10,mu=(0,5),lam=(5,12),b=4:", co, " expected (1,4,7,10,11,10,7,4,1)")
# refusal probe on the concavity test
print("  is_concave_pospart((3,6,6,3,1)) must be False:", is_concave_pospart([3,6,6,3,1]))
print("  is_concave_pospart((2,1,1)) must be False:", is_concave_pospart([2,1,1]))
print("  is_concave_pospart((1,3,6,10,6,3,1)) must be True:", is_concave_pospart([1,3,6,10,6,3,1]))

# ---- smallest witness search, by (n, d, lam) ----
print("\nSEARCH for smallest m=2 witness with the INDEPENDENT instrument")
found=[]
for n in range(2,13):
    for (mu,lam) in gen.pairs(n,2,n+6):
        d=sum(lam)-sum(mu)
        for b in range(0,d+1):
            lo,M=M_direct(mu,lam,n,2,b)
            if M is None: continue
            if not is_concave_pospart(M):
                found.append((n,d,mu,lam,b,lo,M))
    if found: break
print(f"  first failing level n={n}: {len(found)} witnesses")
found.sort(key=lambda t:(t[1],t[3]))
for w in found[:6]: print("   ",w)

# ---- full detail of the smallest ----
n,d,mu,lam,b,lo,M = found[0]
print(f"\nDETAIL  n={n} mu={mu} lam={lam} d={d} b={b}  c=(d-b)/2={(d-b)/2}")
for nu in gen.slice_nus(mu,lam,n,2,b):
    lr=gen.LR(nu,lam,n,2); off,co=gen.f_nu(nu,lam,n,2)
    w=[R-L+1 for (L,R) in lr]
    print(f"   nu={nu}  (L,R)={lr}  w={w}  Lambda={(len(co)-1)/2 if co else None}  f_nu={co} on [{off},{off+len(co)-1}]" if co else f"   nu={nu}  (L,R)={lr}  w={w}  EMPTY")
print(f"   M = {M}  (r from {lo})   concave_pospart = {is_concave_pospart(M)}")
off,co=gen.slice_sum(mu,lam,n,2,b)
print(f"   G = k(.,{b}) = {co} on [{off},{off+len(co)-1}]   PF2 = {gen.is_pf2(co)}")
