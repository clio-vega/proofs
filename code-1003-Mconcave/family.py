"""The explicit all-m counterexample family, and the closed-form data.

For m>=3, n=m+6, b=3:
    mu = (0,1,...,m-2, m+2),        lam = (0,1,...,m-3, m+3, m+5)
giving  g = (0,...,0,5,1),  delta = (0,...,0,-2),  d=8,  sigma=1,
        P = (0,...,0,0,-2),  Q = (0,...,0,3,1),
effective slice {(0,..,0,y_{m-1},y_m)}: (0,1),(1,0),(2,-1),(3,-2) -> k=0,0,1,2,
so N=(2,1,1) and M=(1,1,2) at r=1/2,3/2,5/2 with c=(d-b)/2=5/2.
m=2: n=8, mu=(0,4), lam=(3,9), b=3 (the same data cyclically rotated).
Checked here with the RAW f_nu instrument (gen.LR / gen.f_nu), not the box."""
import gen, tent
from layers import is_concave_pospart
print(f"{'m':>3} {'n':>3} {'mu':>28} {'lam':>28} {'d':>3} {'b':>2}  {'g':>14} {'delta':>14} {'M':>12} concave  G=PF2")
for m in range(2,13):
    n=m+6; b=3
    if m==2: mu,lam=(0,4),(3,9)
    else:
        mu=tuple(list(range(0,m-1))+[m+2])
        lam=tuple(list(range(0,m-2))+[m+3,m+5])
    import cyl
    assert cyl.is_shape(mu,n,m) and cyl.is_shape(lam,n,m), (m,mu,lam)
    assert all(mu[i]<=lam[i] for i in range(m)), (m,mu,lam)
    d=sum(lam)-sum(mu)
    g=[tent.lam_at(lam,i,n,m)-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
    dl=[mu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
    hs=[]
    for nu in gen.slice_nus(mu,lam,n,m,b):
        off,co=gen.f_nu(nu,lam,n,m)
        if co: hs.append((len(co)-1)/2)
    lo,hi=min(hs),max(hs)
    M=[sum(1 for h in hs if abs(h-(lo+j))<1e-9) for j in range(int(round(hi-lo))+1)]
    off,co=gen.slice_sum(mu,lam,n,m,b)
    print(f"{m:>3} {n:>3} {str(mu):>28} {str(lam):>28} {d:>3} {b:>2}  {str(g):>14} {str(dl):>14} {str(M):>12} {str(is_concave_pospart(M)):>7}  {gen.is_pf2(co)}")
