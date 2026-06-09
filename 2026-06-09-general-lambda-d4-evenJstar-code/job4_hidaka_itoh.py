"""
JOB 4 — Hidaka-Itoh (arXiv 2403.10817) orientation probe.  [numerical, fast]

Their statement: the principal specialization s_lambda evaluated at primitive n-th
roots of unity lies in {-1,0,1} whenever n has <= 2 distinct odd prime factors.

We (a) confirm my fiber value G_lambda(zeta_4) equals the fake-degree / q-hook
polynomial  f^lambda(q) = sum_T q^maj(T) = [n]_q!/prod_cells [h]_q  at q=i, then
(b) tabulate f^lambda(zeta_d) over all lambda |- N for a range of d, classifying d by
(#distinct odd prime factors) and by d mod 4, to compare the {-1,0,1} bound with my
'RICH iff d == 2 (mod 4)' trichotomy.

Note: f^lambda is the CSP COUNT normalization (can be large); Hidaka-Itoh's s_lambda
is a DIFFERENT (principal-specialization) normalization.  The probe is about the
arithmetic GATES (<=2 odd prime factors  vs  d==2 mod4), not the literal value set.
"""
import math, cmath
from sympy import primefactors
from job1_tie_census import partitions, M_vector, G_gaussian

def hook_lengths(lam):
    conj=[sum(1 for r in lam if r>c) for c in range(lam[0])]
    hs=[]
    for i,row in enumerate(lam):
        for jcol in range(row):
            hs.append((row-jcol-1)+(conj[jcol]-i-1)+1)
    return hs

def fake_degree_coeffs(lam):
    """f^lambda(q) integer coeff list via q-factorials, exact integer poly arithmetic."""
    n=sum(lam)
    def qint(k): return [1]*k                  # [k]_q = 1+q+..+q^{k-1}
    def pmul(a,b):
        c=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b): c[i+j]+=x*y
        return c
    def pdiv(a,b):                              # exact poly division a/b
        a=a[:]; q=[0]*(len(a)-len(b)+1)
        for i in range(len(q)-1,-1,-1):
            coef=a[i+len(b)-1]//b[-1]; q[i]=coef
            for j,y in enumerate(b): a[i+j]-=coef*y
        assert all(x==0 for x in a), "non-exact"
        return q
    num=[1]
    for k in range(1,n+1): num=pmul(num,qint(k))
    den=[1]
    for h in hook_lengths(lam): den=pmul(den,qint(h))
    return pdiv(num,den)

def eval_poly_root(coeffs,d):
    z=cmath.exp(2j*math.pi/d)
    s=0j; zp=1
    for c in coeffs:
        s+=c*zp; zp*=z
    return s

def round_gauss(z,tol=1e-6):
    re=round(z.real); im=round(z.imag)
    if abs(z.real-re)<tol and abs(z.imag-im)<tol: return (re,im)
    return None  # not a Gaussian integer at this precision

def confirm_equals_G():
    bad=0; tot=0
    for m in range(1,8):
        for lam in partitions(2*m):
            reA,imA=G_gaussian(lam,m,M_vector(lam,m))
            v=eval_poly_root(fake_degree_coeffs(lam),4)
            g=round_gauss(v); tot+=1
            if g!=(reA,imA):
                bad+=1
                if bad<=8: print("  G!=f at",lam,(reA,imA),g,v)
    print(f"confirm G_lambda(i)==f^lambda(i), m<=7: {tot-bad}/{tot} match")

if __name__=="__main__":
    confirm_equals_G()
    print()
    for N in (8,10):
        print(f"=== value behaviour of f^lambda(zeta_d), lambda |- {N} ===")
        print(f"{'d':>3} {'oddpf':>5} {'dmod4':>5}  {'real-int?':>9}  {'|val| range':>12}  {'#zeros':>6}")
        for d in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,21,22,30,33,35,45]:
            mags=[]; allrealint=True; zeros=0; realvals=[]
            for lam in partitions(N):
                v=eval_poly_root(fake_degree_coeffs(lam),d)
                g=round_gauss(v)
                if g is None: allrealint=False; mags.append(abs(v)); continue
                re,im=g
                if im!=0: allrealint=False
                else: realvals.append(re)
                if (re,im)==(0,0): zeros+=1
                mags.append(math.hypot(re,im))
            opf=len([p for p in primefactors(d) if p%2==1])
            rng=f"{min(mags):.0f}..{max(mags):.0f}"
            print(f"{d:>3} {opf:>5} {d%4:>5}  {str(allrealint):>9}  {rng:>12}  {zeros:>6}")
        print()
