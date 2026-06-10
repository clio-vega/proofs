"""
SHARP box via mod-2.  G = Psi(pi), Psi(Y)=sum_r C(m,r)R_r Y^r, basis Y=pi=1+i.
val(r)=r+2v2(C(m,r)R_r).  Within the parity class of J*, the min-locus (=J*) is
the bottom-mod-2 support of the scaled polynomial.  Test if that support is a
shifted power of (1+y).

For even class r=2s: scaled poly = sum_s c_{2s} 2^s y^s, c_r=C(m,r)R_r.
For odd  class r=2s+1: factor out one pi; scaled = sum_s c_{2s+1} 2^s y^s.
bottom-mod-2 support should be y^{l0}(1+y)^g.
"""
import math
from job1_tie_census import partitions, M_vector, v2
from job_hierarchy import R_vector
from job_mod2_engine import is_shift_pow_1pz, poly_mod2

def sharp_box_test(lam,m):
    R=R_vector(lam,m)
    c=[math.comb(m,r)*R[r] for r in range(m+1)]
    val={r:r+2*v2(c[r]) for r in range(m+1) if c[r]!=0}
    if not val: return None
    mu=min(val.values()); J=sorted(r for r in val if val[r]==mu)
    par=J[0]%2
    # build scaled poly over the parity class: index s, coeff c_{2s+par}*2^s
    coeffs=[]
    s=0
    while 2*s+par<=m:
        r=2*s+par
        coeffs.append(c[r] if r<=m else 0)
        s+=1
    # multiply coeff_s by 2^s
    scaled=[(coeffs[s]*(1<<s)) for s in range(len(coeffs))]
    nz=[v2(x) for x in scaled if x!=0]
    if not nz: return None
    e=min(nz)
    red=[x//(1<<e) for x in scaled]
    sup=poly_mod2(red)
    ok,l0,g=is_shift_pow_1pz(sup)
    # also recover J in s-coords to confirm sup==Jlocus
    Js=set((r-par)//2 for r in J)
    return ok, sup, Js, l0, g

def main():
    tot=0; ok_c=0; bad=[]
    for m in range(1,13):
        for lam in partitions(2*m):
            res=sharp_box_test(lam,m)
            if res is None: continue
            ok,sup,Js,l0,g=res
            # only count genuine ties
            R=R_vector(lam,m); c=[math.comb(m,r)*R[r] for r in range(m+1)]
            val={r:r+2*v2(c[r]) for r in range(m+1) if c[r]!=0}
            mu=min(val.values()); J=[r for r in val if val[r]==mu]
            if len(J)<2: continue
            tot+=1
            # sup should equal Js (the min locus), and be a shifted (1+y)^g
            if ok and sup==Js: ok_c+=1
            else: bad.append((m,lam,sorted(sup),sorted(Js),ok))
    print(f"SHARP ties (m<=12): {tot}")
    print(f"  J* (sharp) bottom-mod2 == y^l0 (1+y)^g AND == min-locus: {ok_c}/{tot}")
    for x in bad[:15]: print("   bad:",x)

if __name__=="__main__":
    main()
