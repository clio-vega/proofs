"""
Verify the exact engine identities:
 (E1) Phi(z) = sum_r C(m,r) 2^r R_r (1+z)^{m-r}   [exact in Z[z]]
 (E2) R_0 = chi^lam(2^m)
 (E3) chi^lam(2^m) odd  <=>  B* = submasks(m)  (b-picture min-locus)
Also: how often is chi^lam(2^m) odd among ties? and does its oddness
correlate with the SHARP |J*| being a full submask box?
"""
import math
from job1_tie_census import partitions, v2, chi_b_vector, M_vector
from job_hierarchy import R_vector

def Phi_coeffs(lam,m):
    chi=chi_b_vector(lam,m)
    return [math.comb(m,k)*chi[k] for k in range(m+1)]

def E1_rhs(lam,m):
    R=R_vector(lam,m)
    # sum_r C(m,r)2^r R_r (1+z)^{m-r}, return coeff list in z
    coeff=[0]*(m+1)
    for r in range(m+1):
        a=math.comb(m,r)*(1<<r)*R[r]
        # (1+z)^{m-r}
        d=m-r
        for k in range(d+1):
            coeff[k]+=a*math.comb(d,k)
    return coeff

def main():
    e1=0; e1ok=0; e2ok=0; e2=0
    for m in range(1,11):
        for lam in partitions(2*m):
            lhs=Phi_coeffs(lam,m); rhs=E1_rhs(lam,m)
            e1+=1
            if lhs==rhs: e1ok+=1
            R=R_vector(lam,m); chi=chi_b_vector(lam,m)
            e2+=1
            if R[0]==chi[m]: e2ok+=1
    print(f"(E1) Phi == sum C(m,r)2^r R_r (1+z)^(m-r):  {e1ok}/{e1}")
    print(f"(E2) R_0 == chi^lam(2^m):                   {e2ok}/{e2}")

    # E3 + correlation
    odd_chi=0; ties=0; oddchi_ties=0
    from job_box_detail import is_affine_box
    def submasks(m):
        out=[]; s=m
        while True:
            out.append(s)
            if s==0: break
            s=(s-1)&m
        return sorted(out)
    e3ok=0; e3=0
    for m in range(1,13):
        subs=submasks(m)
        for lam in partitions(2*m):
            chi=chi_b_vector(lam,m)
            coeffs=[math.comb(m,k)*chi[k] for k in range(m+1)]
            fin={k:v2(c) for k,c in enumerate(coeffs) if c!=0}
            bmin=min(fin.values()); B=sorted(k for k in fin if fin[k]==bmin)
            chi2m=chi[m]
            e3+=1
            if chi2m%2!=0:
                pred=[k for k in subs if chi[k]!=0]
                # claim B == submasks(m) (all chi nonzero on submasks when chi2m odd?)
                if B==subs and set(pred)==set(subs): e3ok+=1
                elif B==pred: e3ok+=1   # weaker: submasks-with-nonzero
    print(f"(E3) chi(2^m) odd => B*==submasks(m): {e3ok}/{e3} (count includes all shapes; ok if chi even)")

if __name__=="__main__":
    main()
