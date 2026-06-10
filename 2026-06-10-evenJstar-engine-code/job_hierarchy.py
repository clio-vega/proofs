"""
Compare the various pictures' polygon-minima against the TRUE v_pi(G).
- j-picture (M_j, basis (i-1)^j):  G=sum_j C(m,j)M_j (i-1)^j;  val_j = j+2v2(C(m,j)M_j)
- r-picture (R_r, basis pi^r):     G=sum_r C(m,r)R_r pi^r;     val_r = r+2v2(C(m,r)R_r)
- b-picture (chi_k, basis (-i)^k):  T=sum_k C(m,k)chi_k(-i)^k;  beta-min*2 (+(-m) for G)
True v_pi(G) computed exactly in Z[i].
R_r = <s_lam, p2^{m-r} e2^r>.
"""
import math
from job1_tie_census import partitions, M_vector, v2, chi_b_vector

def R_vector(lam,m):
    """R_r = <s_lam, p2^{m-r} e2^r>, r=0..m.  e2=(p1^2-p2)/2 => expand to chi."""
    chi=chi_b_vector(lam,m)  # chi_b = chi^lam(2^b 1^{2m-2b})
    R=[]
    for r in range(m+1):
        # p2^{m-r} e2^r = 2^{-r} sum_a C(r,a)(-1)^a p2^{m-r+a} p1^{2(r-a)}
        # <s_lam, p2^{b} p1^{2m-2b}> = chi_b
        s=0
        for a in range(r+1):
            b=m-r+a
            s+= math.comb(r,a)*(-1)**a * chi[b]
        assert s % (1<<r)==0, (lam,r,s)
        R.append(s//(1<<r))
    return R

def vpi_G(lam,m):
    """exact v_pi(G), G=sum_r C(m,r)R_r pi^r, pi=1+i; v_pi via Z[i]."""
    R=R_vector(lam,m)
    # build G in Z[i]
    # pi^r in Z[i]
    from functools import reduce
    G=0  # complex integer as (re,im)
    re,im=0,0
    pr_re,pr_im=1,0  # pi^0
    for r in range(m+1):
        c=math.comb(m,r)*R[r]
        re+=c*pr_re; im+=c*pr_im
        # multiply pi=(1+i)
        nre=pr_re-pr_im; nim=pr_re+pr_im
        pr_re,pr_im=nre,nim
    # v_pi(re+im i) = v2(re^2+im^2)
    if re==0 and im==0: return None,(re,im)
    N=re*re+im*im
    return v2(N),(re,im)

def main():
    print("shape           | true vpiG | j-min mu | r-min | (j==true?) (r==true?)")
    for lam,m in [((4,1,1),3),((4,2,2),4),((6,3,1,1,1),6),((8,2,2),6),((6,6),6),
                  ((4,4,2,2),6),((3,3,1,1),4),((2,2),2),((6,3,2,1),6),((5,4,3),6)]:
        Ms=M_vector(lam,m); R=R_vector(lam,m)
        valj=[j+2*v2(math.comb(m,j)*Ms[j]) for j in range(m+1) if Ms[j]!=0]
        valr=[r+2*v2(math.comb(m,r)*R[r]) for r in range(m+1) if R[r]!=0]
        muj=min(valj); mur=min(valr) if valr else None
        vp,_=vpi_G(lam,m)
        print(f"{str(lam):15s} | {str(vp):9s} | {muj:8d} | {str(mur):5s} | {vp==muj!s:5s} {vp==mur!s:5s}")

    # full check: is r-picture min-locus an affine box?  is j==r min-loci?
    from job_box_detail import is_affine_box
    tot=0; rbox=0; jeqr=0; jsharp=0; rsharp=0
    for m in range(1,12):
        for lam in partitions(2*m):
            Ms=M_vector(lam,m); R=R_vector(lam,m)
            vj={j:j+2*v2(math.comb(m,j)*Ms[j]) for j in range(m+1) if Ms[j]!=0}
            vr={r:r+2*v2(math.comb(m,r)*R[r]) for r in range(m+1) if R[r]!=0}
            if not vj or not vr: continue
            muj=min(vj.values()); J=sorted(k for k in vj if vj[k]==muj)
            mur=min(vr.values()); Rl=sorted(k for k in vr if vr[k]==mur)
            if len(J)<2: continue
            tot+=1
            ok,_=is_affine_box(Rl)
            if ok: rbox+=1
            vp,_=vpi_G(lam,m)
            if vp==muj: jsharp+=1
            if vp==mur: rsharp+=1
            if J==Rl: jeqr+=1
    print(f"\nties(m<12): {tot}")
    print(f"  r-picture min-locus is affine box: {rbox}/{tot}")
    print(f"  j-min-locus == r-min-locus:        {jeqr}/{tot}")
    print(f"  j-min == true vpiG (j sharp):      {jsharp}/{tot}")
    print(f"  r-min == true vpiG (r sharp):      {rsharp}/{tot}")

if __name__=="__main__":
    main()
