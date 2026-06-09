"""
JOB 3 — the surviving higher-order coefficient.

For each tie, the leading pi^mu coefficient cancels (|J*| even).  We:
  - record d := v_pi(G) - mu  (order of the surviving cancellation), confirm finite;
  - compute L0 = sum_{j in J*} (unit_j)  in Z[i], where term_j/(1+i)^mu reduces, for
    j in J*, to unit_j = u_j (-i)^{v_j} i^j  (u_j = odd part of C(m,j)M_j);
    its v_pi tells how deep the *within-J* * cancellation goes;
  - decide whether the survivor is WITHIN J* (v_pi(G)-mu == v_pi(L0), all J* terms
    carry the full surviving order) or needs NEXT-val indices (v_pi(G)-mu > v_pi(L0)).
"""
import math, csv
from job1_tie_census import partitions, M_vector, v2, G_gaussian, vpi
from job2_mechanism import val_pieces, jstar_of

def odd_part(n):
    while n%2==0: n//=2
    return n

def unit_sum_Jstar(lam,m,Ms,J):
    """L0 = sum_{j in J*} u_j (-i)^{v_j} i^j  in Z[i] (re,im)."""
    re,im=0,0
    for j in J:
        c=math.comb(m,j)*Ms[j]
        vj=v2(c); u=odd_part(c)
        # (-i)^{v_j}
        a,b=[(1,0),(0,-1),(-1,0),(0,1)][vj%4]
        # i^j
        ar,br=[(1,0),(0,1),(-1,0),(0,-1)][j%4]
        # (a+bi)(ar+br i)
        tr=a*ar-b*br; ti=a*br+b*ar
        re+=u*tr; im+=u*ti
    return re,im

if __name__=="__main__":
    rows=[]
    from collections import Counter
    ddist=Counter(); within=0; needs_next=0; tot=0
    for m in range(1,13):
        for lam in partitions(2*m):
            Ms,pj,pbin,pM,valv=val_pieces(lam,m)
            J,mu=jstar_of(valv)
            if len(J)<2: continue
            tot+=1
            re,im=G_gaussian(lam,m,Ms)
            vp=vpi(re,im)  # None if G=0
            d = None if vp is None else vp-mu
            L0re,L0im=unit_sum_Jstar(lam,m,Ms,J)
            vL0=vpi(L0re,L0im)   # >=1 since |J*| even; None if J*-sum itself 0
            ddist[('inf' if d is None else d)]+=1
            # classify
            if vp is None:
                cls='G=0'
            elif vL0 is None:
                cls='J*-sum=0(needs-next)'; needs_next+=1
            elif d==vL0:
                cls='within-J*'; within+=1
            else:
                cls='needs-next'; needs_next+=1
            rows.append(dict(m=m,lam='|'.join(map(str,lam)),Jsize=len(J),
                mu=mu, vpi=('inf' if vp is None else vp),
                d=('inf' if d is None else d),
                vL0=('inf' if vL0 is None else vL0), cls=cls))
    print(f"ties (m<=12): {tot}")
    print("distribution of d = v_pi(G)-mu :", dict(sorted(ddist.items(),key=lambda kv:(str(kv[0])))))
    print(f"survivor WITHIN J* (d==v_pi(L0)): {within}")
    print(f"survivor NEEDS next-val indices : {needs_next}")
    with open('results/job3-survivor.csv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['m','lam','Jsize','mu','vpi','d','vL0','cls'])
        w.writeheader(); w.writerows(rows)
    print("wrote results/job3-survivor.csv")
