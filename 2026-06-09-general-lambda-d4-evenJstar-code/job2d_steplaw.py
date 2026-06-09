"""
JOB 2d — the crisp lever for PROVE: the valuation-step law on involution pairs.

If J* is the affine box j0 + {subset sums of 2^a, a in S}, the smallest-generator
toggle sigma(j) = j +/- 2^{amin} is a fixed-point-free involution of J*.  For it
to preserve val(j)=j+2 v2(C(m,j) M_j), we need on each pair {j, j+2^a} (a in S):

      v2( C(m,j)   M_j )  -  v2( C(m,j+2^a) M_{j+2^a} )  =  2^{a-1}.

We verify this exact step law for EVERY generator a in S of EVERY tie, and also
test whether it already holds at the level of v2(M_j) alone (binomial-free) when
v2(C(m,.)) is constant across the pair.
"""
import math, csv
from job1_tie_census import partitions, M_vector, v2
from job2_mechanism import val_pieces, jstar_of, box_generators

if __name__=="__main__":
    tot=0; steplaw_ok=0; pairs=0; pairs_ok=0
    Mlaw_pairs=0; Mlaw_ok=0
    fails=[]
    for m in range(1,13):
        for lam in partitions(2*m):
            Ms,pj,pbin,pM,valv=val_pieces(lam,m)
            J,mu=jstar_of(valv)
            if len(J)<2: continue
            tot+=1
            j0,S=box_generators(J)
            Jset=set(J)
            this_ok=True
            for a in S:
                g=1<<a
                for j in J:
                    if ((j-j0)>>a)&1: continue  # pair from the 'off' side only
                    jp=j+g
                    if jp not in Jset: this_ok=False; continue
                    pairs+=1
                    vj=v2(math.comb(m,j)*Ms[j]); vjp=v2(math.comb(m,jp)*Ms[jp])
                    if vj-vjp==2**(a-1): pairs_ok+=1
                    else:
                        this_ok=False
                        fails.append((m,lam,j,jp,a,vj,vjp))
                    # binomial-free sub-test
                    if v2(math.comb(m,j))==v2(math.comb(m,jp)):
                        Mlaw_pairs+=1
                        if v2(Ms[j])-v2(Ms[jp])==2**(a-1): Mlaw_ok+=1
            if this_ok: steplaw_ok+=1
    print(f"ties (m<=12): {tot}")
    print(f"  shapes where step law holds for ALL generators: {steplaw_ok}/{tot}")
    print(f"  involution pairs satisfying  Dv2 = 2^(a-1):     {pairs_ok}/{pairs}")
    print(f"  binomial-constant pairs satisfying  Dv2(M)=2^(a-1): {Mlaw_ok}/{Mlaw_pairs}")
    if fails:
        print("  FAILURES:")
        for f in fails[:20]: print("   ",f)
    else:
        print("  NO FAILURES -- step law is exact.")
