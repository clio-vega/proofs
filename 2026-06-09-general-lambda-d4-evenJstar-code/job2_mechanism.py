"""
JOB 2 — the mechanism behind even-|J*|.

Empirical structure (Job 1): for every tie shape (lam != unique-min),
    J* = j0 + { sum_{a in A} 2^a : A subset S }   (an affine 2-adic box),
so |J*| = 2^{|S|} is even.  Here we locate WHICH piece of
    val(j) = j + 2*v2(C(m,j)) + 2*v2(M_j)
produces the box, hunt for the generating exponent set S, and test
a fixed-point-free involution j <-> j XOR 2^{a_min} on J*.
"""
import math, csv
from job1_tie_census import partitions, M_vector, v2, analyze

def val_pieces(lam, m):
    Ms = M_vector(lam, m)
    pj=[]; pbin=[]; pM=[]; valv=[]
    for j in range(m+1):
        if Ms[j]==0:
            pj.append(j); pbin.append(None); pM.append(None); valv.append(None); continue
        vb=v2(math.comb(m,j)); vM=v2(Ms[j])
        pj.append(j); pbin.append(vb); pM.append(vM); valv.append(j+2*vb+2*vM)
    return Ms,pj,pbin,pM,valv

def jstar_of(valv):
    fin=[v for v in valv if v is not None]; mu=min(fin)
    return [j for j,v in enumerate(valv) if v==mu], mu

def box_generators(J):
    """Given J (an affine pow2 box), return (j0, sorted exponent set S)."""
    j0=J[0]; D=sorted(j-j0 for j in J)
    gens=[d for d in D if d>0 and (d&(d-1))==0]
    from itertools import combinations
    k=len(J); nbits=k.bit_length()-1
    for combo in combinations(gens,nbits):
        ss=set()
        for r in range(len(combo)+1):
            for c in combinations(combo,r): ss.add(sum(c))
        if ss==set(D):
            return j0, sorted(int(g).bit_length()-1 for g in combo)
    return j0, None

if __name__=="__main__":
    # Q1: does the BINOMIAL part alone (w(j)=j+2 v2(C(m,j)), restricted to M_j!=0)
    #     already give a box? Compare argmin of w vs argmin of full val on ties.
    out=[]
    binbox=0; fullbox=0; bineq=0; nties=0
    invol_ok=0
    for m in range(1,13):
        for lam in partitions(2*m):
            Ms,pj,pbin,pM,valv=val_pieces(lam,m)
            J,mu=jstar_of(valv)
            if len(J)<2: continue
            nties+=1
            # binomial-only argmin among j with M_j != 0
            w=[ (j+2*pbin[j]) if pbin[j] is not None else None for j in range(m+1)]
            wf=[x for x in w if x is not None]; wmin=min(wf)
            Jw=[j for j in range(m+1) if w[j]==wmin]
            # full box generators
            j0,S=box_generators(J)
            # canonical involution: toggle the smallest generator 2^{min S} in/out
            # of the subset A (i.e. +g if its bit is absent, -g if present).
            # Fixed-point-free, val-invariant => |J*| even.
            inv_ok=True
            if S:
                g=1<<min(S); Jset=set(J)
                for j in J:
                    partner = j+g if ((j-j0)&g)==0 else j-g
                    if partner not in Jset or partner==j: inv_ok=False; break
            if inv_ok and S: invol_ok+=1
            out.append(dict(m=m,lam='|'.join(map(str,lam)),
                            Jstar='|'.join(map(str,J)),Jsize=len(J),
                            j0=j0,S='|'.join(map(str,S)) if S else 'NONE',
                            binJ='|'.join(map(str,Jw)),binsize=len(Jw),
                            box_in_bin=(len(Jw)==len(J) and Jw==J)))
            if len(Jw)==len(J) and Jw==J: bineq+=1
    print(f"ties (m<=12): {nties}")
    print(f"  full J* == binomial-only argmin J : {bineq}/{nties}")
    print(f"  involution j<->j XOR 2^min(S) maps J*->J* : {invol_ok}/{nties}")
    with open('results/job2-box-structure.csv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['m','lam','Jstar','Jsize','j0','S','binJ','binsize','box_in_bin'])
        w.writeheader(); w.writerows(out)
    print("wrote results/job2-box-structure.csv")
