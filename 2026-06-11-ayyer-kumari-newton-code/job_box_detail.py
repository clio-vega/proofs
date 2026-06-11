"""
Detail on B* box structure, esp. H1 failures.  For each shape print
chi_k, v2(chi_k), v2(C(m,k)), beta_k, and B*.  Also test refined hypotheses:
  H2: B* is the argmin of beta_k and is always an affine 2-adic box.
  H3: v2(chi_k) itself, as fn of k, is governed by binary digits.
"""
import math
from job1_tie_census import partitions, v2, chi_b_vector

def submasks(m):
    out=[]; sub=m
    while True:
        out.append(sub)
        if sub==0: break
        sub=(sub-1)&m
    return sorted(out)

def is_affine_box(S):
    S=sorted(S); n=len(S)
    if n==0: return False,None
    if n&(n-1)!=0: return False,None
    j0=S[0]
    gens=[s-j0 for s in S if (s-j0)>0 and ((s-j0)&(s-j0-1))==0]
    from itertools import combinations
    box=set()
    for r in range(len(gens)+1):
        for c in combinations(gens,r):
            box.add(j0+sum(c))
    return (box==set(S)), gens

def detail(lam,m):
    chi=chi_b_vector(lam,m)
    beta={k:(v2(math.comb(m,k)*chi[k]) if chi[k]!=0 else None) for k in range(m+1)}
    fin={k:v for k,v in beta.items() if v is not None}
    bmin=min(fin.values()); B=sorted(k for k in fin if fin[k]==bmin)
    ok,gens=is_affine_box(B)
    print(f"lam={lam} m={m}={m:b}b  B*={B} box={ok} gens={gens}  submasks(m)={submasks(m)}")
    print(f"   k     : "+" ".join(f"{k:3d}" for k in range(m+1)))
    print(f"   chi_k : "+" ".join(f"{chi[k]:3d}" for k in range(m+1)))
    print(f"   v2chi : "+" ".join(f"{(v2(chi[k]) if chi[k]!=0 else 99):3d}" for k in range(m+1)))
    print(f"   v2Cmk : "+" ".join(f"{v2(math.comb(m,k)):3d}" for k in range(m+1)))
    print(f"   beta  : "+" ".join(f"{(beta[k] if beta[k] is not None else 99):3d}" for k in range(m+1)))
    print()

if __name__=="__main__":
    # H1 failures
    for lam,m in [((5,2,1),4),((7,2,1),5),((9,2,1),6),((5,4,3),6),((11,2,1),7),
                  ((6,3,2,1),6),((5,2,1,1,1),5)]:
        detail(lam,m)
    # verify: is B* ALWAYS an affine 2-adic box over all shapes m<=12?
    tot=0; boxok=0; bad=[]
    for m in range(1,13):
        for lam in partitions(2*m):
            chi=chi_b_vector(lam,m)
            fin={k:v2(math.comb(m,k)*chi[k]) for k in range(m+1) if chi[k]!=0}
            if not fin: continue
            bmin=min(fin.values()); B=sorted(k for k in fin if fin[k]==bmin)
            tot+=1
            ok,_=is_affine_box(B)
            if ok: boxok+=1
            else: bad.append((m,lam,B))
    print(f"B* is an affine 2-adic box: {boxok}/{tot}")
    for x in bad[:20]: print("  NOT BOX:",x)
