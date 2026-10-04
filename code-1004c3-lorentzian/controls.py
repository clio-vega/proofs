"""Controls C1-C3 + regression gate + the M-concavity hedge, for the 1004c3 Lorentzian
route to (Q).  Every block prints the NUMBER OF OBJECTS ENUMERATED and the SIZE of the
objects, not only the number of failures."""
import itertools, random, sys
from lor import *

random.seed(31004)

# =====================================================================  REGRESSION
print("=== R. REGRESSION GATE (one pass; this is NOT the session's evidence) ===")
n=0; bad=0; sz=[]
for h in range(0,11):
    for l in range(0,h+1):
        c=Ptilde(l,h); n+=1; sz.append(len(c))
        if not is_N_lorentzian(c): bad+=1
print(f"  N(Ptilde) Lorentzian: pairs 0<=l<=h<=10 enumerated {n}, failures {bad}; "
      f"support sizes {min(sz)}..{max(sz)}")
n=0; bad=0; sz=[]
for h1,h2 in itertools.product(range(0,6),repeat=2):
    for l1 in range(0,h1+1):
        for l2 in range(0,h2+1):
            f=poly_mult(Ptilde(l1,h1),Ptilde(l2,h2)); n+=1; sz.append(len(f))
            if not is_N_lorentzian(f): bad+=1
print(f"  m=2 products: enumerated {n}, failures {bad}; support sizes {min(sz)}..{max(sz)}")
n=0; bad=0; sz=[]
for hs in itertools.product(range(0,4),repeat=3):
    for ls in itertools.product(*[range(0,h+1) for h in hs]):
        f=Ptilde(ls[0],hs[0])
        for i in (1,2): f=poly_mult(f,Ptilde(ls[i],hs[i]))
        n+=1; sz.append(len(f))
        if not is_N_lorentzian(f): bad+=1
print(f"  m=3 products: enumerated {n}, failures {bad}; support sizes {min(sz)}..{max(sz)}")

# =====  cross-check of Prop 4 (product form) against an INDEPENDENT enumeration
print("=== P4. product form  k(a) = [X^a E^{D-a} W^{H-D}] prod_i Ptilde_i ===")
def G_family(m,l,h,cross=None):
    """3-variable homogenisation of ANY constraint family, graded by (sum x, sum e, rest)."""
    H=sum(h); c={}
    for x in itertools.product(*[range(0,h[i]+1) for i in range(m)]):
        for e in itertools.product(*[range(0,h[i]+1) for i in range(m)]):
            if any(not(l[i]<=x[i]+e[i]<=h[i]) for i in range(m)): continue
            if cross is not None:
                clo,chi=cross
                if any(not(clo<=e[i]+x[(i+1)%m]<=chi) for i in range(m)): continue
            a,b=sum(x),sum(e)
            key=(a,b,H-a-b); c[key]=c.get(key,0)+1
    return c
tot=ok=0
for m in (1,2,3):
    for h in itertools.product(range(0,4),repeat=m):
        for l in itertools.product(*[range(0,hi+1) for hi in h]):
            f=({(0,0,0):1})
            for i in range(m): f=poly_mult(f,Ptilde(l[i],h[i]))
            g=G_family(m,list(l),list(h))
            tot+=1; ok+= (f==g)
            if f!=g and ok<2: print("   MISMATCH",l,h)
print(f"  prod_i Ptilde_i vs direct graded enumeration: instances {tot}, agreements {ok}")
tot=ok=0; Ds=[]
for m in (1,2,3):
    for h in itertools.product(range(1,4),repeat=m):
        for l in itertools.product(*[range(0,hi+1) for hi in h]):
            f=({(0,0,0):1})
            for i in range(m): f=poly_mult(f,Ptilde(l[i],h[i]))
            H=sum(h)
            for D in range(sum(l),H+1):
                ks=k_sequence(list(l),list(h),D)
                pred=[f.get((a,D-a,H-D),0) for a in range(0,D+1)]
                tot+=1; ok+=(ks==pred); Ds.append(D)
                if ks!=pred and ok<2: print("   MISMATCH",l,h,D,ks,pred)
print(f"  k(a) vs the extracted coefficient: instances {tot}, agreements {ok}; D range {min(Ds)}..{max(Ds)}")

# =====================================================================  C1
print("=== C1. delete ONE interior support point of Ptilde -- chain must REFUSE ===")
tot=ref=0; detail=[]
for (l,h) in [(0,3),(1,4),(0,4),(2,5),(1,5)]:
    c=Ptilde(l,h)
    pts=sorted(c)
    interior=[p for p in pts if all(p[i]>0 for i in range(3))] or pts[1:-1]
    for p in interior:
        c2=dict(c); del c2[p]
        tot+=1
        okL, w = is_N_lorentzian(c2, return_witness=True)
        if not okL:
            ref+=1
            if len(detail)<3: detail.append((l,h,p,str(w)[:60]))
print(f"  deletions enumerated {tot}, refusals {ref}  (must equal {tot})")
for d in detail: print("    e.g. (l,h)=",d[0],d[1]," deleted",d[2]," reason:",d[3])

# =====================================================================  C2
print("=== C2. NON-LAMINAR crossing family: chain must FAIL TO CERTIFY, and the ===")
print("===     conclusion must actually be FALSE somewhere (hypothesis load-bearing) ===")
def k_cross(m,l,h,D,cross):
    out={}
    for x in itertools.product(*[range(0,h[i]+1) for i in range(m)]):
        for e in itertools.product(*[range(0,h[i]+1) for i in range(m)]):
            if sum(x)+sum(e)!=D: continue
            if any(not(l[i]<=x[i]+e[i]<=h[i]) for i in range(m)): continue
            clo,chi=cross
            if any(not(clo<=e[i]+x[(i+1)%m]<=chi) for i in range(m)): continue
            out[sum(x)]=out.get(sum(x),0)+1
    return [out.get(a,0) for a in range(0,D+1)] if out else []
for cross in [(0,1),(1,2),(0,2)]:
    tot=notlc=nMconv=notLor=0; wit=None
    for m in (2,3):
        for h in itertools.product(range(1,4),repeat=m):
            for l in itertools.product(*[range(0,hi+1) for hi in h]):
                g=G_family(m,list(l),list(h),cross)
                if not g: continue
                nMconv += 0 if is_Mconvex([a for a,v in g.items() if v]) else 1
                notLor += 0 if is_N_lorentzian(g) else 1
                tot+=1
                for D in range(sum(l),sum(h)+1):
                    ks=k_cross(m,list(l),list(h),D,cross)
                    if ks and not is_logconcave_pf2(ks):
                        notlc+=1
                        if wit is None: wit=(m,l,h,D,ks)
    print(f"  cross={cross}: families enumerated {tot}; support NOT M-convex in {nMconv}; "
          f"N(G) NOT Lorentzian in {notLor}; counts NOT PF2 in {notlc} (D,instance pairs)")
    if wit: print(f"    witness that the conclusion is FALSE there: m={wit[0]} l={wit[1]} h={wit[2]} D={wit[3]} k={wit[4]}")

# =====================================================================  C3
print("=== C3. all-1 coefficients on a NON-M-convex support: normalizedcoefficients ===")
print("===     must not apply, i.e. N(f) must not be certified Lorentzian ===")
cands=[ [(2,0,0),(0,2,0),(0,0,2)],
        [(2,0,0),(0,2,0),(1,1,0),(0,0,2)],
        [(3,0,0),(0,3,0),(0,0,3),(1,1,1)],
        [(2,1,0),(0,1,2),(1,1,1)],
        [(2,0,1),(1,2,0),(0,1,2)] ]
for S in cands:
    c={a:1 for a in S}
    print(f"  support {S}: M-convex={is_Mconvex(S)}, N(f) Lorentzian={is_N_lorentzian(c)}, "
          f"ignoring support={is_N_lorentzian(c,check_support=False)}")

# ============================================================  HEDGE
print("=== HEDGE. is nu = log c M-concave for the PRODUCT f = prod_i Ptilde_i? ===")
print("===        if yes, normalizedcoefficients alone suffices, dropping CorollaryConvolution ===")
def is_M_concave(c):
    supp=[a for a,v in c.items() if v!=0]; S=set(supp)
    if not is_Mconvex(supp): return False
    n=len(supp[0])
    for a in supp:
        for b in supp:
            for i in range(n):
                if a[i]>b[i]:
                    ok=False
                    for j in range(n):
                        if a[j]<b[j]:
                            a2=list(a);a2[i]-=1;a2[j]+=1;a2=tuple(a2)
                            b2=list(b);b2[i]+=1;b2[j]-=1;b2=tuple(b2)
                            if a2 in S and b2 in S and c[a]*c[b]<=c[a2]*c[b2]: ok=True;break
                    if not ok: return False
    return True
tot=yes=0; firstno=None
for m in (1,2,3):
    for h in itertools.product(range(1,5),repeat=m):
        for l in itertools.product(*[range(0,hi+1) for hi in h]):
            f={(0,0,0):1}
            for i in range(m): f=poly_mult(f,Ptilde(l[i],h[i]))
            tot+=1
            if is_M_concave(f): yes+=1
            elif firstno is None: firstno=(m,l,h)
print(f"  products enumerated {tot}; nu M-concave in {yes}; first failure {firstno}")
