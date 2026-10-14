"""Which (y,sigma,S) attain the top t-exponent at k>=2, and does the top coefficient vanish?"""
import sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1014-q458')
from itertools import product, permutations
from collections import defaultdict
from censusk import perm_sign
from walk import walk_counts

def analyse(n,k,ell,show=3):
    M=k+ell+1; N=2*M
    W=walk_counts(n,k)
    def Wn(p): return W.get(tuple(p),0)
    terms=defaultdict(int)
    detail=defaultdict(list)
    ymax=(2*n+2*k+4)//N+1
    for ys in product(range(-ymax,ymax+1),repeat=k):
        for sg in permutations(range(1,k+1)):
            s=perm_sign(sg)
            E=sum(M*ys[i]**2 - sg[i]*ys[i] for i in range(k))
            tot=0
            for S in product((0,1),repeat=k):
                p=[sg[i]+(1 if S[i] else -1)*(i+1)-N*ys[i] for i in range(k)]
                w=Wn(p)
                if w: tot += (-1)**sum(S)*w
            if tot:
                terms[E]+=s*tot
                detail[E].append((ys,sg,s*tot))
    terms={e:v for e,v in terms.items() if v}
    if not terms: return None
    D=max(terms)
    # also: what is the max exponent over all (y,sigma) with ANY nonzero walk support,
    # i.e. before the alternating sum over S is taken?
    return D, terms[D], detail[D], len(terms)

if __name__=='__main__':
    for k,ell in [(2,0),(2,1),(3,0),(3,1),(4,0)]:
        M=k+ell+1
        print(f"### k={k} ell={ell} M={M}",flush=True)
        for n in range(2, (20 if k<4 else 14)+1, 2):
            r=analyse(n,k,ell)
            if r is None: print(f"  n={n}: c=0",flush=True); continue
            D,coef,det,nt=r
            ys_list=sorted(set(d[0] for d in det))
            print(f"  n={n:3d}: deg={D:4d} topcoef={coef:<12d} nterms={nt:3d} "
                  f"#(y,sigma) at top={len(det)} y's={ys_list[:4]}",flush=True)
