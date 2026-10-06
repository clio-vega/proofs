"""The search LS Corollary 4's remark asks for, plus the planted controls.

A candidate object is a map f : Gamma(u,w) -> disjoint-union_v Gamma(w_0 v, w_0)
that is (i) TYPE-PRESERVING and (ii) has |f^{-1}(C)| = c^w_{u,v} for every v and
every C in Gamma(w_0 v, w_0).

Because f is type-preserving it decomposes over types, and within a type the
data is exactly a labelled set partition of Gamma_alpha(u,w) into blocks
B_C of size c^w_{u,v}.  So:
  - FEASIBLE  iff  I_alpha(u,w) = sum_v c^w_{u,v} I_alpha(w_0 v, w_0) for all alpha
                   (= Corollary 4), and
  - the number of such f is  prod_alpha  multinomial( I_alpha(u,w) ; block sizes ).
"""
import collections, math
from ls import *
from cor4 import cor4_data

def targets(al, cs, Itab):
    """list of (v, chain, required fibre size) for one type alpha"""
    out=[]
    for v,c in sorted(cs.items()):
        for ch in Itab[v].get(al,[]):
            out.append((v,ch,c))
    return out

def build(u,w,n,W0,Itab,C, cmod=None, alpha_override=None):
    """construct a map; return (ok, nmaps, detail).  cmod perturbs c."""
    lhs=I_table(u,w,n); cs=dict(C.get((u,w),{}))
    if cmod: cs=cmod(cs)
    alphas=set(lhs)|{a for v in cs for a in Itab[v]}
    nmaps=1; fmap={}
    for al in alphas:
        src=[l for l,_ in lhs.get(al,[])]
        tg = targets(alpha_override(al) if alpha_override else al, cs, Itab)
        need=sum(sz for _,_,sz in tg)
        if len(src)!=need:
            return False, 0, ('size mismatch', u,w,al,len(src),need)
        # greedy assignment: fill blocks in order
        i=0
        for (v,ch,sz) in tg:
            for _ in range(sz):
                fmap[(al,src[i])]=(v,ch); i+=1
        # multinomial count
        num=math.factorial(len(src)); den=1
        for (_,_,sz) in tg: den*=math.factorial(sz)
        nmaps *= num//den
    return True, nmaps, fmap

def verify_map(u,w,n,W0,Itab,C,fmap):
    """independently re-check the constructed map against BOTH constraints"""
    lhs=I_table(u,w,n); cs=C.get((u,w),{})
    # (i) totality and type-preservation
    for al,ch in lhs.items():
        for l,_ in ch:
            if (al,l) not in fmap: return 'not total'
            v,tgt=fmap[(al,l)]
            if chain_type(tgt[0],n)!=al: return 'type not preserved'
    # (ii) fibre sizes
    fib=collections.Counter()
    for key,(v,tgt) in fmap.items(): fib[(v,tgt[0])]+=1
    for v,c in cs.items():
        for al,ch in Itab[v].items():
            for l,_ in ch:
                if fib.get((v,l),0)!=c: return 'fibre size wrong for v=%s'%(v,)
    # no stray image
    for (v,l),k in fib.items():
        if v not in cs: return 'image outside the prescribed components'
    return 'OK'

if __name__=='__main__':
    for n in (3,4):
        W0,Itab,C = cor4_data(n)
        tot=feas=0; counts=collections.Counter(); verif=collections.Counter()
        for u in perms(n):
            for w in perms(n):
                if length(w)<length(u): continue
                tot+=1
                ok,nm,d = build(u,w,n,W0,Itab,C)
                if ok:
                    feas+=1; counts[nm]+=1
                    verif[verify_map(u,w,n,W0,Itab,C,d)]+=1
        print('n=%d: %d pairs (u,w) with u<=w in length; %d FEASIBLE' % (n,tot,feas))
        print('      independent re-verification of each constructed map:', dict(verif))
        print('      number of valid maps N(u,w), distribution:',
              dict(sorted(counts.items())))
        print('      max N(u,w) =', max(counts))
