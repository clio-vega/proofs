"""The half-width multiplicity profile M(r) = #{nu in effective slice : Lambda(nu)=r}.
The smallest abstract counterexample to 'gap-free suffices' has M=(3,1,1), which is NOT
log-concave.  Question: what does M look like on REAL cylindric slices?"""
import gen, layers, tent
from collections import Counter
st=Counter(); wit=[]
for m in range(2,6):
    for n in range(m, m+6):
        for (mu,lam) in gen.pairs(n,m,11):
            for b in range(0,sum(lam)-sum(mu)+1):
                nus=[nu for nu in gen.slice_nus(mu,lam,n,m,b) if gen.f_nu(nu,lam,n,m)[1]]
                if not nus: continue
                hs=[tent.Lambda_formula(nu,lam,n,m) for nu in nus]
                lo,hi=min(hs),max(hs)
                M=[sum(1 for h in hs if abs(h-(lo+k))<1e-9) for k in range(int(hi-lo)+1)]
                st[f'm{m}_slices']+=1
                st[f'm{m}_M_lc' if layers.is_lc(M) else f'm{m}_M_NOTLC']+=1
                st[f'm{m}_M_concave' if layers.is_concave_pospart(M) else f'm{m}_M_NOTconcave']+=1
                st[f'm{m}_M_nondecr' if all(M[i]<=M[i+1] for i in range(len(M)-1)) else f'm{m}_M_not_nondecr']+=1
                st[f'm{m}_M_nonincr' if all(M[i]>=M[i+1] for i in range(len(M)-1)) else f'm{m}_M_not_nonincr']+=1
                if not layers.is_lc(M) and len(wit)<5: wit.append((m,n,mu,lam,b,M))
for k in sorted(st): print(f"{k:24s} {st[k]}")
print("\nM not log-concave e.g.:")
for w in wit: print("  ",w)
