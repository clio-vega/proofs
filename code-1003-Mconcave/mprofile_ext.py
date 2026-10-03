"""mprofile.py from 10-02, VERBATIM logic, with only the n-range widened.
10-02 ran n in range(m, m+6) -> at m=2 that is n<=7.  Here n goes further.
Nothing else is changed: same gen, same tent.Lambda_formula, same
layers.is_concave_pospart, same M construction."""
import gen, layers, tent, sys
from collections import Counter
st=Counter(); wit=[]
NEXTRA=int(sys.argv[1]) if len(sys.argv)>1 else 6
DMAX=int(sys.argv[2]) if len(sys.argv)>2 else 11
MS=[int(x) for x in sys.argv[3].split(',')] if len(sys.argv)>3 else [2,3,4,5]
for m in MS:
    for n in range(m, m+NEXTRA):
        for (mu,lam) in gen.pairs(n,m,DMAX):
            for b in range(0,sum(lam)-sum(mu)+1):
                nus=[nu for nu in gen.slice_nus(mu,lam,n,m,b) if gen.f_nu(nu,lam,n,m)[1]]
                if not nus: continue
                hs=[tent.Lambda_formula(nu,lam,n,m) for nu in nus]
                lo,hi=min(hs),max(hs)
                M=[sum(1 for h in hs if abs(h-(lo+k))<1e-9) for k in range(int(hi-lo)+1)]
                st[f'm{m}_slices']+=1
                st[f'm{m}_M_concave' if layers.is_concave_pospart(M) else f'm{m}_M_NOTconcave']+=1
                st[f'm{m}_M_lc' if layers.is_lc(M) else f'm{m}_M_NOTLC']+=1
                if not layers.is_concave_pospart(M) and len(wit)<6: wit.append((m,n,mu,lam,b,M))
for k in sorted(st): print(f"{k:24s} {st[k]}")
print("M NOT positive-part-of-concave, e.g.:")
for w in wit: print("  ",w)
