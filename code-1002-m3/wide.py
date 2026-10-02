"""RANGE EXTENSION for the strengthening 'beta log-concave', deliberately far outside
the range that bore it (m<=5, n<=10, d<=12).  Also (A) itself.
a-strengthening-is-confirmed-by-the-sweep-that-bore-it: a perfect score on the discovery
set is the tell, not the evidence."""
import gen, layers, sys
from collections import Counter

def scan(m,nmax,dmax):
    st=Counter(); wit=[]
    for n in range(m,nmax+1):
        for (mu,lam) in gen.pairs(n,m,dmax):
            for b in range(0,sum(lam)-sum(mu)+1):
                off,co=gen.slice_sum(mu,lam,n,m,b)
                if not co: continue
                st['slices']+=1
                if co!=co[::-1]: st['NOTSYM']+=1; continue
                bet=layers.beta_of(layers.radial(co))
                if min(bet)<0: st['BETANEG']+=1; continue
                st['A_ok' if gen.is_pf2(co) else 'A_FAIL']+=1
                st['beta_lc' if layers.is_lc(bet) else 'beta_NOTLC']+=1
                if not (gen.is_pf2(co) and layers.is_lc(bet)) and len(wit)<8:
                    wit.append((n,mu,lam,b,co,bet))
                if layers.is_concave_pospart(bet): st['beta_cplus']+=1
        print(f"   [m={m} n={n} done] {dict(st)}", flush=True)
    return st,wit

for (m,nmax,dmax) in [(3,16,20),(4,13,16),(5,12,13),(6,11,12),(7,11,11)]:
    st,wit=scan(m,nmax,dmax)
    print(f"=== m={m}, {m}<=n<={nmax}, d<={dmax}: {dict(st)}")
    for x in wit: print("   WITNESS",x)
    sys.stdout.flush()
