"""(a) The 162 slices at m=3 where beta is NOT the positive part of a concave function:
    minimise and record WITH incidental properties.
(b) Wide sweep of the strengthening 'beta log-concave', far outside its birth range."""
import gen, layers, sys
from collections import Counter

def scan(m, nmax, dmax, collect_nonconcave=True):
    st=Counter(); nc=[]; nlc=[]
    for n in range(m,nmax+1):
        for (mu,lam) in gen.pairs(n,m,dmax):
            d=sum(lam)-sum(mu)
            for b in range(0,d+1):
                off,co=gen.slice_sum(mu,lam,n,m,b)
                if not co: continue
                st['slices']+=1
                assert co==co[::-1], (n,mu,lam,b,co)
                bet=layers.beta_of(layers.radial(co))
                assert min(bet)>=0
                if gen.is_pf2(co): st['A_ok']+=1
                else:
                    st['A_FAIL']+=1
                    nlc.append(('A',n,mu,lam,b,co))
                if layers.is_lc(bet): st['beta_lc']+=1
                else:
                    st['beta_NOTLC']+=1
                    nlc.append(('beta',n,mu,lam,b,co,bet))
                if layers.is_concave_pospart(bet): st['beta_cplus']+=1
                else:
                    st['beta_not_cplus']+=1
                    if collect_nonconcave: nc.append((d,n,mu,lam,b,co,bet))
    return st,nc,nlc

st,nc,nlc = scan(3,10,12)
print("m=3 n<=10 d<=12:", dict(st))
nc.sort()
print(f"\n{len(nc)} slices where beta is NOT c_+ for concave c.  Smallest by d:")
for x in nc[:6]:
    d,n,mu,lam,b,co,bet=x
    print(f"  d={d} n={n} mu={mu} lam={lam} b={b}")
    print(f"     G={co}   beta={bet}   LC={layers.is_lc(bet)}  G PF2={gen.is_pf2(co)}")
    nus=gen.slice_nus(mu,lam,n,3,b)
    for nu in nus:
        lr=gen.LR(nu,lam,n,3); w=[R-L+1 for (L,R) in lr]
        print(f"       nu={nu} LR={lr} w={w} f_nu={gen.f_nu(nu,lam,n,3)}")
print("\nany (A) or beta-LC failures:", nlc[:3])
