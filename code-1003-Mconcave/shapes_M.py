"""What M profiles actually occur on real slices?  And at m=2, what are the arm
lengths L,R and the flat count N(0)?  (At m=2 the abstract analysis says
N = (N(0),2,...,2,1,...,1) with min(L,R) twos and |L-R| ones, so concavity
needs |L-R|<=1 AND N(0)<=2.  Is that what real slices do?)"""
import gen, tent
from collections import Counter
from itertools import product

def ybox(mu, lam, n, m):
    P=[];Q=[]
    for i in range(m):
        a=tent.lam_at(lam,i-1,n,m)+1
        P.append(max(mu[i], tent.lam_at(lam,i-2,n,m)+2)-a)
        Q.append(min(gen.nu_next(mu,i,n,m)-1, lam[i])-a)
    return P,Q

prof=Counter(); st=Counter(); arms=Counter(); flat=Counter()
for m in range(2,6):
    for n in range(m, m+6):
        for (mu,lam) in gen.pairs(n,m,9):
            d=sum(lam)-sum(mu); P,Q=ybox(mu,lam,n,m)
            if any(Q[i]<P[i] for i in range(m)): continue
            for b in range(0,d+1):
                sigma=n-m-d+b
                pts=[y for y in product(*[range(P[i],Q[i]+1) for i in range(m)]) if sum(y)==sigma]
                if not pts: continue
                c=Counter(sum(max(-t,0) for t in y) for y in pts)
                kmin,kmax=min(c),max(c)
                N=[c.get(k,0) for k in range(kmin,kmax+1)]
                st[f'm{m}']+=1
                prof[(m,tuple(N))]+=1
                if kmin!=0: st[f'm{m}_kmin_nonzero']+=1
                if m==2:
                    flat[N[0] if kmin==0 else 'kmin>0']+=1
                    # arm lengths
                    ts=sorted(y[0] for y in pts)
                    A,Z=ts[0],ts[-1]
                    L=max(0,-A); R=max(0,Z-sigma)
                    arms[abs(L-R)]+=1
print("slices:",dict(st))
print("\nm=2: |L-R| distribution:",dict(arms))
print("m=2: N(0) distribution:",dict(flat))
print("\nmost common profiles per m:")
for m in range(2,6):
    sub=[(v,p) for (mm,p),v in prof.items() if mm==m]
    sub.sort(reverse=True)
    print(f" m={m}: ", sub[:8])
    print(f"   longest: ", sorted(set(p for v,p in sub), key=len)[-4:])
