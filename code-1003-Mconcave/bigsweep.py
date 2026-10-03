"""Extended census (brief hazard 2: the 135,422 are the discovery set).

Reduced form, verified exactly in l1form.py / realbox.py:
  y_i = nu_i - (lam_{i-1}+1),  effective slice = {y in prod[P_i,Q_i] : sum y_i = sigma},
  P_i = max(mu_i, lam_{i-2}+2)-lam_{i-1}-1,  Q_i = min(mu_{i+1}-1,lam_i)-lam_{i-1}-1,
  sigma = n-m-d+b,   k(y) = sum_i max(-y_i,0),   Lambda = (d-b)/2 - k.
M(r) concave on its support  <=>  N(k)=#{y : k(y)=k} concave on its support.

The box does not depend on b, so compute F(x,q) = prod_i sum_{t=P_i}^{Q_i} x^t q^{max(-t,0)}
ONCE per (mu,lam) and read off every slice as [x^sigma]F.
"""
import gen, tent, sys
from collections import Counter

def ybox(mu, lam, n, m):
    P=[];Q=[]
    for i in range(m):
        a=tent.lam_at(lam,i-1,n,m)+1
        P.append(max(mu[i], tent.lam_at(lam,i-2,n,m)+2)-a)
        Q.append(min(gen.nu_next(mu,i,n,m)-1, lam[i])-a)
    return P,Q

def profiles(P,Q):
    """returns {sigma: [N(k) for k=0..]} (trailing/leading zeros kept in k-index 0..)"""
    m=len(P)
    cur={0:{0:1}}          # x-deg -> {q-deg: count}
    for i in range(m):
        new={}
        for xd,qd in cur.items():
            for t in range(P[i],Q[i]+1):
                dq=max(-t,0); tgt=new.setdefault(xd+t,{})
                for q,c in qd.items(): tgt[q+dq]=tgt.get(q+dq,0)+c
        cur=new
    return cur

def is_cc(N):
    """N: list over k=kmin..kmax with no leading/trailing zeros. concave on support?"""
    nz=[i for i,c in enumerate(N) if c]
    if not nz: return True
    if nz[-1]-nz[0]+1!=len(nz): return False   # gap -> would also refute tent(3)
    s=N[nz[0]:nz[-1]+1]
    inc=[s[i+1]-s[i] for i in range(len(s)-1)]
    return all(inc[i]>=inc[i+1] for i in range(len(inc)-1))

def sweep(m,nmax,dmax,label):
    st=Counter(); wit=[]; lens=Counter()
    for n in range(m,nmax+1):
        for (mu,lam) in gen.pairs(n,m,dmax):
            d=sum(lam)-sum(mu); P,Q=ybox(mu,lam,n,m)
            if any(Q[i]<P[i] for i in range(m)): continue
            F=profiles(P,Q)
            for b in range(0,d+1):
                sigma=n-m-d+b
                if sigma not in F: continue
                qd=F[sigma]
                kmin,kmax=min(qd),max(qd)
                N=[qd.get(k,0) for k in range(kmin,kmax+1)]
                st['slices']+=1; lens[len(N)]+=1
                if is_cc(N): st['ok']+=1
                else:
                    st['FAIL']+=1
                    if len(wit)<8: wit.append((m,n,mu,lam,b,P,Q,sigma,N))
    print(f"[{label}] m={m} n<={nmax} d<={dmax}: {dict(st)}")
    print(f"   profile length distribution: {dict(sorted(lens.items()))}")
    for w in wit: print("   FAIL:",w)
    sys.stdout.flush()
    return st

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(); p.add_argument('m',type=int); p.add_argument('nmax',type=int); p.add_argument('dmax',type=int)
    a=p.parse_args(); sweep(a.m,a.nmax,a.dmax,'sweep')
