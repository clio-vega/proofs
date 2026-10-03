"""What do the REAL boxes look like?  In y_i = nu_i - a_i, a_i = lam_{i-1}+1:
   P_i = max(mu_i, lam_{i-2}+2) - lam_{i-1} - 1,   Q_i = min(mu_{i+1}-1, lam_i) - lam_{i-1} - 1
   sigma = sum y_i on the slice = |mu|+b - sum a_i = n - m - d + b.
   k(y) = sum_i max(-y_i,0),  Lambda = c - k with c=(d-b)/2.
Collect the observed (P,Q,sigma) and test candidate structural facts."""
import gen, tent
from collections import Counter
from itertools import product

def ybox(mu, lam, n, m):
    P=[]; Q=[]
    for i in range(m):
        a = tent.lam_at(lam, i-1, n, m) + 1
        P.append(max(mu[i], tent.lam_at(lam, i-2, n, m) + 2) - a)
        Q.append(min(gen.nu_next(mu, i, n, m) - 1, lam[i]) - a)
    return P, Q

st=Counter(); wit=Counter(); ex={}
for m in range(2,6):
    for n in range(m, m+5):
        for (mu,lam) in gen.pairs(n,m,8):
            d = sum(lam)-sum(mu)
            P,Q = ybox(mu,lam,n,m)
            if any(Q[i]<P[i] for i in range(m)): st['emptybox']+=1; continue
            # candidate facts about the box itself
            st['box']+=1
            for i in range(m):
                st['P_le_0' if P[i]<=0 else 'P_pos']+=1
                st['Q_ge_0' if Q[i]>=0 else 'Q_neg']+=1
            # sum constraints
            if sum(Q) == n-m: st['sumQ_eq_n-m']+=1
            else:
                st['sumQ_NE']+=1; ex.setdefault('sumQ',(m,n,mu,lam,P,Q,sum(Q),n-m))
            if sum(P) <= 0: st['sumP_le0']+=1
            for b in range(0,d+1):
                sigma = n-m-d+b
                nus=[nu for nu in gen.slice_nus(mu,lam,n,m,b) if gen.f_nu(nu,lam,n,m)[1]]
                if not nus: continue
                st['slice']+=1
                ys=[tuple(nu[i]-(tent.lam_at(lam,i-1,n,m)+1) for i in range(m)) for nu in nus]
                if all(sum(y)==sigma for y in ys): st['sigma_ok']+=1
                else: st['sigma_BAD']+=1; ex.setdefault('sigma',(m,n,mu,lam,b,ys[:3],sigma))
                # does the y-box exactly describe the effective slice?
                full=[y for y in product(*[range(P[i],Q[i]+1) for i in range(m)]) if sum(y)==sigma]
                if set(full)==set(ys): st['box_exact']+=1
                else: st['box_BAD']+=1; ex.setdefault('boxex',(m,n,mu,lam,b,P,Q,sigma,sorted(set(full)^set(ys))[:4]))
print(dict(st))
for k,v in ex.items(): print(' EX',k,v)
