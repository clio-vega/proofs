"""Verify the structural lemmas behind the refutation.

(S1)  k(nu) = sum_i (lam_{i-1}+1-nu_i)_+  and  k=0  <=>  lam/nu is a horizontal strip.
(S2)  With g_i = lam_i-lam_{i-1}-1 >= 0 (sum g_i = n-m =: G) and
      delta_i = mu_i - lam_{i-1} - 1,
          P_i = max(delta_i, -g_{i-1}),   Q_i = g_i - (delta_{i+1})_-,
          sigma = sum_i delta_i + b,      d = G - sum_i delta_i.
(S3)  sum_i (Q_i - P_i) = 2G - sum_i (delta_i)_- - sum_i (delta_i + g_{i-1})_+  <= 2G.
(S4)  m=2 structure:  k(t) = (sigma)_- + dist(t, I_0),  I_0=[min(0,sigma),max(0,sigma)],
      so N = (F, 2^J, 1^D) with F=#(window cap I_0), J=min(L,R), D=|L-R|.
(S5)  the exact concavity criterion for that shape.
"""
import gen, tent, cyl
from collections import Counter
from itertools import product
from layers import is_concave_pospart

def gd(mu,lam,n,m):
    g=[tent.lam_at(lam,i,n,m)-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
    dl=[mu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
    return g,dl

st=Counter(); ex={}
for m in range(2,6):
    for n in range(m,m+7):
        for (mu,lam) in gen.pairs(n,m,8):
            g,dl=gd(mu,lam,n,m); G=n-m
            st['sum_g_ok' if sum(g)==G else 'sum_g_BAD']+=1
            st['g_nonneg_ok' if all(x>=0 for x in g) else 'g_BAD']+=1
            d=sum(lam)-sum(mu)
            st['d_ok' if d==G-sum(dl) else 'd_BAD']+=1
            # (S2) box
            P=[max(dl[i], -g[(i-1)%m]) for i in range(m)]
            Q=[g[i]-max(-dl[(i+1)%m],0) for i in range(m)]
            P0=[max(mu[i], tent.lam_at(lam,i-2,n,m)+2)-(tent.lam_at(lam,i-1,n,m)+1) for i in range(m)]
            Q0=[min(gen.nu_next(mu,i,n,m)-1,lam[i])-(tent.lam_at(lam,i-1,n,m)+1) for i in range(m)]
            if P==P0 and Q==Q0: st['S2_ok']+=1
            else: st['S2_BAD']+=1; ex.setdefault('S2',(m,n,mu,lam,P,P0,Q,Q0))
            # (S3)
            lhs=sum(Q[i]-P[i] for i in range(m))
            rhs=2*G - sum(max(-dl[i],0) for i in range(m)) - sum(max(dl[i]+g[(i-1)%m],0) for i in range(m))
            if lhs==rhs and lhs<=2*G: st['S3_ok']+=1
            else: st['S3_BAD']+=1; ex.setdefault('S3',(m,n,mu,lam,lhs,rhs,2*G))
            # (S1) on every nu of every slice
            for b in range(0,d+1):
                sigma=sum(dl)+b
                nus=[nu for nu in gen.slice_nus(mu,lam,n,m,b) if gen.f_nu(nu,lam,n,m)[1]]
                if not nus: continue
                st['slice']+=1
                st['sigma_ok' if sigma==n-m-d+b else 'sigma_BAD']+=1
                c=(d-b)/2
                for nu in nus:
                    k=sum(max(tent.lam_at(lam,i-1,n,m)+1-nu[i],0) for i in range(m))
                    Lam=tent.Lambda_formula(nu,lam,n,m)
                    if abs(Lam-(c-k))<1e-9: st['S1a_ok']+=1
                    else: st['S1a_BAD']+=1; ex.setdefault('S1a',(m,n,mu,lam,b,nu,Lam,c,k))
                    hs = cyl.hstrip(nu,lam,n,m)
                    if (k==0)==bool(hs): st['S1b_ok']+=1
                    else: st['S1b_BAD']+=1; ex.setdefault('S1b',(m,n,mu,lam,nu,k,hs))
print(dict(st))
for k,v in ex.items(): print(' EX',k,v)

# ---------- (S4)+(S5): m=2 structure theorem ----------
print("\n(S4)/(S5) m=2 structure theorem, exhaustive over boxes+sigma")
def predict(P,Q,sigma):
    A=max(P[0],sigma-Q[1]); Z=min(Q[0],sigma-P[1])
    if Z<A: return None
    u,v=min(0,sigma),max(0,sigma)
    F=max(0,min(Z,v)-max(A,u)+1)
    L=max(0,min(Z,u-1)-A+1); R=max(0,Z-max(A,v+1)+1)
    if F==0: return [1]*(Z-A+1)
    J=min(L,R); D=abs(L-R)
    return [F]+[2]*J+[1]*D
def brute(P,Q,sigma):
    c=Counter()
    for t in range(P[0],Q[0]+1):
        s=sigma-t
        if not (P[1]<=s<=Q[1]): continue
        c[max(-t,0)+max(-s,0)]+=1
    if not c: return None
    return [c.get(k,0) for k in range(min(c),max(c)+1)]
def crit(F,J,D):
    """closed-form concavity criterion for N=(F,2^J,1^D), F>=1"""
    if J==0: return D<=1 or F<=1
    if D>=2: return False
    if J==1 and D==1: return F<=3
    if J==1 and D==0: return True
    return F<=2            # J>=2, D in {0,1}
bad=0; badc=0; tot=0
R0=range(-6,7)
for P in product(R0,repeat=2):
    for Q in product(R0,repeat=2):
        if Q[0]<P[0] or Q[1]<P[1]: continue
        for sigma in range(P[0]+P[1],Q[0]+Q[1]+1):
            b1=brute(P,Q,sigma)
            if b1 is None: continue
            tot+=1
            p1=predict(P,Q,sigma)
            if b1!=p1: bad+=1
            # criterion check
            A=max(P[0],sigma-Q[1]); Z=min(Q[0],sigma-P[1])
            u,v=min(0,sigma),max(0,sigma)
            F=max(0,min(Z,v)-max(A,u)+1)
            L=max(0,min(Z,u-1)-A+1); R=max(0,Z-max(A,v+1)+1)
            pred = True if F==0 else crit(F,min(L,R),abs(L-R))
            if pred != is_concave_pospart(b1): badc+=1
print(f"  {tot} (box,sigma) cases:  shape formula wrong {bad};  criterion wrong {badc}")
