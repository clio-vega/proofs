"""Test the flagged claim: on Sigma^circ_b, is  w_i(nu) = g_i - (y_{i+1})_- - (y_i)_+ + 1 ?
y_i = nu_i - lam_{i-1} - 1,  g_i = lam_i - lam_{i-1} - 1."""
import gen, tent
from collections import Counter
st=Counter(); ex={}
for m in range(2,6):
    for n in range(m,m+6):
        for (mu,lam) in gen.pairs(n,m,8):
            g=[tent.lam_at(lam,i,n,m)-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
            for nu in gen.box(mu,n,m):
                off,co=gen.f_nu(nu,lam,n,m)
                if not co: continue
                lr=gen.LR(nu,lam,n,m); w=[R-L+1 for (L,R) in lr]
                y=[nu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
                pred=[g[i]-max(-y[(i+1)%m],0)-max(y[i],0)+1 for i in range(m)]
                if w==pred: st['ok']+=1
                else:
                    st['BAD']+=1
                    if len(ex)<3: ex[len(ex)]=(m,n,mu,lam,nu,y,g,w,pred)
                # and the sum shadow: sum w_i = 2*Lambda + m
                if sum(w)==2*tent.Lambda_formula(nu,lam,n,m)+m: st['sum_ok']+=1
                else: st['sum_BAD']+=1
print(dict(st))
for v in ex.values(): print(' EX',v)
