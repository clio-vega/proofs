"""The two sharp box bounds, and the resulting FINITE check.

With g_i = lam_i-lam_{i-1}-1 >= 0 (sum = G = n-m), delta_i = mu_i-lam_{i-1}-1,
    P_i = -g_{i-1} + (delta_i+g_{i-1})_+ ,   Q_i = g_i - (delta_{i+1})_- ,
so with X = sum (delta_i)_-,  Y = sum (delta_i+g_{i-1})_+ ,
    (B1)  sum_i u_i = 2G - X - Y  and  Y >= G - X + sum (delta_i)_+ ,
          hence   sum_i u_i + sum_i (delta_i)_+ <= G,  in particular  sum_i u_i <= G.
    (B2)  sum_i Q_i + sum_i (P_i)_-  =  sum_i u_i + sum_i (delta_i)_+  <= G.
Consequences: |{i:u_i>0}| <= G,  sum (Q_i)_+ <= G,  sum (P_i)_- <= G,  P_i in [-G,G].
Frozen coordinates (u_i=0) only translate sigma and k, so they do not affect the
SHAPE of N.  Hence concavity for a given G is a FINITE question."""
import gen, tent, sys
from collections import Counter
from itertools import product
from bigsweep import ybox, is_cc

import os
SKIP=os.environ.get("SKIPCAL")
# ---------- verify (B1),(B2) on real shapes, + plant a violation ----------
st=Counter(); ex={}
for m in ([] if SKIP else range(2,6)):
    for n in range(m,m+7):
        G=n-m
        for (mu,lam) in gen.pairs(n,m,min(3*G+3,14)):
            P,Q=ybox(mu,lam,n,m)
            if any(Q[i]<P[i] for i in range(m)): continue
            g=[tent.lam_at(lam,i,n,m)-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
            dl=[mu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
            U=sum(Q[i]-P[i] for i in range(m)); DP=sum(max(x,0) for x in dl)
            st['B1_ok' if U+DP<=G else 'B1_BAD']+=1
            if U+DP>G: ex.setdefault('B1',(m,n,mu,lam,P,Q,U,DP,G))
            lhs=sum(Q)+sum(max(-x,0) for x in P)
            st['B2id_ok' if lhs==U+DP else 'B2id_BAD']+=1
            st['B2_ok' if lhs<=G else 'B2_BAD']+=1
            st['nfroz_ok' if sum(1 for i in range(m) if Q[i]>P[i])<=G else 'nfroz_BAD']+=1
            st['Prange_ok' if all(-G<=P[i]<=G for i in range(m)) else 'Prange_BAD']+=1
print("real-shape check:",dict(st)); [print(' EX',k,v) for k,v in ex.items()]
# refusal probe: the bounds must REFUSE a perturbed box
badref=0; tot=0
for m in ([] if SKIP else range(2,5)):
    for n in range(m,m+5):
        G=n-m
        for (mu,lam) in gen.pairs(n,m,2*G+2):
            P,Q=ybox(mu,lam,n,m)
            if any(Q[i]<P[i] for i in range(m)): continue
            Qp=[Q[0]+1]+Q[1:]            # widen one coordinate by 1
            tot+=1
            if sum(Qp)+sum(max(-x,0) for x in P)>G: badref+=1
print(f"  refusal probe (Q_1 -> Q_1+1 must break B2): refused {badref}/{tot}")
sys.stdout.flush()

# ---------- the finite check ----------
def Nprof(P,Q,sigma):
    cur={0:{0:1}}
    for i in range(len(P)):
        new={}
        for xd,qd in cur.items():
            for t in range(P[i],Q[i]+1):
                dq=max(-t,0); tgt=new.setdefault(xd+t,{})
                for q,c in qd.items(): tgt[q+dq]=tgt.get(q+dq,0)+c
        cur=new
    if sigma not in cur: return None
    qd=cur[sigma]; return [qd.get(k,0) for k in range(min(qd),max(qd)+1)]

def finite_check(G):
    """enumerate every non-frozen box allowed by (B1)+(B2) with budget G."""
    fails=[]; tot=0
    def rec(boxes, usum, budget):
        nonlocal tot
        if boxes:
            P=[b[0] for b in boxes]; Q=[b[1] for b in boxes]
            for sigma in range(sum(P), sum(Q)+1):
                N=Nprof(P,Q,sigma)
                if N is None: continue
                tot+=1
                if not is_cc(N): fails.append((tuple(P),tuple(Q),sigma,N))
        if usum>=G or len(boxes)>=G: return
        for u in range(1, G-usum+1):
            for p in range(-G, G+1):
                q=p+u
                if q>G: continue
                # (B2) running budget:  sum Q_i + sum (P_i)_-  <= G
                cost=q+max(-p,0)
                if cost>budget: continue
                rec(boxes+[(p,q)], usum+u, budget-cost)
    rec([], 0, G)
    return tot, fails

import os
for G in [int(x) for x in os.environ.get("GS","0,1,2,3,4,5,6").split(",")]:
    tot,fails=finite_check(G)
    print(f"  G={G}: {tot} reduced (box,sigma) cases in the allowed class, FAIL={len(fails)}")
    for f in fails[:3]: print("     ",f)
    sys.stdout.flush()
