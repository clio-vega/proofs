"""REFUSAL CONTROL (brief control 1): are (B1) and (B2) load-bearing, and tight?

The G<=5 theorem spends exactly two bounds:
   (B1) sum_i u_i <= G            [u_i = Q_i - P_i]
   (B2) sum_i Q_i + sum_i (P_i)_- <= G
Decouple their budgets.  If relaxing either one alone already produces a failure at
budget 5, the bound was not load-bearing / not tight.  Prediction, from the witness
(sum u = 6 and sum Q + sum(P)_- = 6): the LEAST budgets admitting a failure are
U=6 AND B=6 -- neither alone suffices.
"""
import sys
from bigsweep import is_cc

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

def check(U,B,R=None):
    """enumerate non-frozen boxes with sum u_i <= U and sum Q_i + sum (P_i)_- <= B;
    coordinate range |P_i|,|Q_i| <= R (default max(U,B))."""
    if R is None: R=max(U,B)
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
        if usum>=U: return
        for u in range(1, U-usum+1):
            for p in range(-R, R+1):
                q=p+u
                if q>R: continue
                cost=q+max(-p,0)
                if cost>budget: continue
                rec(boxes+[(p,q)], usum+u, budget-cost)
    rec([],0,B)
    return tot,fails

print("panel: least (U,B) admitting a failure   [prediction: U=6 and B=6, neither alone]")
for U,B in [(5,5),(6,5),(5,6),(6,6),(5,8),(8,5),(5,12),(12,5)]:
    tot,f=check(U,B)
    tag = "FAIL" if f else "clean"
    print(f"   U={U:>2} B={B:>2}: {tot:>8} cases -> {tag}" + (f"   first {f[0]}" if f else ""))
    sys.stdout.flush()
