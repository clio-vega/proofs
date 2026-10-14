"""Check the generalised census against the dict-Laurent ground truth, then the
whole general-k c_{(1^n)} against engine M1."""
import sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1014-q458')
from itertools import product
from collections import defaultdict
from censusk import L_poly, poly_pow

def edots(n):
    """Full edot_r as Laurent dicts (ground truth)."""
    vs = []
    for i in range(n):
        v=[0]*n; v[i]=1; vs.append(tuple(v))
        w=[0]*n; w[i]=-1; vs.append(tuple(w))
    e=[defaultdict(int) for _ in range(2*n+1)]; e[0][tuple([0]*n)]=1
    for v in vs:
        for r in range(2*n,0,-1):
            for mon,c in list(e[r-1].items()):
                e[r][tuple(a+b for a,b in zip(mon,v))] += c
    return [dict(x) for x in e]

def direct(n,k,b,rs):
    """[x^{b(1..1)}] prod edot_{r_i} by brute Laurent multiplication."""
    E = edots(n)
    for r in rs:
        if r<0 or r>2*n: return 0
    cur = {tuple([0]*n):1}
    for r in rs:
        nxt=defaultdict(int)
        for m1,c1 in cur.items():
            for m2,c2 in E[r].items():
                nxt[tuple(a+bb for a,bb in zip(m1,m2))]+=c1*c2
        cur=dict(nxt)
    return cur.get(tuple([b]*n),0)

ok=fail=0
for n in range(1,6):
  for k in (1,2,3):
    b=1-k
    Ln = poly_pow(L_poly(k,b), n, k, 2*n)
    for rs in product(range(0,2*n+1), repeat=k):
        got = Ln.get(tuple(rs),0); want = direct(n,k,b,list(rs))
        if got==want: ok+=1
        else:
            fail+=1
            if fail<10: print("MISMATCH",n,k,rs,got,want)
print(f"GENERALISED CENSUS: ok={ok} fail={fail}")
