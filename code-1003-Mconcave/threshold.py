"""Exact threshold in G = n-m.

Claim (S6): for a NONEMPTY effective slice, d <= 3G.   [since sum_i u_i >= 0 in (S3)
gives sum_i (delta_i)_- <= 2G, hence sum delta_i >= -2G, hence d = G - sum delta_i <= 3G.]
So sweeping d <= 3G is COMPLETE for a given (m, n): no slice is missed.

Then: for each m, find the least G with a counterexample.
"""
import gen, tent, sys
from collections import Counter
from bigsweep import ybox, profiles, is_cc

def run(m, Gmax):
    out={}
    for G in range(0, Gmax+1):
        n=m+G
        fails=[]; tot=0; dmaxseen=0
        for (mu,lam) in gen.pairs(n,m,3*G if G else 0):
            d=sum(lam)-sum(mu); P,Q=ybox(mu,lam,n,m)
            if any(Q[i]<P[i] for i in range(m)): continue
            dmaxseen=max(dmaxseen,d)
            F=profiles(P,Q)
            for b in range(0,d+1):
                sigma=n-m-d+b
                if sigma not in F: continue
                qd=F[sigma]; kmin,kmax=min(qd),max(qd)
                N=[qd.get(k,0) for k in range(kmin,kmax+1)]
                tot+=1
                if not is_cc(N): fails.append((mu,lam,b,P,Q,sigma,N))
        out[G]=(tot,len(fails),fails[:2],dmaxseen)
        print(f"  m={m} G={G} (n={n}, d<={3*G}): slices={tot} FAIL={len(fails)} maxd_nonempty={dmaxseen}")
        if fails: print(f"     first: mu={fails[0][0]} lam={fails[0][1]} b={fails[0][2]} P={fails[0][3]} Q={fails[0][4]} sigma={fails[0][5]} N={fails[0][6]}")
        sys.stdout.flush()
    return out

# calibration of (S6): d<=3G on nonempty slices
st=Counter()
for m in range(2,5):
    for n in range(m,m+6):
        G=n-m
        for (mu,lam) in gen.pairs(n,m,3*G+5 if G else 5):
            d=sum(lam)-sum(mu); P,Q=ybox(mu,lam,n,m)
            if any(Q[i]<P[i] for i in range(m)): st['emptybox']+=1; continue
            st['d_le_3G_ok' if d<=3*G else 'd_le_3G_BAD']+=1
print("(S6) calibration:",dict(st)); sys.stdout.flush()

for m in [int(x) for x in sys.argv[1].split(',')]:
    run(m, int(sys.argv[2]))
