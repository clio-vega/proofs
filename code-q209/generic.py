"""Does (A)+(B') imply M-convexity?  Exhaustive search for counterexamples."""
from itertools import combinations_with_replacement, product

def simplex(r, d):
    pts=[]
    def rec(pre, rem, k):
        if k==1: pts.append(tuple(pre+[rem])); return
        for x in range(rem+1): rec(pre+[x], rem-x, k-1)
    rec([], d, r); return pts

def satisfies_B(W):
    S=set(W); r=len(next(iter(S)))
    for a in S:
        for i in range(r):
            for j in range(r):
                if i>=j: continue
                tot=a[i]+a[j]
                vals=sorted(b[i] for b in S
                            if all(b[s]==a[s] for s in range(r) if s not in (i,j))
                            and b[i]+b[j]==tot)
                if vals != list(range(vals[0], vals[-1]+1)): return False   # not an interval
                if vals[0]+vals[-1] != tot: return False                    # not palindromic
    return True

def is_M(W):
    S=set(W); r=len(next(iter(S)))
    for a in S:
        for b in S:
            for i in range(r):
                if a[i]<=b[i]: continue
                if not any(a[j]<b[j] and tuple(a[k]-(k==i)+(k==j) for k in range(r)) in S
                           for j in range(r)): return False
    return True

for r,d in ((3,2),(3,3),(3,4),(4,2),(4,3),(2,5)):
    pts=simplex(r,d); N=len(pts); bad=[]; nB=0
    for mask in range(1, 1<<N):
        W=[pts[i] for i in range(N) if mask>>i & 1]
        if not satisfies_B(W): continue
        nB += 1
        if not is_M(W): bad.append(W)
    print(f"r={r} d={d}: |simplex|={N}  sets satisfying (A)+(B')={nB}  "
          f"of which NOT M-convex = {len(bad)}")
    for W in bad[:3]: print("     counterexample:", sorted(W))
