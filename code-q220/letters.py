"""Letter-level comparison: E(S,T) vs B(S,T) = { b_x : x in X }."""
import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *

def B_letters(S,T,n):
    out={}
    for x in X_set(S,T,n):
        r=ms_etilde(S,T,x,n)
        if r is not None: out.setdefault(r[0],[]).append(x)
    return out

def scan(n):
    rows=[]
    for S,T,wv in additive_pairs(n):
        X=X_set(S,T,n)
        if not X: continue
        E=clio_letters(S,T,n); B=B_letters(S,T,n)
        rows.append((S,T,X,E,set(B)))
    return rows

tot=collections.Counter()
for n in range(3,8):
    rows=scan(n)
    eq=sum(1 for r in rows if r[3]==r[4])
    sub=sum(1 for r in rows if r[4]<=r[3])
    print(f"n={n} (X nonempty): {len(rows)} pairs, B subset E on {sub}, E==B on {eq}")
    tot['rows']+=len(rows); tot['sub']+=sub; tot['eq']+=eq
print("TOTAL", dict(tot))

# what do excess letters look like?  n=6, tabulate
print("\n--- excess examples, n=5, |S|>|T| ---")
for S,T,X,E,B in scan(5):
    if E-B and len(S)>len(T):
        print(f"S={sorted(S)} T={sorted(T)} X={X} E={sorted(E)} B={sorted(B)} excess={sorted(E-B)}")
