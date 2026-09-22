"""Final table for the write-up.  Counts only -- no fractions (Q217 unresolved)."""
import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *
from emptyX import case3_bs, ext_move

def MS_letters(S,T,n):
    X=X_set(S,T,n)
    out=set()
    if X:
        for x in X:
            r=ms_etilde(S,T,x,n)
            if r is not None: out.add(r[0])
    else:
        for b,t in case3_bs(S,T,n):
            r=ext_move(S,T,b,n,'A')
            if r: out.add(r[0])
    return out

print(f"{'n':>2} {'pairs':>6} {'X=0':>5} {'|E|':>6} {'|B|':>6} {'B<=E':>6} {'B=E':>6} {'excess moves':>13} {'pairs w/ excess':>15}")
first_excess=None
for n in range(3,8):
    N=Xe=SE=EQ=0; tE=tB=exc=pexc=0
    for S,T,_ in additive_pairs(n):
        N+=1
        if not X_set(S,T,n): Xe+=1
        E=clio_letters(S,T,n); B=MS_letters(S,T,n)
        tE+=len(E); tB+=len(B)
        if B<=E: SE+=1
        if B==E: EQ+=1
        else:
            exc+=len(E-B); pexc+=1
            if first_excess is None and E-B:
                first_excess=(n,sorted(S),sorted(T),sorted(X_set(S,T,n)),sorted(E),sorted(B))
    print(f"{n:>2} {N:>6} {Xe:>5} {tE:>6} {tB:>6} {SE:>6} {EQ:>6} {exc:>13} {pexc:>15}")
print("\nsmallest pair with excess:", first_excess)

# smallest n with excess, enumerated fully
n=5
print(f"\nall excess instances at n={n}:")
for S,T,_ in additive_pairs(n):
    E=clio_letters(S,T,n); B=MS_letters(S,T,n)
    if E-B:
        print(f"  S={sorted(S)} T={sorted(T)} X={sorted(X_set(S,T,n))} E={sorted(E)} B={sorted(B)} excess={sorted(E-B)}")
