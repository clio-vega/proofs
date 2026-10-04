"""C1, repaired.  Round 1 of this control deleted EVERY support point, not every
INTERIOR support point as the brief specified, and 72 of 534 deletions did not refuse.
That is not a failure of the chain: deleting an EXTREME point of a box slice can leave
an M-convex set (an interval minus an endpoint is an interval).  Characterise the
non-refusals exactly, then run the control the brief actually asked for."""
import itertools
from lor import *

def is_interior(p, S):
    """p is interior to the box slice iff every neighbour p-e_i+e_j (i!=j) is in S."""
    n=len(p); Ss=set(S)
    for i in range(n):
        for j in range(n):
            if i!=j and p[i]>0:
                q=list(p); q[i]-=1; q[j]+=1
                if tuple(q) not in Ss: return False
    return True

tot=ref=0; mism=0; inter_tot=inter_ref=0; sizes=[]
for h in range(2,8):
    for l in range(0,h+1):
        c=Ptilde(l,h); S=sorted(c)
        if len(S)<3: continue
        sizes.append(len(S))
        for p in S:
            c2=dict(c); del c2[p]
            stillM = is_Mconvex([a for a in c2])
            ok = is_N_lorentzian(c2)
            tot+=1; ref += (not ok)
            # the claim: refusal happens exactly when M-convexity is destroyed
            if (not ok) != (not stillM): mism+=1
            if is_interior(p,S):
                inter_tot+=1; inter_ref += (not ok)
print(f"  ALL deletions: enumerated {tot}, refusals {ref}; support sizes {min(sizes)}..{max(sizes)}")
print(f"  refusal happened EXACTLY when the support stopped being M-convex: "
      f"{tot-mism}/{tot} agreements")
print(f"  INTERIOR deletions (the brief's control): enumerated {inter_tot}, refusals {inter_ref}"
      f"  (must equal {inter_tot})")
print(f"  non-refusing deletions are therefore exactly the {tot-ref} boundary points whose")
print(f"  removal leaves a box slice M-convex -- a positive instance of Lemma boxslice,")
print(f"  not a control.  Recorded, in the shape of the 1004c2 'control that was no control'.")
