from math import comb as C
from functools import lru_cache

# Murnaghan-Nakayama: character chi^lambda(rho)
# chi^lambda(rho): remove border strips of sizes rho (in any fixed order), sign (-1)^(height-1)
def border_strips(lam, k):
    # yield (resulting partition, height-1) for removing a border strip of size k from lam
    lam = list(lam)
    n = len(lam)
    # represent by beta-set / use the standard "remove rim hook" by scanning cells
    # Use first-column hook lengths? Simpler: iterate over starting cell.
    # We'll use the boundary-path method via partition profile.
    # Use the classic algorithm with "diagonals".
    res = []
    # A rim hook of length k: choose cells. Use the method: for each cell (i,j) that is a
    # "top-right" start... implement via removing from the outer boundary.
    # Standard approach: a rim hook corresponds to choosing rows i<=i' such that it's connected.
    # Use beta-numbers (first-column hook lengths): lam_i + (n - i).
    beta = sorted([lam[i] + (n-1-i) for i in range(n)], reverse=True)
    bset = set(beta)
    for b in list(bset):
        if b-k >= 0 and (b-k) not in bset:
            newbeta = (bset - {b}) | {b-k}
            nb = sorted(newbeta, reverse=True)
            # convert back to partition
            mu = [nb[i] - (len(nb)-1-i) for i in range(len(nb))]
            mu = [x for x in mu if x>0]
            # height = number of rows spanned = b - (b-k) ... actually height-1 = (number of
            # beta values strictly between b-k and b that are present) 
            ht = sum(1 for x in bset if (b-k) < x < b)
            res.append((tuple(mu), ht))
    return res

@lru_cache(maxsize=None)
def chi(lam, rho):
    lam = tuple(x for x in lam if x>0)
    if sum(lam)==0 and len(rho)==0:
        return 1
    if len(rho)==0:
        return 1 if sum(lam)==0 else 0
    k = rho[0]
    rest = rho[1:]
    total = 0
    for mu, ht in border_strips(lam, k):
        total += ((-1)**ht) * chi(mu, rest)
    return total

def Mj_char(m,a,b,j):
    # M_j = 2^{-j} sum_t C(j,t)(-1)^t chi^lambda(2^t 1^{2(m-t)})
    lam = (a,b) if b>0 else (a,)
    tot = 0
    for t in range(0,j+1):
        rho = tuple([2]*t + [1]*(2*(m-t)))
        tot += C(j,t)*((-1)**t)*chi(lam, rho)
    assert tot % (2**j) == 0, (m,a,b,j,tot)
    return tot // (2**j)

ok=True
for m in range(1,9):
    for b in range(0,m+1):
        a=2*m-b
        for j in range(0,b+1):
            lhs = Mj_char(m,a,b,j)
            rhs = C(2*(m-j), b-j) - (C(2*(m-j), b-j-1) if b-j-1>=0 else 0)
            if lhs!=rhs:
                ok=False
                print("MISMATCH",m,a,b,j,lhs,rhs)
print("M_j (def via characters) == C(2(m-j),b-j)-C(2(m-j),b-j-1):", ok)
