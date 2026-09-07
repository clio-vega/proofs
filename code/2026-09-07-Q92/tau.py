"""Check: for a connected ribbon mu/lam of size g arising from bead move b->b+g,
   m(b+j) == tau_j := [the ribbon turns (changes row) between its j-th and
   (j+1)-st cell, read from the SW tail to the NE head],  1<=j<=g-1."""
import sys, sympy as sp
sys.path.insert(0,'/home/clio/projects/scratch/q92')
import engineA as A, engineB as B

def ribbon_cells_in_order(lam, mu):
    sk = sorted(A.cells(mu)-A.cells(lam))          # (row,col), row 0 = top
    # tail = SW end = largest row, smallest col; walk along the strip
    start = max(sk, key=lambda p:(p[0], -p[1]))
    order=[start]; seen={start}
    while len(order)<len(sk):
        i,j = order[-1]
        nxt = None
        for cand in ((i,j+1),(i-1,j)):             # head is NE: go right or up
            if cand in sk and cand not in seen: nxt=cand; break
        if nxt is None: return None
        order.append(nxt); seen.add(nxt)
    return order

def bead_move(lam, mu, g):
    Ml, Mm = B.maya(lam), B.maya(mu)
    rem = sorted(Ml-Mm); add = sorted(Mm-Ml)
    if len(rem)!=1 or len(add)!=1 or add[0]-rem[0]!=g: return None
    return rem[0]

ok=bad=0
for n in range(0,8):
    for lam in A.partitions(n):
        for g in range(1,7):
            for mu in A.R_e(lam,g):
                b = bead_move(lam,mu,g)
                order = ribbon_cells_in_order(lam,mu)
                if b is None or order is None: bad+=1; print("STRUCT FAIL",lam,mu,g); continue
                M = B.maya(lam)
                tau = [1 if order[j-1][0]!=order[j][0] else 0 for j in range(1,g)]
                mm  = [B.occ(M,b+j) for j in range(1,g)]
                if tau==mm: ok+=1
                else:
                    bad+=1
                    if bad<4: print("FAIL",lam,mu,g,"tau",tau,"m",mm)
print(f"tau_j == m(b+j):  {ok} agree, {bad} fail   (|lam|<=7, ribbon size <=6)")
