"""First-order obstruction at the anchor t=-1.

M(t)_{lambda,gamma} = t^{ht(lambda/gamma)}.  At t=-1 the system is solvable with the
UNIQUE solution beta^0_gamma = (-1)^{ht(mu/gamma)}/n  (Euler identity sum_L p_L p_L^perp = n.id).
If rank M(-1) = n, the candidate solution w(t) is regular at t=-1 and

      r'(-1)  in   M_1 beta^0  +  im M_0,      M_1 = dM/dt|_{t=-1}.

So  Omega(mu) := M_1 beta^0  not in  W(mu) := im M_0  ==>  INCONSISTENT over Q(t).
Everything below is over Q -- no t anywhere.
"""
from rim import *

def M0_M1(mu):
    n = sum(mu)
    rows = list(partitions(n))
    P = predecessors(mu)
    cols = sorted(P.keys(), key=lambda g: (sum(g), g))
    M0 = sp.zeros(len(rows), len(cols)); M1 = sp.zeros(len(rows), len(cols))
    for i, lam in enumerate(rows):
        for j, gam in enumerate(cols):
            L = n - sum(gam)
            for g, h in remove_rim_abacus(lam, L):
                if g == gam:
                    M0[i, j] = (-1)**h
                    M1[i, j] = h * (-1)**(h-1)
    beta0 = sp.Matrix([(-1)**P[g] for g in cols])   # x n dropped: scalar
    return rows, cols, M0, M1, beta0

def in_colspan(A, v):
    return sp.Matrix(A).rank() == sp.Matrix(A.row_join(v)).rank()

print(" n | mu                | rankM0 (=n?) | Omega in W(mu)? | generic consistent?")
from scan import consistent  # reuse
for n in range(2, 9):
    for mu in partitions(n):
        rows, cols, M0, M1, beta0 = M0_M1(mu)
        r0 = M0.rank()
        Om = M1 * beta0
        inW = in_colspan(M0, Om)
        _,_,M,b = local_system(mu)
        gc = consistent(M, b)
        flag = "" if (not inW) == (not gc) or True else ""
        print(f"{n:2d} | {str(mu):18s} | {r0:2d} {'OK' if r0==n else '**DROP**':8s} | "
              f"{str(inW):5s} | {gc}")
