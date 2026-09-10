from rim import *

# FACT 1: |C(mu)| = n  (border strips of mu <-> cells of mu)
bad=[]
for n in range(1,10):
    for mu in partitions(n):
        if len(predecessors(mu)) != n: bad.append((mu,len(predecessors(mu))))
print("FACT1  |C(mu)|=n :", "OK for all mu, n<=9" if not bad else bad, flush=True)

# FACT 2: at t=-1 the solution is wtB(mu,gamma) = (-1)^ht(mu/gamma) / n
bad=[]
for n in range(1,9):
    for mu in partitions(n):
        rows, cols, M, b = local_system(mu)
        P = predecessors(mu)
        w = sp.Matrix([sp.Rational((-1)**P[g], n) for g in cols])
        r = sp.expand(M.subs(t,-1)*w - b.subs(t,-1))
        if any(x != 0 for x in r): bad.append(mu)
print("FACT2  wtB=(-1)^ht/n solves at t=-1 :", "OK for all mu, n<=8" if not bad else bad, flush=True)

# FACT 3: generic consistency, n=8,9
def gen_consistent(mu):
    rows, cols, M, b = local_system(mu)
    A = M.row_join(b)
    R, piv = A.rref(simplify=True)
    return (A.cols-1) not in piv, len(piv)
for n in (8,9):
    sol=[]
    for mu in partitions(n):
        ok, r = gen_consistent(mu)
        if ok: sol.append(mu)
    print(f"FACT3  n={n}: generically consistent mu = {sol}   ({len(sol)}/{len(list(partitions(n)))})", flush=True)
