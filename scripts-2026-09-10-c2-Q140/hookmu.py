from rim import *
print("mu a HOOK: does {hooks != mu} + row (2,2,1^{n-4}) have rank n?")
for n in range(4, 10):
    hooksn = [tuple([a]+[1]*(n-a)) for a in range(1, n+1)]
    extra = tuple([2,2]+[1]*(n-4))
    for a in range(1, n+1):
        mu = tuple([a]+[1]*(n-a))
        rows, cols, M, b = local_system(mu)
        idx = [rows.index(h) for h in hooksn if h != mu] + [rows.index(extra)]
        S = M[idx, :]
        print(f"  n={n} mu={str(mu):20s} rank {S.rank()}/{n} {'FULL' if S.rank()==n else '<<< SHORT'}")

print("\nGeneral mu: rank of the family F1 = {lambda : lambda_3 <= 1}, row mu removed")
for n in range(4, 9):
    for mu in partitions(n):
        rows, cols, M, b = local_system(mu)
        idx = [i for i,l in enumerate(rows) if (len(l)<3 or l[2]<=1) and l != tuple(mu)]
        S = M[idx, :]
        print(f"  n={n} mu={str(mu):20s} |F1|={len(idx):2d} rank {S.rank()}/{n} {'FULL' if S.rank()==n else '<<<'}")
