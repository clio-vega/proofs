from rim import *
# Does deleting row mu leave full column rank n?  (=> w=0 forced => row mu gives 0=1)
print(" n | mu                | rank(M minus row mu) | =n? | generic consistent?")
for n in range(3, 9):
    for mu in partitions(n):
        rows, cols, M, b = local_system(mu)
        i = rows.index(tuple(mu))
        Mm = M[[k for k in range(M.rows) if k != i], :]
        r = Mm.rank()
        A = M.row_join(b); R,piv = A.rref(simplify=True)
        gc = (A.cols-1) not in piv
        print(f"{n:2d} | {str(mu):18s} | {r:3d} | {'YES' if r==n else 'no ':4s}| {gc}")
