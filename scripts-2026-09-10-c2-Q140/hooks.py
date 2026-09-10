from rim import *
print(" n | mu                  | det(hook rows) [mu-row deleted if mu is a hook]")
for n in range(4, 9):
    hooksn = [tuple([a]+[1]*(n-a)) for a in range(1, n+1)]
    for mu in partitions(n):
        rows, cols, M, b = local_system(mu)
        idx = [rows.index(h) for h in hooksn if h != tuple(mu)]
        S = M[idx, :]
        r = S.rank()
        extra = ""
        if S.rows == S.cols:
            extra = " det=" + str(sp.factor(S.det()))
        print(f"{n:2d} | {str(mu):19s} | rank {r}/{n} {'FULL' if r==n else '    '}{extra}")
