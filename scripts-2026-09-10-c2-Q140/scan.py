from rim import *

def consistent(M, b, val=None):
    """rank test over Q(t) (val=None) or at t=val"""
    if val is not None:
        M = M.subs(t, val); b = b.subs(t, val)
    Mb = M.row_join(b)
    return sp.Matrix(M).rank() == sp.Matrix(Mb).rank()

# calibrate: the paper's certificate
rows, cols, M, b = local_system((4,))
c = sp.Matrix([[t**2*(t+1), -t**2*(t+1), -t*(3*t-1), -(t-1)**2, 2*(t-1)]])
print("c*M =", sp.simplify(c*M), "   c*b =", sp.simplify((c*b)[0,0]))

print("\n n | mu           | rows x cols | rank | t=-1 | t=0 | t=1 | generic")
tally = {}
for n in range(2, 8):
    for mu in partitions(n):
        rows, cols, M, b = local_system(mu)
        r = sp.Matrix(M).rank()
        g = consistent(M, b)
        res = {v: consistent(M, b, v) for v in (-1, 0, 1)}
        tally.setdefault(n, []).append((mu, g, res))
        print(f"{n:2d} | {str(mu):13s} | {M.rows}x{M.cols} | {r:2d} | "
              f"{str(res[-1]):5s}|{str(res[0]):5s}|{str(res[1]):5s}| {g}")
print()
for n in sorted(tally):
    L = tally[n]
    print(f"n={n}: t=-1 {sum(1 for _,g,r in L if r[-1])}/{len(L)}  "
          f"t=0 {sum(1 for _,g,r in L if r[0])}/{len(L)}  "
          f"t=1 {sum(1 for _,g,r in L if r[1])}/{len(L)}  "
          f"generic {sum(1 for _,g,r in L if g)}/{len(L)}"
          f"   [generic-solvable mu: {[mu for mu,g,r in L if g]}]")
