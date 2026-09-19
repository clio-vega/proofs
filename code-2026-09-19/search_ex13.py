"""Find the reading of Le-Nguyen 2608.13544 Example 1.3 that reproduces 6 and 5.

Stated: nu=(9,5,3), pi=(6,5,4,2), x=(6,4,3), y=(5,3,1),
        f(lam) = c^nu_{lam,(pi-lam)},
        f(x)f(y) = 6 > 5 = f(ceil((x+y)/2)) f(floor((x+y)/2)).
The lengths of pi (4) and nu,x,y (3) do not match, so we search.
"""
from lr import lr_coeff
from itertools import product

NU3 = (9, 5, 3)
x3, y3 = (6, 4, 3), (5, 3, 1)
m1 = tuple(-((-(a + b)) // 2) for a, b in zip(x3, y3))   # ceil
m2 = tuple((a + b) // 2 for a, b in zip(x3, y3))          # floor
print("x =", x3, " y =", y3, " ceil =", m1, " floor =", m2)
print("sums:", sum(x3), sum(y3), sum(m1), sum(m2), "| x+y =", sum(x3)+sum(y3),
      " m1+m2 =", sum(m1)+sum(m2))

def f(nu, pi, lam):
    mu = tuple(p - l for p, l in zip(pi, lam))
    if any(v < 0 for v in mu):
        return 0
    if any(mu[i] < mu[i+1] for i in range(len(mu)-1)):
        return 0
    return lr_coeff(nu, lam, mu)

print("\n=== Reading A: n=3, pi in Z^3 with |pi|=17, brute force ===")
hits = []
for p1 in range(18):
    for p2 in range(p1 + 1):
        p3 = 17 - p1 - p2
        if p3 < 0 or p3 > p2:
            continue
        pi = (p1, p2, p3)
        fx, fy, fm1, fm2 = (f(NU3, pi, l) for l in (x3, y3, m1, m2))
        if fx*fy == 6 and fm1*fm2 == 5:
            hits.append((pi, fx, fy, fm1, fm2, "EXACT 6>5"))
        elif fx*fy > fm1*fm2 and fx*fy > 0:
            hits.append((pi, fx, fy, fm1, fm2, "violation"))
for h in hits:
    print("  pi=%-12s f(x)=%d f(y)=%d f(m1)=%d f(m2)=%d  %s  [%d vs %d]"
          % (h[0], h[1], h[2], h[3], h[4], h[5], h[1]*h[2], h[3]*h[4]))
if not hits:
    print("  (none)")

print("\n=== Reading B: labels swapped, nu=(6,5,4,2), pi=(9,5,3) ===")
NU4, PI3 = (6,5,4,2), (9,5,3)
vals = {}
for name, lam in [("x",x3),("y",y3),("ceil",m1),("floor",m2)]:
    mu = tuple(p-l for p,l in zip(PI3,lam))
    c = lr_coeff(NU4, lam, mu)
    vals[name] = c
    print("  f(%-5s=%s) : mu = %-10s c^%s_{%s,%s} = %d" % (name,lam,mu,NU4,lam,mu,c))
lhs = vals["x"]*vals["y"]; rhs = vals["ceil"]*vals["floor"]
print("  f(x)f(y) = %d   f(ceil)f(floor) = %d   -> %s"
      % (lhs, rhs, "MATCHES 6 > 5" if (lhs,rhs)==(6,5) else "does NOT match"))
