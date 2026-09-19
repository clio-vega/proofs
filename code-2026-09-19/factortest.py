import sympy as sp, itertools
from sixvertex import Z, schur_weights, ff_weights

x = sp.symbols('x1:6'); eps = sp.Symbol('eps')

print("CONJECTURE:  Z_eps = Z_0 * prod_i (1 + eps*x_i)^k   for some k.")
print("%-4s %-4s %-14s %-14s %-6s %s" % ("n","m","top","bot","k","Z_eps == Z_0*prod(1+eps x_i)^k ?"))
for n in (2, 3):
    for m in range(n+1, n+5):
        w0 = [schur_weights(x[i]) for i in range(n)]
        wd = [ff_weights(x[i], eps) for i in range(n)]
        top = [1 if j < n else 0 for j in range(m)]
        for botset in itertools.combinations(range(m), n):
            bot = [1 if j in botset else 0 for j in range(m)]
            z0 = Z(n, m, top, bot, [0]*n, [0]*n, w0)
            if z0 == 0:
                continue
            ze = Z(n, m, top, bot, [0]*n, [0]*n, wd)
            # find k by comparing total eps-degree
            q = sp.simplify(sp.cancel(sp.expand(ze) / sp.expand(z0)))
            found = None
            for k in range(0, 2*m+1):
                cand = sp.prod([(1 + eps*x[i])**k for i in range(n)])
                if sp.simplify(sp.expand(ze - z0*cand)) == 0:
                    found = k; break
            print("%-4d %-4d %-14s %-14s %-6s %s" %
                  (n, m, str(tuple(j for j in range(m) if top[j])), str(botset),
                   str(found), "YES" if found is not None else "NO  (quotient=%s)" % q))
