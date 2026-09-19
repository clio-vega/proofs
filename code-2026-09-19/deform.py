import sympy as sp, itertools
from sixvertex import Z, schur_weights, ff_weights, check_ff

x = sp.symbols('x1:6'); eps = sp.Symbol('eps')

print("=== n=3 control: five-vertex FF point gives Schur ===")
n, m = 3, 7
ws = [schur_weights(x[i]) for i in range(n)]
top = [1 if j in (0,1,2) else 0 for j in range(m)]
for botset in itertools.combinations(range(m), n):
    if min(botset) < 3: continue
    bot = [1 if j in botset else 0 for j in range(m)]
    z = Z(n, m, top, bot, [0]*n, [0]*n, ws)
    if z != 0:
        # lambda from bot positions: bot sorted ascending = lambda_{n-i} + n + i
        lam = tuple(sorted([b - n - i for i, b in enumerate(sorted(botset))], reverse=True))
        s = sp.expand(sp.Matrix(n,n, lambda r,c: x[c]**(lam[r]+n-1-r)).det() /
                      sp.Matrix(n,n, lambda r,c: x[c]**(n-1-r)).det())
        s = sp.simplify(sp.expand(s))
        match = sp.simplify(sp.expand(z - s)) == 0
        print("  bot=%s  lambda=%s  Z=%s   == s_lambda ? %s" % (botset, lam, sp.factor(z), match))

print("\n=== deformed free-fermion Z (eps != 0) ===")
n, m = 2, 5
wd = [ff_weights(x[i], eps) for i in range(n)]
print("  FF check:", all(check_ff(w) for w in wd))
for topset, botset in [((0,1),(2,3)), ((0,1),(2,4)), ((0,1),(3,4)), ((0,2),(3,4))]:
    top = [1 if j in topset else 0 for j in range(m)]
    bot = [1 if j in botset else 0 for j in range(m)]
    z = Z(n, m, top, bot, [0]*n, [0]*n, wd)
    print("  top=%s bot=%s :" % (topset, botset))
    print("      Z =", sp.expand(z))
