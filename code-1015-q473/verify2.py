import engine as E, sympy as sp, random, sys, itertools
from sympy import Rational as R
q,a = E.q, E.a
rr = random.Random(7)
s = {q: R(5,97), a: R(37,101)}
for j in range(1,5): s[E.ysym(j)] = R(rr.randint(3,29), 103+j)

print("=== CONTROL C: stability lemma  G_kappa(x_1..x_{L-1},1) = G_kappa(x_1..x_{L-1}) ===")
for kap in [(), (1,), (2,), (2,1), (3,1), (3,2,1)]:
    for L in range(max(1,len(kap)), 4):
        if L < len(kap): continue
        lhs = E.G(kap, L).subs({E.xsym(L): 1})
        rhs = E.G(kap, L-1) if L-1 >= len(kap) else sp.Integer(0)
        r = sp.simplify(sp.cancel(sp.together(lhs - rhs)))
        print("   kappa=%-10s L=%d -> L-1=%d   residual %s" % (str(kap), L, L-1, r))
    sys.stdout.flush()

print("=== CONTROL D: G_kappa = 0 when |kappa| > L  (literal subset formula = empty sum) ===")
# engine asserts L>=d; emulate the literal empty sum
def G_literal(kappa, L):
    d=len(kappa)
    if L < d: return sp.Integer(0)
    return E.G(kappa,L)
for kap,L in [((2,1),1), ((3,2,1),2), ((3,2,1),1)]:
    print("   kappa=%s L=%d -> %s" % (str(kap), L, G_literal(kap,L)))

print("=== CONTROL E: rank of the G-matrix, exact rational AND mod p (two engines) ===")
p = 2147483647
for L in (2,3):
    kaps = [k for k in E.subsets_upto(4,L) if len(k)<=L]
    Gs = [sp.cancel(E.G(k,L).subs(s)) for k in kaps]
    # evaluate at many random x-points
    pts=[]
    for t in range(len(kaps)+4):
        pts.append({E.xsym(i): R(rr.randint(7,400), 1009+7*t+i) for i in range(1,L+1)})
    M = sp.Matrix([[sp.Rational(sp.cancel(g.subs(pt))) for g in Gs] for pt in pts])
    rk_exact = M.rank()
    Mp = sp.Matrix([[sp.Integer(M[i,j].p) * pow(int(M[i,j].q), p-2, p) % p
                     for j in range(M.cols)] for i in range(M.rows)])
    rk_modp = Mp.rank(iszerofunc=lambda z: sp.Integer(z) % p == 0)
    # planted dependency: append column G_{k1}+G_{k2}
    Mdep = M.row_join(sp.Matrix([[M[i,0]+M[i,1]] for i in range(M.rows)]))
    print("   L=%d  #kappa=%d   rank_exact=%d  rank_modp=%d   rank_with_planted_dep=%d (must equal rank_exact)"
          % (L, len(kaps), rk_exact, rk_modp, Mdep.rank()))
    sys.stdout.flush()
