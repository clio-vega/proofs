import engine as E, sympy as sp, random
from sympy import Rational as R
q,a = E.q, E.a

def rnd_sub(L, maxy, seed):
    rr = random.Random(seed)
    s = {q: R(rr.randint(2,40), rr.randint(41,80)), a: R(rr.randint(2,40), rr.randint(41,80))}
    for j in range(1, maxy+1):
        s[E.ysym(j)] = R(rr.randint(2,40), rr.randint(41,80))
    for i in range(1, L+1):
        s[E.xsym(i)] = R(rr.randint(2,40), rr.randint(41,80))
    return s

print("=== CONTROL A: symmetry of G_kappa in x_1..x_L (NOT manifest in the subset formula) ===")
for kap, L in [((1,),2), ((2,),2), ((2,1),2), ((1,),3), ((3,1),3), ((3,2,1),3)]:
    for plant in (False, True):
        Gk = E.G(kap, L, sign_plant=plant)
        worst = 0
        for sd in range(3):
            s = rnd_sub(L, max(kap) if kap else 1, sd)
            v0 = sp.nsimplify(Gk.subs(s))
            # swap x1 <-> x2
            s2 = dict(s); s2[E.xsym(1)], s2[E.xsym(2)] = s[E.xsym(2)], s[E.xsym(1)]
            v1 = sp.nsimplify(Gk.subs(s2))
            worst = max(worst, abs(sp.simplify(v0-v1)))
        tag = "PLANTED" if plant else "clean  "
        print("  %s kappa=%-9s L=%d  max|G(x)-G(swap)| = %s" % (tag, str(kap), L, worst))
