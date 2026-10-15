"""Control B: subset formula  vs  contour-integral formula, both at c -> infinity.
Two independently-derived printed formulas for the same object (`G nu sum expression thm eq`
line 3386 vs `eq: G nu integral forumla thm-cinf` line 3850).  The integral is evaluated as a
sum of residues at w_i = x_{j_i} over INJECTIVE tuples (the paper's own prescription:
"only the residues at w_i=x_j should be evaluated", l.3802, and the (w_j-w_i) numerators make
repeated choices vanish, l.3821)."""
import engine as E, sympy as sp, random, itertools, sys
from sympy import Rational as R
q, a = E.q, E.a

def integral_formula(kappa, L, subs, sign_plant=False):
    n = len(kappa)
    xv = [subs[E.xsym(i)] for i in range(1, L+1)]
    ws = [sp.Symbol('w%d' % i) for i in range(1, n+1)]
    pref = sp.Integer(1)
    for xi in xv:
        pref *= xi*(1-subs[a]*xi)/(xi-subs[a])
    if n == 0:
        return pref
    I = sp.Integer(1)
    for i in range(n):
        for j in range(L):
            I *= (subs[q]*ws[i]-xv[j])/(ws[i]-xv[j]) * (1-ws[i]*xv[j])/(1-subs[q]*ws[i]*xv[j])
    for i in range(n):
        for j in range(i+1, n):
            num = ws[j]-ws[i]
            if sign_plant:
                num = ws[j]+ws[i]       # planted error, mirroring engine.sign_plant
            I *= num/(subs[q]*ws[j]-ws[i]) * (1-subs[q]*ws[i]*ws[j])/(1-ws[i]*ws[j])
    for i in range(n):
        m = kappa[i]
        ym = subs[E.ysym(m)]
        I *= subs[a]*(1-subs[q]*ws[i]**2)/(ws[i]*(1-subs[a]*ws[i])) * ym/(1-subs[q]*ws[i]*ym)
        for j in range(1, m):
            I *= (1-ws[i]*subs[E.ysym(j)])/(1-subs[q]*ws[i]*subs[E.ysym(j)])
    tot = sp.Integer(0)
    for tup in itertools.permutations(range(L), n):
        Rf = I
        for i in range(n):
            Rf = sp.cancel((ws[i]-xv[tup[i]])*Rf)
            Rf = Rf.subs(ws[i], xv[tup[i]])
        tot += sp.cancel(Rf)
    return sp.cancel(pref*tot)

def rnd(L, maxy, seed):
    rr = random.Random(seed)
    s = {q: R(rr.randint(3,29), 97), a: R(rr.randint(31,59), 101)}
    for j in range(1, maxy+1):
        s[E.ysym(j)] = R(rr.randint(3,29), 103+j)
    for i in range(1, L+1):
        s[E.xsym(i)] = R(rr.randint(7,61), 107+i)
    return s

cases = [((1,),1), ((1,),2), ((2,),2), ((2,1),2), ((3,1),3), ((2,1),3), ((3,2,1),3)]
for plant in (False, True):
    print("=== %s ===" % ("PLANTED sign error" if plant else "clean"))
    for kap, L in cases:
        ok = []
        for sd in (0,1):
            s = rnd(L, max(kap), sd)
            sub = E.G(kap, L, sign_plant=plant).subs(s)
            itg = integral_formula(kap, L, s, sign_plant=plant)
            ok.append(sp.nsimplify(sp.cancel(sub-itg)))
        print("  kappa=%-10s L=%d  residuals %s" % (str(kap), L, ok))
        sys.stdout.flush()
