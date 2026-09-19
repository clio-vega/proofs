"""Box theorem on the FREE-FERMION side, reproducing ff_enum.py's 352 exactly,
plus a strengthened five-vertex Schur control.

Two questions:
  (A) Does M-convex <=> A-is-a-box hold on the free-fermion families too?
  (B) Q186's premise: is the free-fermion condition CORRELATED with M-convexity
      failure?  Compare the FF failure rate against the non-FF rate (328/769).
"""
import sympy as sp, itertools
from fractions import Fraction
from sixvertex import Z, schur_weights
from lorentzian import is_lorentzian, is_M_convex, multi_factorial

X, Y = sp.symbols('X Y')
CO = [(p,q) for p in range(3) for q in range(3) if (p,q) != (0,0)]
def lf(pq, xi, yi): return pq[0]*xi + pq[1]*yi

def is_box(A):
    A = set(A); n = len(next(iter(A)))
    lo = [min(a[i] for a in A) for i in range(n)]
    hi = [max(a[i] for a in A) for i in range(n)]
    size = 1
    for i in range(n): size *= hi[i] - lo[i] + 1
    return size == len(A)

def cdict(expr, vars_):
    P = sp.Poly(sp.expand(expr), *vars_)
    return {tuple(m): Fraction(int(c)) for m, c in zip(P.monoms(), P.coeffs())}

fams = []
for pat in itertools.product(CO, repeat=6):
    a1,a2,b1,b2,c1,c2 = [lf(p, X, Y) for p in pat]
    if sp.expand(a1*a2 + b1*b2 - c1*c2) == 0:
        fams.append(pat)
six = [p for p in fams if all(q != (0,0) for q in p)]
print("free-fermion families, all six weights nonzero: %d" % len(six))

tot=mono=nontriv=notMconv=disagree=xlev1=0
notLor=mconv_fail=eig_fail=0
sigma_const_passers=0
for pat in six[:40]:
    for n, m in [(2,2), (2,3)]:
        xs = sp.symbols('x1:%d'%(n+1)); ys = sp.symbols('y1:%d'%(n+1))
        w = [{k: lf(p, xs[i], ys[i]) for k,p in zip(['a1','a2','b1','b2','c1','c2'], pat)}
             for i in range(n)]
        vars_ = list(xs) + list(ys)
        for top in itertools.product([0,1], repeat=m):
            for bot in itertools.product([0,1], repeat=m):
                z = Z(n, m, list(top), list(bot), [0]*n, [0]*n, w)
                if z == 0: continue
                cd = cdict(z, vars_); tot += 1
                for a in cd:
                    for i in range(n): assert a[i]+a[n+i] == m
                if len(cd) == 1: mono += 1; continue
                nontriv += 1
                mc = is_M_convex(set(cd))
                A = {tuple(a[:n]) for a in cd}
                bx = is_box(A)
                if mc != bx:
                    disagree += 1
                    if disagree <= 5: print("  MISMATCH", pat, n, m, sorted(A), mc, bx)
                if not mc: notMconv += 1
                if len({sum(t) for t in A}) == 1:
                    xlev1 += 1
                    if mc: sigma_const_passers += 1
                nd = {a: c/multi_factorial(a) for a,c in cd.items()}
                ok, why = is_lorentzian(nd, len(vars_))
                if not ok:
                    notLor += 1
                    if '(L2)' in why: mconv_fail += 1
                    else: eig_fail += 1

print("\n--- FREE-FERMION side (reproducing ff_enum.py) ---")
print("  nonzero Z / monomial / NON-MONOMIAL : %d / %d / %d" % (tot, mono, nontriv))
print("  not Lorentzian : %d  (M-convexity %d, eigenvalue %d)" % (notLor, mconv_fail, eig_fail))
print("  not M-convex   : %d" % notMconv)
print("  M-convex <=> A is a box : DISAGREE %d / %d   <-- must be 0" % (disagree, nontriv))
print("  constant Sigma : %d  (of which M-convex: %d, (T2) refuted if >0)" % (xlev1, sigma_const_passers))
print("\n--- (B) Q186's premise: does free-fermion CORRELATE with failure? ---")
print("  free-fermion     M-convexity failure rate : %d/%d = %.1f%%" % (notMconv, nontriv, 100.0*notMconv/nontriv))
print("  NON-free-fermion M-convexity failure rate : 328/769 = %.1f%%" % (100.0*328/769))

print("\n--- STRENGTHENED five-vertex Schur control (supports of size > 2) ---")
xv = sp.symbols('x1:7')
found = 0
for (n, m) in [(2,5),(2,6),(3,6),(3,7)]:
    ws = [schur_weights(xv[i]) for i in range(n)]
    vars_ = list(xv[:n])
    for topset in itertools.combinations(range(m), n):
        top = [1 if j in topset else 0 for j in range(m)]
        for botset in itertools.combinations(range(m), n):
            bot = [1 if j in botset else 0 for j in range(m)]
            z = Z(n, m, top, bot, [0]*n, [0]*n, ws)
            if z == 0: continue
            cd = cdict(z, vars_)
            if len(cd) < 3: continue
            mc = is_M_convex(set(cd))
            # does lem:rowdeg's conclusion hold here?  (it must NOT: no y variables)
            print("    n=%d m=%d top=%s bot=%s |supp|=%d M-convex=%s  Z=%s"
                  % (n,m,topset,botset,len(cd),mc,sp.factor(z)))
            found += 1
            assert mc, "INSTRUMENT WRONG: HMMS Thm 2 says this must be M-convex"
            if found >= 8: break
        if found >= 8: break
    if found >= 8: break
print("    Schur controls with |supp|>=3 : %d, all M-convex" % found)
