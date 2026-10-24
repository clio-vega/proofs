"""T3. THE DECISIVE INSTRUMENT.  Factor the numerator of Psi^C over QQ(w_1..w_n).

P1 (pre-registered in the brief): for n > r+1 there is NO cancellation locus, so the
zero locus is exactly the structural one.  The brief's own test: if every irreducible
factor is (w_a w_b - 1) or (w_a^2 - 1), P1 confirmed; any other factor refutes it and
IS the new locus in closed form.
"""
import sympy as sp
from itertools import permutations, combinations

def psi_sym(w, theta):
    n = len(w); r = len(theta)
    m = [theta[j] + j for j in range(r)]
    total = sp.Integer(0)
    for T in combinations(range(n), r):
        x = [w[i] for i in T]
        Bt = sp.Integer(1)
        for xi in x: Bt *= (1 - 1/xi**2)
        for a in range(r):
            for b in range(a+1, r): Bt *= (1 - 1/(x[a]*x[b]))
        for t in range(n):
            if t in T: continue
            for xi in x: Bt *= (1 - 1/(w[t]*xi))/(1 - w[t]/xi)
        Mx = sp.Matrix(r, r, lambda i, j: x[i]**m[j])
        V = sp.Integer(1)
        for i in range(r):
            for k in range(i+1, r): V *= (x[k] - x[i])
        total += Bt * Mx.det()/V
    return sp.together(sp.simplify(total))

def recip_shape(fac, w):
    """is fac (up to unit) of the form w_a w_b - 1 (a<=b)?  i.e. a structural factor"""
    n = len(w)
    for a in range(n):
        for b in range(a, n):
            q = sp.simplify(sp.expand(fac) - sp.expand(w[a]*w[b] - 1))
            if q == 0: return f"w{a+1}w{b+1}-1"
            q = sp.simplify(sp.expand(fac) + sp.expand(w[a]*w[b] - 1))
            if q == 0: return f"-(w{a+1}w{b+1}-1)"
    return None

cases = []
for n in range(1, 5):
    for r in range(1, min(n, 2)+1):
        for theta in ([tuple(range(r)), tuple(1+j for j in range(r)), (0, 3)[:r] if r == 2 else (2,)]):
            if len(theta) != r: continue
            m = [theta[j]+j for j in range(r)]
            if len(set(m)) != r: continue
            cases.append((n, r, theta))

seen = set()
for (n, r, theta) in cases:
    key = (n, r, theta)
    if key in seen: continue
    seen.add(key)
    w = sp.symbols(f'w1:{n+1}', nonzero=True)
    P = psi_sym(list(w), list(theta))
    num, den = sp.fraction(sp.cancel(P))
    num = sp.factor(num)
    m = [theta[j]+j for j in range(r)]
    print("="*78)
    print(f"n={n} r={r} theta={theta}  m={m}")
    print(f"  numerator factored: {num}")
    facs = sp.factor_list(num)
    unit, flist = facs[0], facs[1]
    verdict = []
    for f, e in flist:
        if f.free_symbols and f.is_Symbol:
            verdict.append((str(f), e, "monomial (unit on U)")); continue
        tag = recip_shape(f, list(w))
        verdict.append((str(f), e, tag if tag else "*** NOT OF RECIPROCAL SHAPE ***"))
    for f, e, t in verdict:
        print(f"    [{e}]  {f:40}  {t}")
