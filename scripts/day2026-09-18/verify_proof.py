"""
Verification of every step of the proof of Rick's (L1)-(L4).

Chain:
  Step 1  T_i (x_i^c) = x_{i+1}^c + (1-t) sum_{p=1}^{c-1} x_i^p x_{i+1}^{c-p}   (c>=1)
  Step 2  rho_m(x_1^a) = x_m^a + (1-t) sum_{b=1}^{a-1} qt_b(x_1..x_{m-1}) x_m^{a-b}
          where rho_m = T_{m-1}...T_1 and qt_b = q_b/(1-t)
  Step 3  sigma_m(x_1^a) = qt_a(x_1..x_m)  for a>=1;  sigma_m(1) = [m]_t
  Step 4  Lambda_m-linearity: sigma_m(g f) = g sigma_m(f) for g symmetric
  Step 5  Master Lemma  M(a,k) := sigma_m[x_1^a e_k(hat x_1)] = sum_l (-1)^l e_{k-l} qt_{a+l}
  Step 6  closed forms of M(1,k), M(2,k), M(3,k)
  Step 7  Product Lemma  sigma_m[x_1^a e_k(hat) e_1(hat)] = e_1 M(a,k) - M(a+1,k)
  Step 8  (L1)-(L4)
"""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/scripts/day2026-09-18')
from sigma_bruteforce import T, sigma, s_act, e, make, tail, qint, t, transposition_formula
from sympy import symbols, simplify, cancel, expand, series, Poly, together

def qtilde(n, vars_):
    """qt_n = q_n/(1-t), where sum_n q_n u^n = prod_j (1-t x_j u)/(1 - x_j u).
       Computed as a truncated power series coefficient."""
    if n == 0: return 1
    u = symbols('u')
    P = 1
    for x in vars_:
        P *= (1 - t*x*u)
    D = 1
    for x in vars_:
        D *= (1 - x*u)
    # series expand P/D to order n
    ser = (P * sum((sum(x*u for x in vars_))**0 for _ in [0]))  # placeholder
    from sympy import series as sser
    ex = sser(P/D, u, 0, n+1).removeO()
    qn = expand(ex.coeff(u, n))
    return cancel(expand(qn/(1-t)))

def check(label, lhs, rhs):
    d = simplify(cancel(expand(lhs - rhs)))
    print("  %-52s %s" % (label, "OK" if d == 0 else "FAIL  diff=%s" % d))
    return d == 0

ok = True
print("STEP 1: T_i(x_i^c) formula")
for m in [3, 4]:
    X = make(m)
    for i in [1, 2]:
        for c in [1, 2, 3, 4]:
            lhs = T(X[i-1]**c, i, X)
            rhs = X[i]**c + (1-t)*sum(X[i-1]**p * X[i]**(c-p) for p in range(1, c))
            ok &= check("m=%d i=%d c=%d" % (m, i, c), lhs, rhs)

print("STEP 2/3: sigma_m(x_1^a) = qtilde_a   (a>=1);  sigma_m(1)=[m]_t")
for m in [2, 3, 4, 5]:
    X = make(m)
    ok &= check("m=%d a=0 -> [m]_t" % m, sigma(X[0]**0, X), qint(m))
    for a in [1, 2, 3, 4]:
        ok &= check("m=%d a=%d" % (m, a), sigma(X[0]**a, X), qtilde(a, list(X)))

print("STEP 2 alone: rho_m(x_1^a) closed form")
for m in [2, 3, 4]:
    X = make(m)
    for a in [1, 2, 3]:
        g = X[0]**a
        for i in range(1, m):
            g = T(g, i, X)                      # rho_m = T_{m-1} ... T_1
        rhs = X[m-1]**a + (1-t)*sum(qtilde(b, list(X[:m-1]))*X[m-1]**(a-b)
                                    for b in range(1, a))
        ok &= check("m=%d a=%d" % (m, a), g, rhs)

print("STEP 4: Lambda_m-linearity")
for m in [3, 4]:
    X = make(m); allv = list(X)
    g = e(2, allv)
    f = X[0]**2 * e(1, list(X[1:]))
    ok &= check("m=%d  sigma(g f) = g sigma(f)" % m, sigma(expand(g*f), X), expand(g*sigma(f, X)))

print("STEP 5/6: Master Lemma M(a,k) and its closed forms")
for m in [3, 4, 5]:
    X = make(m); allv = list(X); E = lambda j: e(j, allv)
    for a in [1, 2, 3, 4]:
        for k in [0, 1, 2]:
            lhs = sigma(expand(X[0]**a * e(k, list(X[1:]))), X)
            rhs = sum((-1)**l * E(k-l) * qtilde(a+l, allv) for l in range(0, k+1))
            ok &= check("m=%d master a=%d k=%d" % (m, a, k), lhs, expand(rhs))
    for k in [0, 1, 2, 3]:
        M1 = qint(k+1)*E(k+1)
        M2 = -qint(k+2)*E(k+2) + E(1)*E(k+1)
        M3 = qint(k+3)*E(k+3) - E(1)*E(k+2) + (E(1)**2 - qint(2)*E(2))*E(k+1)
        for a, Mc in [(1, M1), (2, M2), (3, M3)]:
            lhs = sum((-1)**l * E(k-l) * qtilde(a+l, allv) for l in range(0, k+1))
            ok &= check("m=%d closed-form a=%d k=%d" % (m, a, k), expand(lhs), expand(Mc))

print("STEP 7/8: Product Lemma and (L1)-(L4) vs Rick")
for m, rs in [(4, [1, 2]), (5, [1, 2, 3]), (6, [2, 3, 4])]:
    X = make(m); allv = list(X); E = lambda j: e(j, allv); Tl = list(X[1:])
    for r in rs:
        rick = {
          'L1': (1, [r, 1],   qint(r+2)*E(r+2) + t*qint(r)*E(r+1)*E(1)),
          'L2': (2, [r],     -qint(r+2)*E(r+2) + E(r+1)*E(1)),
          'L3': (2, [r-1,1], -qint(r+2)*E(r+2) - t*qint(r)*E(r+1)*E(1) + qint(2)*E(r)*E(2)),
          'L4': (3, [r-1],    qint(r+2)*E(r+2) - E(r+1)*E(1) - qint(2)*E(r)*E(2) + E(r)*E(1)**2),
        }
        for nm, (a, mu, rhs) in rick.items():
            f = X[0]**a
            for p in mu: f *= e(p, Tl)
            ok &= check("m=%d r=%d %s" % (m, r, nm), sigma(expand(f), X), expand(rhs))

print()
print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
