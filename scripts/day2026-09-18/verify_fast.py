"""
Fast re-verification of every step of the proof of (L1)-(L4).
T_i is computed by EXACT polynomial division rather than rational-function
cancellation; q_n from the Newton recursion for h_l rather than sympy.series.
Nothing from the proof is assumed: sigma_m composes the T_i literally.
"""
import sys
from sympy import symbols, expand, Poly, div, symmetric_poly

t = symbols('t')

def make(m): return symbols('x1:%d' % (m+1))

def s_act(f, i, X):
    a, b = X[i-1], X[i]
    tmp = symbols('TMP_')
    return f.subs({a: tmp, b: a}).subs({tmp: b})

def T(f, i, X):
    """T_i f = t f + (t x_i - x_{i+1})/(x_i - x_{i+1}) * (s_i f - f)."""
    a, b = X[i-1], X[i]
    g = expand(s_act(f, i, X) - f)
    if g == 0:
        return expand(t*f)
    gens = list(X) + [t]
    quo, rem = div(Poly(g, *gens), Poly(a - b, *gens))
    assert rem.is_zero, "non-exact division in T_%d" % i
    return expand(t*f + (t*a - b)*quo.as_expr())

def rho(f, j, X):
    g = f
    for i in range(1, j):
        g = T(g, i, X)
    return g

def sigma(f, X):
    return expand(sum(rho(f, j, X) for j in range(1, len(X)+1)))

def e(j, vars_):
    if j < 0 or j > len(vars_): return 0
    if j == 0: return 1
    return expand(symmetric_poly(j, *vars_))

def h(l, vars_):
    H = [1]
    for L in range(1, l+1):
        H.append(expand(sum((-1)**(k-1)*e(k, vars_)*H[L-k] for k in range(1, L+1))))
    return H[l]

_qc = {}
def qtilde(n, vars_):
    """q_n/(1-t) with q_n = sum_{k+l=n} (-t)^k e_k h_l"""
    if n == 0: return 1
    key = (n, tuple(vars_))
    if key in _qc: return _qc[key]
    qn = expand(sum((-t)**k * e(k, vars_) * h(n-k, vars_) for k in range(0, n+1)))
    gens = list(vars_) + [t]
    quo, rem = div(Poly(qn, *gens), Poly(1-t, *gens))
    assert rem.is_zero, "q_%d not divisible by (1-t)" % n
    _qc[key] = expand(quo.as_expr())
    return _qc[key]

def qint(n):
    return 0 if n <= 0 else sum(t**i for i in range(n))

nok = nfail = 0
def check(label, lhs, rhs):
    global nok, nfail
    d = expand(lhs - rhs)
    if d == 0:
        nok += 1; print("  OK    %s" % label)
    else:
        nfail += 1; print("  FAIL  %s   diff=%s" % (label, d))

print("STEP 1: T_i(x_i^c) = x_{i+1}^c + (1-t) sum_{p=1}^{c-1} x_i^p x_{i+1}^{c-p}")
for m in [3,4]:
    X = make(m)
    for i in range(1, m):
        for c in [1,2,3,4,5]:
            check("m=%d i=%d c=%d" % (m,i,c), T(X[i-1]**c, i, X),
                  expand(X[i]**c + (1-t)*sum(X[i-1]**p*X[i]**(c-p) for p in range(1,c))))

print("STEP 2 (Prop 2.6, eq (rho)): rho_j(x_1^a) closed form")
for m in [2,3,4,5]:
    X = make(m)
    for j in range(1, m+1):
        for a in [1,2,3,4]:
            check("m=%d j=%d a=%d" % (m,j,a), rho(X[0]**a, j, X),
                  expand(X[j-1]**a + (1-t)*sum(qtilde(b, list(X[:j-1]))*X[j-1]**(a-b)
                                               for b in range(1,a))))

print("STEP 3: sigma_m(x_1^a) = qtilde_a  (a>=1),  sigma_m(1)=[m]_t")
for m in [2,3,4,5,6]:
    X = make(m)
    check("m=%d a=0" % m, sigma(X[0]**0, X), qint(m))
    for a in [1,2,3,4,5]:
        check("m=%d a=%d" % (m,a), sigma(X[0]**a, X), qtilde(a, list(X)))

print("STEP 4: Lambda_m-linearity sigma(g f) = g sigma(f)")
for m in [3,4,5]:
    X = make(m); A = list(X)
    for gi in [1,2]:
        g = e(gi, A)
        f = expand(X[0]**2 * e(1, list(X[1:])))
        check("m=%d g=e_%d" % (m,gi), sigma(expand(g*f), X), expand(g*sigma(f, X)))

print("STEP 5: Master Lemma  sigma[x_1^a e_k(tail)] = sum_l (-1)^l e_{k-l} qtilde_{a+l}")
for m in [3,4,5]:
    X = make(m); A = list(X); Tl = list(X[1:]); E = lambda j: e(j, A)
    for a in [1,2,3,4]:
        for k in [0,1,2,3]:
            check("m=%d a=%d k=%d" % (m,a,k),
                  sigma(expand(X[0]**a * e(k, Tl)), X),
                  expand(sum((-1)**l*E(k-l)*qtilde(a+l, A) for l in range(0,k+1))))

print("STEP 6: closed forms M(1,k), M(2,k), M(3,k)")
for m in [3,4,5]:
    X = make(m); A = list(X); E = lambda j: e(j, A)
    for k in [0,1,2,3,4]:
        M = {1: qint(k+1)*E(k+1),
             2: -qint(k+2)*E(k+2) + E(1)*E(k+1),
             3: qint(k+3)*E(k+3) - E(1)*E(k+2) + (E(1)**2 - qint(2)*E(2))*E(k+1)}
        for a in [1,2,3]:
            check("m=%d closed a=%d k=%d" % (m,a,k),
                  expand(sum((-1)**l*E(k-l)*qtilde(a+l, A) for l in range(0,k+1))),
                  expand(M[a]))

print("STEP 7: Product Lemma  sigma[x_1^a e_k(tail) e_1(tail)] = e_1 M(a,k) - M(a+1,k)")
for m in [3,4,5]:
    X = make(m); A = list(X); Tl = list(X[1:]); E = lambda j: e(j, A)
    Mf = lambda a,k: expand(sum((-1)**l*E(k-l)*qtilde(a+l, A) for l in range(0,k+1)))
    for a in [1,2]:
        for k in [0,1,2,3]:
            check("m=%d a=%d k=%d" % (m,a,k),
                  sigma(expand(X[0]**a * e(k,Tl) * e(1,Tl)), X),
                  expand(E(1)*Mf(a,k) - Mf(a+1,k)))

print("STEP 8: (L1)-(L4) exactly as Rick states them")
for m, rs in [(3,[1,2]), (4,[1,2,3]), (5,[1,2,3,4])]:
    X = make(m); A = list(X); Tl = list(X[1:]); E = lambda j: e(j, A)
    for r in rs:
        tag = "m=%d r=%d%s" % (m, r, "  [degenerate m<r+2]" if m < r+2 else "")
        spec = {
          'L1': (1, [r,1],    qint(r+2)*E(r+2) + t*qint(r)*E(r+1)*E(1)),
          'L2': (2, [r],     -qint(r+2)*E(r+2) + E(r+1)*E(1)),
          'L3': (2, [r-1,1], -qint(r+2)*E(r+2) - t*qint(r)*E(r+1)*E(1) + qint(2)*E(r)*E(2)),
          'L4': (3, [r-1],    qint(r+2)*E(r+2) - E(r+1)*E(1) - qint(2)*E(r)*E(2) + E(r)*E(1)**2),
        }
        for nm in ['L1','L2','L3','L4']:
            a, mu, rhs = spec[nm]
            f = X[0]**a
            for p in mu: f = expand(f*e(p,Tl))
            check("%s %s" % (tag, nm), sigma(f, X), expand(rhs))

print()
print("TOTAL: %d passed, %d failed" % (nok, nfail))
