"""
Brute-force check of Rick's (L1)-(L4) (Day 204 S4a).

Instrument: Demazure-Lusztig operators
    T_i f = t f + (t x_i - x_{i+1})/(x_i - x_{i+1}) * (s_i f - f)
which satisfy (T_i - t)(T_i + 1) = 0 and T_i f = t f for f symmetric in x_i, x_{i+1}.

sigma_m = sum_{k=0}^{m-1} T_k T_{k-1} ... T_1   (k=0 term is the identity)

We check three things:
  (A) the transposition formula
        sigma_m f = sum_i (s_{1i} f) * prod_{j != i} (x_i - t x_j)/(x_i - x_j)
      for tail-symmetric f  -- this is the engine of the proof;
  (B) Rick's (L2),(L4) -- these are the UNTUNED POSITIVE CONTROL. If the
      instrument fails these, the instrument is broken and nothing else it
      says counts;
  (C) Rick's (L1),(L3) and Clio's corrected (L1),(L3).
"""
import sys
from sympy import symbols, simplify, cancel, together, expand, Poly, Rational, factor

t = symbols('t')

def make(m):
    return symbols('x1:%d' % (m+1))

def s_act(f, i, X):
    """swap x_i, x_{i+1}  (1-indexed i)"""
    a, b = X[i-1], X[i]
    tmp = symbols('TMP_')
    return f.subs({a: tmp, b: a}).subs({tmp: b})

def T(f, i, X):
    a, b = X[i-1], X[i]
    return cancel(t*f + (t*a - b)/(a - b)*(s_act(f, i, X) - f))

def sigma(f, X):
    m = len(X)
    total = 0
    for k in range(0, m):
        g = f
        for i in range(1, k+1):          # T_k...T_1 f  =  apply T_1 first, then T_2, ...
            g = T(g, i, X)
        total = total + g
    return cancel(together(total))

def e(j, vars_):
    """elementary symmetric polynomial e_j in the given variables"""
    from sympy import symmetric_poly
    if j < 0: return 0
    if j == 0: return 1
    if j > len(vars_): return 0
    return symmetric_poly(j, *vars_)

def tail(X, i):
    return [x for k, x in enumerate(X) if k != i]

def transposition_formula(a, mu, X):
    """sum_i x_i^a * e_mu(hat x_i) * prod_{j!=i} (x_i - t x_j)/(x_i - x_j)"""
    m = len(X)
    tot = 0
    for i in range(m):
        xi = X[i]
        rest = tail(X, i)
        term = xi**a
        for p in mu:
            term *= e(p, rest)
        pr = 1
        for xj in rest:
            pr *= (xi - t*xj)/(xi - xj)
        tot += term*pr
    return cancel(together(tot))

def lhs_input(a, mu, X):
    """x_1^a * prod_p e_p(x_2..x_m)"""
    rest = list(X[1:])
    f = X[0]**a
    for p in mu:
        f *= e(p, rest)
    return expand(f)

def qint(n):
    """[n]_t"""
    if n <= 0: return 0
    return sum(t**i for i in range(n))

def report(name, lhs, rhs):
    d = cancel(expand(lhs - rhs))
    ok = (simplify(d) == 0)
    print("   %-28s %s" % (name, "MATCH" if ok else "DIFFERS"))
    return ok, d

def run(m, r):
    X = make(m)
    allv = list(X)
    E = lambda j: e(j, allv)
    print("=== m=%d, r=%d ===" % (m, r))

    cases = {
        'L1': (1, [r, 1]),
        'L2': (2, [r]),
        'L3': (2, [r-1, 1]),
        'L4': (3, [r-1]),
    }
    rick = {
      'L1':  qint(r+2)*E(r+2) + t*qint(r)*E(r+1)*E(1),
      'L2': -qint(r+2)*E(r+2) + E(r+1)*E(1),
      'L3': -qint(r+2)*E(r+2) - t*qint(r)*E(r+1)*E(1) + qint(2)*E(r)*E(2),
      'L4':  qint(r+2)*E(r+2) - E(r+1)*E(1) - qint(2)*E(r)*E(2) + E(r)*E(1)**2,
    }
    clio = {
      'L1':  qint(r+2)*E(r+2) + qint(r+1)*E(r+1)*E(1),
      'L2': -qint(r+2)*E(r+2) + E(r+1)*E(1),
      'L3': -qint(r+2)*E(r+2) - t*qint(r)*E(r+1)*E(1) + E(r)*E(1)**2,
      'L4':  qint(r+2)*E(r+2) - E(r+1)*E(1) - qint(2)*E(r)*E(2) + E(r)*E(1)**2,
    }

    out = {}
    for nm, (a, mu) in cases.items():
        f = lhs_input(a, mu, X)
        S = sigma(f, X)
        # (A) engine check
        TF = transposition_formula(a, mu, X)
        agree_engine = simplify(cancel(expand(S - TF))) == 0
        print(" %s: sigma_m (brute Hecke)  vs  transposition formula : %s"
              % (nm, "MATCH" if agree_engine else "DIFFERS"))
        okr, _ = report("vs Rick's RHS", S, rick[nm])
        okc, _ = report("vs Clio's RHS", S, clio[nm])
        out[nm] = (agree_engine, okr, okc)
    return out

if __name__ == '__main__':
    m = int(sys.argv[1]); r = int(sys.argv[2])
    run(m, r)
