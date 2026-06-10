"""
fiber_engine.py — exact G_lambda(zeta_d) via the PROVED master identity, plus
principal-specialisation machinery for the RSW probe (Job C).

Master identity (proved, 2026-06-02-zeta3-nonvanishing.tex, Prop. "Master identity"):
  For lambda |- n,  m = floor(n/2),  e = n mod 2,
      A = {1,3,...,n-1}  if n even ;   {2,4,...,n-1}  if n odd,
  and ANY complex zeta,
      G_lambda(zeta) = < s_lambda,  h_1^e  prod_{k in A} (h_2 + e_2 * zeta^{2k-1}) >.
  Specialisations:  zeta=-1 -> chi^lambda(2^m 1^e);  zeta=i -> <s_lambda, (h_2+i e_2)^m>.

We carry zeta as a sympy symbol t, build the product ONCE per n (it is lambda-independent),
pair with each s_lambda to get P_lambda(t) in Z[t], then evaluate at exact roots of unity.

Principal specialisation s_lambda(1,q,...,q^{N-1}) via the q-hook-content product
  s_lambda(1,q,..,q^{N-1}) = q^{n(lambda)} prod_{cells} [N + c(i,j)]_q / [h(i,j)]_q,
  c(i,j) = j-i (content),  h = hook length,  [k]_q = (1-q^k)/(1-q),  n(lambda)=sum (i-1) lambda_i.
Computed as an exact polynomial in q, then evaluated at q = zeta_d.

All exact (sympy).  Cross-checked against job1_tie_census.G_gaussian for d=4.
"""
import sympy as sp
from functools import lru_cache
import symfunc as sf
from job1_tie_census import partitions as _parts_gen

t = sp.Symbol('t')

# ---------------------------------------------------------------- partitions
def partitions(n):
    return list(_parts_gen(n))

# ---------------------------------------------------------------- master-id product
@lru_cache(maxsize=None)
def fiber_product(n):
    """h_1^e * prod_{k in A} (h_2 + e_2 t^{2k-1})  in the power-sum basis (coeffs in Z[t])."""
    m = n // 2
    e = n % 2
    A = list(range(1, n, 2)) if n % 2 == 0 else list(range(2, n, 2))
    assert len(A) == m, (n, A, m)
    h2 = sf.h_n(2)
    e2 = sf.e_n(2)
    prod = sf.power(sf.h1(), e) if e else {(): sp.Integer(1)}
    for k in A:
        factor = sf.add(h2, sf.scale(e2, t**(2*k - 1)))
        prod = sf.mul(prod, factor)
    return prod

@lru_cache(maxsize=None)
def G_poly(lam):
    """P_lambda(t) = G_lambda(zeta) as an exact polynomial in t = zeta.  lam tuple."""
    n = sum(lam)
    prod = fiber_product(n)
    P = sf.inner(sf.schur_p(lam), prod)     # sympy expr in t
    return sp.expand(P)

# zeta_d exact values (primitive d-th root e^{2 pi i / d})
def zeta(d):
    if d == 1: return sp.Integer(1)
    if d == 2: return sp.Integer(-1)
    if d == 3: return sp.Rational(-1, 2) + sp.I*sp.sqrt(3)/2
    if d == 4: return sp.I
    if d == 6: return sp.Rational(1, 2) + sp.I*sp.sqrt(3)/2
    return sp.exp(2*sp.I*sp.pi/d)

def G_value(lam, d):
    """Exact G_lambda(zeta_d) as a sympy complex number (a + b i form when possible)."""
    P = G_poly(tuple(lam))
    return sp.simplify(sp.expand(P.subs(t, zeta(d))))

# ---------------------------------------------------------------- principal specialisation
def hook_lengths(lam):
    lam = [x for x in lam if x > 0]
    conj = [sum(1 for r in lam if r > c) for c in range(lam[0])] if lam else []
    hooks = {}
    for i, row in enumerate(lam):
        for j in range(row):
            arm = row - 1 - j
            leg = conj[j] - 1 - i
            hooks[(i, j)] = arm + leg + 1
    return hooks

def n_stat(lam):
    return sum(i*lam[i] for i in range(len(lam)))   # sum (i-1)*lam_i with 0-indexing -> i*lam_i? careful

def n_lambda(lam):
    # n(lambda) = sum_{i>=1} (i-1) lambda_i  (1-indexed); 0-indexed: sum i*lam_i
    return sum(i*part for i, part in enumerate(lam))

def qint(k, q):
    """[k]_q = 1 + q + ... + q^{k-1} as exact sympy."""
    if k == 0:
        return sp.Integer(0)
    return sum(q**i for i in range(k))

@lru_cache(maxsize=None)
def princ_spec_poly(lam, N):
    """s_lambda(1,q,...,q^{N-1}) as an exact polynomial in q (sympy).  0 if N < len(lam)."""
    lam = tuple(x for x in lam if x > 0)
    if N < len(lam):
        return sp.Integer(0)
    q = sp.Symbol('q')
    hooks = hook_lengths(lam)
    num = sp.Integer(1)
    den = sp.Integer(1)
    for (i, j), h in hooks.items():
        c = j - i
        num *= qint(N + c, q)      # [N + content]_q
        den *= qint(h, q)          # [hook]_q
    expr = q**n_lambda(lam) * num / den
    # the q-hook-content quotient is a genuine polynomial; cancel to kill denominators
    # (else it is 0/0 at roots of unity).  sp.cancel reduces to lowest terms.
    poly = sp.cancel(expr)
    assert sp.denom(poly) == 1, f"princ-spec not a polynomial: lam={lam} N={N}"
    return sp.Poly(poly, q).as_expr(), q

def princ_spec_value(lam, N, d):
    """s_lambda(1,..,q^{N-1}) | q = zeta_d, exact."""
    res = princ_spec_poly(tuple(lam), N)
    if res == 0 or (isinstance(res, sp.Integer) and res == 0):
        return sp.Integer(0)
    poly, q = res
    return sp.simplify(poly.subs(q, zeta(d)))

# ---------------------------------------------------------------- self-test
if __name__ == "__main__":
    from job1_tie_census import G_gaussian, M_vector, mn_character
    # 1) master identity reproduces d=4 G_gaussian (re,im) for all lam |- 2m, m<=5
    print("=== cross-check: master-id G vs G_gaussian (d=4) ===")
    bad = 0; tot = 0
    for m in range(1, 6):
        for lam in partitions(2*m):
            Ms = M_vector(lam, m)
            re, im = G_gaussian(lam, m, Ms)
            gv = G_value(lam, 4)
            gv = sp.expand(gv)
            re2 = sp.re(gv); im2 = sp.im(gv)
            tot += 1
            if sp.simplify(re2 - re) != 0 or sp.simplify(im2 - im) != 0:
                bad += 1
                if bad <= 5:
                    print(f"  MISMATCH lam={lam}: master=({re2},{im2}) gauss=({re},{im})")
    print(f"  d=4: {tot-bad}/{tot} match")

    # 2) d=2: G_lambda(-1) == chi^lambda(2^m 1^e)
    print("=== cross-check: G(-1) == chi^lambda(2^m 1^e) ===")
    bad = 0; tot = 0
    for n in range(1, 11):
        m = n // 2; e = n % 2
        cl = tuple([2]*m + [1]*e)
        for lam in partitions(n):
            g = sp.simplify(G_value(lam, 2))
            chi = mn_character(tuple(lam), cl)
            tot += 1
            if sp.simplify(g - chi) != 0:
                bad += 1
                if bad <= 5:
                    print(f"  MISMATCH lam={lam}: G(-1)={g} chi={chi}")
    print(f"  d=2: {tot-bad}/{tot} match")

    # 3) d=3: vanishing set should be exactly {(3,1),(2,1,1)} among |lam|<=8 (proved)
    print("=== d=3 vanishing set (should be (3,1),(2,1,1) only) ===")
    vanish = []
    for n in range(1, 9):
        for lam in partitions(n):
            g = sp.simplify(G_value(lam, 3))
            if g == 0:
                vanish.append(lam)
    print("  G_lambda(zeta_3)=0 for:", vanish)

    # 4) principal-spec sanity: s_lambda(1,..,1) (N vars, q->1) = dim of GL_N irrep
    print("=== princ-spec sanity: s_lambda(1,q,..,q^{N-1}) at q=1 = #SSYT ===")
    q = sp.Symbol('q')
    for lam in [(2,1), (3,1), (2,2)]:
        poly, qq = princ_spec_poly(lam, 4)
        val1 = poly.subs(qq, 1)
        print(f"  lam={lam} N=4: s(1^4)={val1}")
