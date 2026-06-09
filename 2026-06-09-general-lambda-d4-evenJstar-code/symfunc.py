"""
Exact symmetric-function machinery in sympy (no Sage).

Strategy:
- Represent symmetric functions in the POWER-SUM basis p_mu as dicts {mu(tuple sorted desc): coeff}.
- Power sums are orthogonal: <p_mu, p_nu> = z_mu * delta. So inner products are exact.
- Multiplication of p-basis elements is trivial (concatenate partitions).
- Schur s_lambda in p-basis:  s_lambda = sum_mu (chi^lambda(mu)/z_mu) p_mu,
  with chi^lambda(mu) computed by Murnaghan-Nakayama (exact integers).
- h_n = sum_{mu |- n} p_mu / z_mu ; e_n = sum_{mu|-n} eps(mu) p_mu/z_mu  (eps = sign).
  Actually h_n = sum_{mu} p_mu/z_mu, e_n = sum_{mu} (-1)^{n-len(mu)} p_mu/z_mu.

All coeffs are sympy Rationals / Gaussian integers as needed.
"""
from sympy import Rational, factorial, I, Integer
from functools import lru_cache
from itertools import product as iproduct

# ---------- partitions ----------
def partitions(n):
    """All partitions of n as tuples sorted descending."""
    if n == 0:
        yield ()
        return
    def helper(n, maxpart):
        if n == 0:
            yield ()
            return
        for k in range(min(n, maxpart), 0, -1):
            for rest in helper(n-k, k):
                yield (k,) + rest
    yield from helper(n, n)

def z_mu(mu):
    """z_mu = prod_i i^{m_i} m_i!  where m_i = #parts equal to i."""
    from collections import Counter
    c = Counter(mu)
    z = Integer(1)
    for part, mult in c.items():
        z *= Integer(part)**mult * factorial(mult)
    return z

def sign_mu(mu):
    """epsilon_mu = (-1)^{n - len(mu)} = prod (-1)^{part-1}."""
    n = sum(mu)
    return (-1)**(n - len(mu))

# ---------- Murnaghan-Nakayama for chi^lambda(mu) ----------
@lru_cache(maxsize=None)
def mn_character(lam, mu):
    """
    chi^lambda(mu) via Murnaghan-Nakayama recursion.
    lam, mu are tuples (descending). Returns integer (sympy Integer).
    """
    lam = tuple(x for x in lam if x > 0)
    mu = tuple(x for x in mu if x > 0)
    if sum(lam) != sum(mu):
        return Integer(0)
    if len(mu) == 0:
        return Integer(1) if len(lam) == 0 else Integer(0)
    r = mu[0]
    rest = mu[1:]
    total = Integer(0)
    for (newlam, height) in border_strips(lam, r):
        total += Integer(-1)**height * mn_character(newlam, rest)
    return total

def border_strips(lam, r):
    """
    Yield (lam_minus_strip, height) for all border strips of size r removable from lam.
    Standard beta-set / hook approach.
    height = (#rows spanned) - 1.
    """
    lam = list(lam)
    k = len(lam)
    # beta-set: beta_i = lam_i + (k-1-i)  for i=0..k-1  (strictly decreasing)
    beta = [lam[i] + (k - 1 - i) for i in range(k)]
    beta_set = set(beta)
    results = []
    for b in beta:
        if b - r >= 0 and (b - r) not in beta_set:
            # remove a hook: move bead from b to b-r
            new_beta = sorted([x for x in beta if x != b] + [b - r], reverse=True)
            # height = number of beads strictly between (b-r) and b in original beta_set
            height = sum(1 for x in beta if (b - r) < x < b)
            # reconstruct partition from new_beta
            newlam = tuple(new_beta[i] - (k - 1 - i) for i in range(k))
            newlam = tuple(x for x in newlam if x > 0)
            results.append((newlam, height))
    return results

# ---------- representations: f in power-sum basis as dict mu->coeff ----------
def p_one(part):
    """p_part as single power sum: {(part,):1}."""
    return {(part,): Integer(1)}

def schur_p(lam):
    """s_lambda in power-sum basis: dict mu -> chi^lambda(mu)/z_mu."""
    n = sum(lam)
    d = {}
    for mu in partitions(n):
        chi = mn_character(tuple(lam), mu)
        if chi != 0:
            d[mu] = Rational(chi, 1) / z_mu(mu)
    return d

def h_n(n):
    """complete homogeneous h_n in p-basis: sum_mu p_mu/z_mu."""
    return {mu: Rational(1,1)/z_mu(mu) for mu in partitions(n)}

def e_n(n):
    """elementary e_n in p-basis: sum_mu eps(mu) p_mu/z_mu."""
    return {mu: Rational(sign_mu(mu),1)/z_mu(mu) for mu in partitions(n)}

def add(d1, d2):
    out = dict(d1)
    for k, v in d2.items():
        out[k] = out.get(k, Integer(0)) + v
        if out[k] == 0:
            del out[k]
    return out

def scale(d, c):
    return {k: v*c for k, v in d.items() if v*c != 0}

def mul(d1, d2):
    """Multiply two p-basis elements: concatenate partitions (p_mu * p_nu = p_{mu cup nu})."""
    out = {}
    for mu, a in d1.items():
        for nu, b in d2.items():
            comb = tuple(sorted(mu + nu, reverse=True))
            out[comb] = out.get(comb, Integer(0)) + a*b
    return {k: v for k, v in out.items() if v != 0}

def inner(d1, d2):
    """<f,g> with p_mu orthogonal, <p_mu,p_mu>=z_mu."""
    total = Integer(0)
    for mu, a in d1.items():
        if mu in d2:
            total += a * d2[mu] * z_mu(mu)
    return total

# h1 = p_1
def h1():
    return {(1,): Integer(1)}

def power(d, k):
    if k == 0:
        return {(): Integer(1)}
    res = d
    for _ in range(k-1):
        res = mul(res, d)
    return res
