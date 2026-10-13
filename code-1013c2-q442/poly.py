"""Minimal (Laurent-)polynomial arithmetic over Z, and over Z[t], in n variables.

Monomials are exponent tuples of length n (integers, possibly negative for Laurent).
Coefficients are ints, or dicts {t-degree: int} for the Z[t] layer.

Written for Q442 (PROVE 2026-10-13 c2).  No external dependencies beyond the
standard library, so this engine shares no code path with SymPy.
"""
from itertools import permutations


# ---------------------------------------------------------------- Z coefficients
def zero(): return {}

def const(c, n): return {tuple([0]*n): c} if c else {}

def add(p, q):
    r = dict(p)
    for m, c in q.items():
        v = r.get(m, 0) + c
        if v: r[m] = v
        else: r.pop(m, None)
    return r

def neg(p): return {m: -c for m, c in p.items()}

def sub(p, q): return add(p, neg(q))

def mul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(a + b for a, b in zip(m1, m2))
            v = r.get(m, 0) + c1 * c2
            if v: r[m] = v
            else: r.pop(m, None)
    return r

def scal(c, p):
    if c == 0: return {}
    return {m: c * v for m, v in p.items()}


# ---------------------------------------------------------------- symmetric fns
def elem(n):
    """e_0..e_n in n variables, as polynomials."""
    out = [const(1, n)]
    for i in range(1, n + 1):
        out.append({})
    # e_r = sum over r-subsets
    from itertools import combinations
    for r in range(1, n + 1):
        d = {}
        for S in combinations(range(n), r):
            m = [0] * n
            for i in S: m[i] = 1
            d[tuple(m)] = 1
        out[r] = d
    return out

def e_of(E, i, n):
    """e_i with the convention e_i = 0 for i<0 or i>n."""
    if i < 0 or i > n: return {}
    return E[i]

def f_r(E, r, n):
    """f_r = sum_{i in Z} e_i e_{i+r}.  f_r = f_{-r}."""
    r = abs(r)
    acc = {}
    for i in range(0, n + 1):
        if i + r > n: break
        acc = add(acc, mul(E[i], E[i + r]))
    return acc


# ---------------------------------------------------------------- determinants
def det(M, n):
    """Determinant of a k x k matrix of polynomials, by permutation expansion."""
    k = len(M)
    acc = {}
    for perm in permutations(range(k)):
        sgn = 1
        # count inversions
        inv = sum(1 for a in range(k) for b in range(a + 1, k) if perm[a] > perm[b])
        if inv % 2: sgn = -1
        term = const(1, n)
        for i in range(k):
            term = mul(term, M[i][perm[i]])
            if not term: break
        if term:
            acc = add(acc, scal(sgn, term))
    return acc


# ---------------------------------------------------------------- Z[t] layer
def tzero(): return {}

def tconst(c): return {0: c} if c else {}

def tadd(a, b):
    r = dict(a)
    for d, c in b.items():
        v = r.get(d, 0) + c
        if v: r[d] = v
        else: r.pop(d, None)
    return r

def tneg(a): return {d: -c for d, c in a.items()}

def tmul(a, b):
    r = {}
    for d1, c1 in a.items():
        for d2, c2 in b.items():
            d = d1 + d2
            v = r.get(d, 0) + c1 * c2
            if v: r[d] = v
            else: r.pop(d, None)
    return r

def tshift(a, k):
    """multiply by t^k"""
    return {d + k: c for d, c in a.items()}


# Polynomials over Z[t]: dict monomial -> dict(tdeg -> int)
def Tadd(p, q):
    r = {m: dict(c) for m, c in p.items()}
    for m, c in q.items():
        v = tadd(r.get(m, {}), c)
        if v: r[m] = v
        else: r.pop(m, None)
    return r

def Tneg(p): return {m: tneg(c) for m, c in p.items()}

def Tsub(p, q): return Tadd(p, Tneg(q))

def Tmul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(a + b for a, b in zip(m1, m2))
            v = tadd(r.get(m, {}), tmul(c1, c2))
            if v: r[m] = v
            else: r.pop(m, None)
    return r

def Tlift(p):
    """lift a Z-polynomial to Z[t]-coefficients"""
    return {m: {0: c} for m, c in p.items()}

def Tscal_t(p, k):
    """multiply a Z[t]-polynomial by t^k"""
    return {m: tshift(c, k) for m, c in p.items()}

def Tdet(M, n):
    k = len(M)
    acc = {}
    for perm in permutations(range(k)):
        inv = sum(1 for a in range(k) for b in range(a + 1, k) if perm[a] > perm[b])
        sgn = -1 if inv % 2 else 1
        term = {tuple([0] * n): {0: 1}}
        for i in range(k):
            term = Tmul(term, M[i][perm[i]])
            if not term: break
        if term:
            acc = Tadd(acc, {m: {d: sgn * c for d, c in cc.items()} for m, cc in term.items()})
    return acc

def Tcoeff_t(p, k):
    """extract the t^k coefficient: a Z-polynomial"""
    out = {}
    for m, cc in p.items():
        c = cc.get(k, 0)
        if c: out[m] = c
    return out
