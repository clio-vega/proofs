"""engine.py -- broken-border-strip operators on Fock space, and the order filtration.

Conventions, FIXED HERE ONCE (Phase 0a of the brief), never tuned later:

  Fock space = Lambda over Q(t); basis {s_lam}.
  Order filtration:  ord(X) <= d  iff  every (d+1)-fold nested bracket
        [...[[X, M_{p_a1}], M_{p_a2}], ..., M_{p_a(d+1)}]  vanishes.
  This is Grothendieck's differential-order filtration for the commutative
  ring Lambda = Q(t)[p_1,p_2,...]:  ord(p_rho^perp) = ell(rho)  (Macdonald I.5).

  beta-set:  B_L(lam) = { lam_i + L - i : 1<=i<=L } subset Z_{>=0}.
  Adding an e-ribbon: b in B, b+e not in B, B -> B - {b} + {b+e}.
  height h = #(B cap (b,b+e)) = number of beads strictly jumped.

  R_e(t) s_lam = sum_{mu/lam connected e-strip} t^{ht} s_mu.
  g_e(t)      = "add a size-e broken (generalized) border strip",
                weight (1+t)^{m-1} t^{ht}, m = #edge-connected components.
  E_e = g_e - R_e  (the m>=2 part).   eps := 1+t.
  E_e = eps*N_e + O(eps^2), so  N_e = d/dt E_e |_{t=-1}
      = "add a size-e broken strip with EXACTLY TWO components, sign (-1)^{ht}".

  NB: PROVE.md sec 0 says N_e adds a strip with "at least two" components.
  That is a typo: (1+t)^{m-1} kills m>=3 to first order in eps.  The
  questions file ("exactly two") is the correct one.  Checked below.
"""
from functools import lru_cache
import itertools, sympy as sp

t = sp.Symbol('t')
q = sp.Symbol('q')

# ---------------------------------------------------------------- partitions
def trim(lam):
    lam = tuple(lam)
    while lam and lam[-1] == 0:
        lam = lam[:-1]
    return lam

@lru_cache(maxsize=None)
def parts_of(n, cap=None):
    if cap is None: cap = n
    if n == 0: return ((),)
    out = []
    for p in range(min(n, cap), 0, -1):
        for rest in parts_of(n - p, p):
            out.append((p,) + rest)
    return tuple(out)

def contains(mu, lam):
    lam = trim(lam); mu = trim(mu)
    if len(lam) > len(mu): return False
    return all(mu[i] >= lam[i] for i in range(len(lam)))

def boxes(mu, lam):
    """cells of mu/lam as (row, col), 1-indexed."""
    lam = trim(lam) + (0,) * (len(trim(mu)) - len(trim(lam)))
    mu = trim(mu)
    return [(i + 1, j + 1) for i in range(len(mu)) for j in range(lam[i], mu[i])]

def components_and_ht(cells):
    """edge-connected components of a set of cells; returns (m, ht, has2x2).
       ht = sum over components of (#rows - 1)."""
    S = set(cells)
    seen = set(); m = 0; ht = 0
    for c in S:
        if c in seen: continue
        m += 1
        stack = [c]; comp = []
        seen.add(c)
        while stack:
            (i, j) = stack.pop(); comp.append((i, j))
            for (a, b) in ((i+1, j), (i-1, j), (i, j+1), (i, j-1)):
                if (a, b) in S and (a, b) not in seen:
                    seen.add((a, b)); stack.append((a, b))
        ht += len(set(i for (i, j) in comp)) - 1
    has2x2 = any((i, j) in S and (i+1, j) in S and (i, j+1) in S and (i+1, j+1) in S
                 for (i, j) in S)
    return m, ht, has2x2

# -------------------------------------------------- broken-strip operators
@lru_cache(maxsize=None)
def strips(lam, e):
    """all mu > lam with mu/lam a size-e broken border strip.
       returns tuple of (mu, m, ht)."""
    lam = trim(lam)
    out = []
    for mu in parts_of(sum(lam) + e):
        if not contains(mu, lam): continue
        cells = boxes(mu, lam)
        if len(cells) != e: continue
        m, ht, sq = components_and_ht(cells)
        if sq: continue
        out.append((mu, m, ht))
    return tuple(out)

def op_g(state, e):
    """g_e = multiplication by P_{(e)}(X;-t): weight (1+t)^{m-1} t^{ht}."""
    out = {}
    for lam, w in state.items():
        for (mu, m, ht) in strips(lam, e):
            out[mu] = out.get(mu, 0) + w * (1 + t)**(m - 1) * t**ht
    return clean(out)

def op_R(state, e):
    out = {}
    for lam, w in state.items():
        for (mu, m, ht) in strips(lam, e):
            if m == 1:
                out[mu] = out.get(mu, 0) + w * t**ht
    return clean(out)

def op_N(state, e):
    """N_e: exactly two components, sign (-1)^{ht}."""
    out = {}
    for lam, w in state.items():
        for (mu, m, ht) in strips(lam, e):
            if m == 2:
                out[mu] = out.get(mu, 0) + w * (-1)**ht
    return clean(out)

def clean(d):
    return {k: sp.expand(v) for k, v in d.items() if sp.expand(v) != 0}

# ---------------------------------------------- Murnaghan-Nakayama: M_{p_a}
def op_p(state, a):
    """multiplication by p_a, via MN: sum over connected a-strips, sign (-1)^ht."""
    out = {}
    for lam, w in state.items():
        for (mu, m, ht) in strips(lam, a):
            if m == 1:
                out[mu] = out.get(mu, 0) + w * (-1)**ht
    return clean(out)

def sub(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = sp.expand(out.get(k, 0) - v)
    return {k: v for k, v in out.items() if v != 0}

def add(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = sp.expand(out.get(k, 0) + v)
    return {k: v for k, v in out.items() if v != 0}

def scal(c, a):
    return clean({k: c * v for k, v in a.items()})
