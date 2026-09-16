"""bead.py -- the FERMIONIC side of Conjecture N.  Built from scratch.

NOTHING is imported from the shape engine (probes/2026-09-06-Q84/engine.py).
The only shared ancestor is the definition of a partition as a tuple.

Conventions, fixed here once.

  Maya set   M(lam) = { lam_j - j : j >= 1 } subset Z,  cofinite in Z_{<0}.
             (maya.py, Q59 Phase 1: same normalisation.)

  Wedge      |M> = v_{a_1} ^ v_{a_2} ^ ...   with a_1 > a_2 > ... the elements
             of M in DECREASING order.   |lam> = |M(lam)>.

  Fermions   {psi_a, psi_b^*} = delta_{ab},  {psi_a,psi_b} = {psi_a^*,psi_b^*} = 0.
             psi_a |M> = v_a ^ |M>  (0 if a in M);   psi_a^* = its adjoint.

  Occupation n_j = psi_j psi_j^*,  so  n_j |M> = [j in M] |M>.

CONJECTURE N.   R_e(t) = sum_{b in Z} psi_{b+e} psi_b^* (-t)^{N_(b,b+e)},
                N_(b,b+e) = sum_{b<j<b+e} n_j.

We implement the right-hand side DIRECTLY from the fermionic definition:
the sign of psi_{b+e}psi_b^* is computed by literally sorting the wedge,
not by quoting the formula (-1)^ht.  That keeps the check honest: the
identity sign = (-1)^{#(M cap (b,b+e))} is a *lemma* of the proof, so the
verification must not assume it.
"""
import sympy as sp

t = sp.Symbol('t')


def trim(lam):
    lam = tuple(lam)
    while lam and lam[-1] == 0:
        lam = lam[:-1]
    return lam


def maya(lam, lo):
    """M(lam) intersected with [lo, oo), as a frozenset.  lo must be < -len(lam)."""
    lam = trim(lam)
    n = len(lam)
    assert lo < -n, (lam, lo)
    M = set(lam[j - 1] - j for j in range(1, n + 1))
    for j in range(n + 1, -lo + 1):
        if -j >= lo:
            M.add(-j)
    return frozenset(M)


def from_maya(M, lo):
    """inverse of maya: read lam off a Maya set that agrees with Z_{<0} below lo."""
    a = sorted(M, reverse=True)
    lam = [a[j - 1] + j for j in range(1, len(a) + 1)]
    # the tail must be exactly ...,lo+1,lo  -> contributes zeros
    return trim(tuple(lam))


# ---------------------------------------------------------------- fermions
def psi_star(a, M, sign):
    """psi_a^* on sign*|M>.  Returns (sign', M') or None."""
    if a not in M:
        return None
    # a is the k-th element in decreasing order (1-indexed); moving v_a to the
    # front of the wedge costs (-1)^(k-1).
    k = sum(1 for x in M if x > a)          # = k-1
    return ((-1) ** k * sign, M - {a})


def psi(a, M, sign):
    """psi_a on sign*|M>.  Returns (sign', M') or None."""
    if a in M:
        return None
    # v_a is prepended, then sorted into place past the l-1 larger elements.
    l = sum(1 for x in M if x > a)          # = l-1
    return ((-1) ** l * sign, M | {a})


def bilinear(bplus, b, M):
    """psi_{bplus} psi_b^* |M>.  Returns (sign, M') or None."""
    r = psi_star(b, M, 1)
    if r is None:
        return None
    s, M1 = r
    r = psi(bplus, M1, s)
    if r is None:
        return None
    return r


# ------------------------------------------------- the conjectured operator
def R_bead(lam, e, lo=None, dressing=None, interval='open', shift=0, tsign=-1):
    """RHS of Conjecture N applied to s_lam.  Returns {mu: coeff in Z[t]}.

    The keyword arguments exist ONLY to run the planted-defect controls:
      interval = 'open'   -> j in (b, b+e)          [the conjecture]
                 'closed' -> j in [b, b+e]          [control 1]
      shift    = 0 or 1   -> exponent N + shift     [control 2]
      tsign    = -1 or +1 -> (tsign*t)^N            [control 3: +1 gives t^N]
    """
    lam = trim(lam)
    if lo is None:
        lo = -(len(lam) + 2 * e + 6)
    M = maya(lam, lo)
    hi = (lam[0] if lam else 0) + 2 * e + 6
    out = {}
    for b in range(lo, hi + 1):
        if b not in M:
            continue
        if b + e in M:
            continue
        res = bilinear(b + e, b, M)
        if res is None:
            continue
        sign, Mp = res
        if interval == 'open':
            N = sum(1 for j in range(b + 1, b + e) if j in M)
        elif interval == 'closed':
            N = sum(1 for j in range(b, b + e + 1) if j in M)
        else:
            raise ValueError(interval)
        N += shift
        mu = from_maya(Mp, lo)
        out[mu] = sp.expand(out.get(mu, 0) + sign * (tsign * t) ** N)
    return {k: v for k, v in out.items() if sp.expand(v) != 0}


def open_interval_load(lam, e):
    """max over valid bead moves of #(M cap (b,b+e)) -- the witness statistic."""
    lam = trim(lam)
    lo = -(len(lam) + 2 * e + 6)
    M = maya(lam, lo)
    hi = (lam[0] if lam else 0) + 2 * e + 6
    best = -1
    arg = None
    for b in range(lo, hi + 1):
        if b in M and b + e not in M:
            N = sum(1 for j in range(b + 1, b + e) if j in M)
            if N > best:
                best, arg = N, b
    return best, arg
