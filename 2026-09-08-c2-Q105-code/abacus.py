"""
Q105 primary engine.  Maya/abacus picture ONLY.

R_g(u)|M> = sum over legal g-moves (b, b+g) of u^{#(M cap (b,b+g))} |M.(b,g)>.

This file shares no code path with reviews/2026-09-08-selfreview-code/ribbon.py
(which enumerates border strips on Young diagrams) or hbasis.py.  It is used to
compute [R_e(t),R_f(s)] and to test the closed forms.  ribbon.py is the
independent cross-check, not the source.
"""
import sympy as sp
from itertools import product

t, s = sp.symbols('t s')


def parts(n):
    if n == 0:
        yield ()
        return
    def rec(rem, mx):
        if rem == 0:
            yield ()
            return
        for k in range(min(rem, mx), 0, -1):
            for tail in rec(rem - k, k):
                yield (k,) + tail
    yield from rec(n, n)


def to_beads(lam, L):
    """Explicit window: beads lam_j - j for j=1..L (lam zero-padded).
    All integers < -L are implicitly beads; all > lam_1 are implicitly holes."""
    lam = list(lam) + [0] * (L - len(lam))
    return frozenset(lam[j] - (j + 1) for j in range(L))


def to_partition(M, L):
    b = sorted(M, reverse=True)
    assert len(b) == L
    lam = tuple(b[j] + (j + 1) for j in range(L))
    assert all(lam[i] >= lam[i + 1] for i in range(L - 1)), lam
    assert lam[-1] == 0, ("window too small", lam)
    return tuple(x for x in lam if x > 0)


def moves(M, g):
    """(b, weight_exponent) for every legal g-move."""
    out = []
    for b in M:
        if b + g in M:
            continue
        w = sum(1 for x in M if b < x < b + g)
        out.append((b, w))
    return out


def apply_R(vec, g, u, L):
    """vec: dict frozenset-Maya -> coeff."""
    out = {}
    for M, c in vec.items():
        for b, w in moves(M, g):
            Mp = frozenset((M - {b}) | {b + g})
            # guard: the implicit tail must not be reachable
            assert min(Mp) >= -L, "window too small"
            out[Mp] = out.get(Mp, 0) + c * u ** w
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def commutator(lam, e, f, u_e=t, u_f=s):
    """<M'| [R_e(u_e), R_f(u_f)] |lam>, dict Maya -> coeff, plus (M, L)."""
    L = len(lam) + 2 * (e + f) + 4
    M0 = to_beads(lam, L)
    v = {M0: sp.Integer(1)}
    ef = apply_R(apply_R(v, f, u_f, L), e, u_e, L)   # R_e R_f  (R_f first)
    fe = apply_R(apply_R(v, e, u_e, L), f, u_f, L)   # R_f R_e
    out = dict(ef)
    for k, val in fe.items():
        out[k] = out.get(k, 0) - val
    out = {k: sp.expand(v_) for k, v_ in out.items() if sp.expand(v_) != 0}
    return out, M0, L


def stats(M, a, e, f):
    """N, A, B, m_e, m_f in the Q96 Thm 4.1 convention (A uses f, B uses e)."""
    cnt = lambda lo, hi: sum(1 for x in M if lo < x < hi)
    N = cnt(a, a + e + f)
    A = cnt(a, a + f)
    B = cnt(a, a + e)
    return N, A, B, int(a + e in M), int(a + f in M)


def q96_one_bead(M, a, e, f, u_e=t, u_f=s):
    N, A, B, me, mf = stats(M, a, e, f)
    return sp.expand(u_f**A * u_e**(N - 1 - A) * (u_e - mf * (1 + u_e))
                     - u_e**B * u_f**(N - 1 - B) * (u_f - me * (1 + u_f)))


def kappa(M, a, e, f):
    N, A, B, me, mf = stats(M, a, e, f)
    return 2 * (A + B - N) + me + mf, N
