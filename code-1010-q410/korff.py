"""Korff 1110.6356 cylindric Hall-Littlewood weight, in HKKO/my coordinates.

DICTIONARY (see scratch/q410/DICTIONARY.md):
  Korff n_K = h  = number of cyclic slots      (= 2*k_mine in Warnaar/HKKO)
  Korff k_K = w  = level / total particles     (= 2*ell_mine+2)
A state is x = (x_1..x_h) in R = {x_1>=..>=x_h>=x_1-w}; Korff's particle config is
  m_i = x_i - x_{i+1} (1<=i<h),   m_h = w - (x_1 - x_h),   sum_i m_i = w,  m_i >= 0.
A step is theta in {0,1}^h, x -> x+theta (staying in R); then
  m_i(new) - m_i(old) = theta_i - theta_{i+1}  (cyclically, theta_{h+1}:=theta_1),
which is Korff's own relation in the proof of his Lemma 5.4.
Korff (5.13):  Psi = prod_{i in J} (1 - t^{m_i(mu)}),  J = {i in Z_h : theta_i=0, theta_{i+1}=1}.
Macdonald's ordinary psi is the same with the index set NOT wrapped (i = 1..h-1 only).
"""
import itertools
from functools import lru_cache
import sympy as sp

t = sp.Symbol('t')


def step_vectors(h, a):
    out = []
    for S in itertools.combinations(range(h), a):
        v = [0]*h
        for i in S:
            v[i] = 1
        out.append(tuple(v))
    return out


def in_region(x, w):
    h = len(x)
    return all(x[i] >= x[i+1] for i in range(h-1)) and x[h-1] >= x[0] - w


def mvec(x, w):
    """Korff's particle configuration m in Z_h, sum = w."""
    h = len(x)
    m = [x[i]-x[i+1] for i in range(h-1)]
    m.append(w - (x[0]-x[h-1]))
    return tuple(m)


def J_cyclic(theta):
    h = len(theta)
    return [i for i in range(h) if theta[i] == 0 and theta[(i+1) % h] == 1]


def J_open(theta):
    """Macdonald's index set: no wrap (theta_{h+1} = 0)."""
    h = len(theta)
    return [i for i in range(h-1) if theta[i] == 0 and theta[i+1] == 1]


def I_cyclic(theta):
    h = len(theta)
    return [i for i in range(h) if theta[i] == 1 and theta[(i+1) % h] == 0]


def weight_step(x, theta, w, kind="Psi"):
    """Korff (5.13) Psi  (or (5.12) Phi) for one cylindric horizontal strip x -> x+theta."""
    m = mvec(x, w)
    if kind == "Psi":
        idx, mm = J_cyclic(theta), m
    elif kind == "Psi_open":
        idx, mm = J_open(theta), m
    elif kind == "Phi":
        idx, mm = I_cyclic(theta), mvec(tuple(x[i]+theta[i] for i in range(len(x))), w)
    else:
        raise ValueError(kind)
    out = sp.Integer(1)
    for i in idx:
        assert mm[i] >= 1, (x, theta, w, i, mm)      # region condition must force this
        out *= (1 - t**mm[i])
    return sp.expand(out)


def paths(h, w, alpha):
    """All HKKO lattice paths (prop:cssyt=path): 0 -> endpoint, step a has alpha_a ones.
    Returns list of tuples of states (x^0=0, x^1, ..., x^len(alpha))."""
    cur = [((tuple([0]*h),))]
    for a in alpha:
        nxt = []
        for P in cur:
            x = P[-1]
            for v in step_vectors(h, a):
                q = tuple(x[i]+v[i] for i in range(h))
                if in_region(q, w):
                    nxt.append(P + (q,))
        cur = nxt
    return cur


def weight_path(P, w, kind="Psi"):
    out = sp.Integer(1)
    for a in range(len(P)-1):
        theta = tuple(P[a+1][i]-P[a][i] for i in range(len(P[a])))
        out *= weight_step(P[a], theta, w, kind)
    return sp.expand(out)


def cyl_HL(h, w, nvars, endpoint=None, kind="Psi"):
    """sum over cylindric tableaux (all contents alpha in Z_{>=0}^nvars) of
    weight(T) * x^T, restricted to a fixed endpoint if given.  = Korff's P_{lam/d/mu}
    with mu the state x=0, summed over all shapes if endpoint is None."""
    xs = sp.symbols(f'x1:{nvars+1}', positive=True)
    tot = {}
    for alpha in itertools.product(range(h+1), repeat=nvars):
        for P in paths(h, w, alpha):
            if endpoint is not None and P[-1] != tuple(endpoint):
                continue
            wt = weight_path(P, w, kind)
            mon = sp.prod([xs[i]**alpha[i] for i in range(nvars)])
            tot[mon] = tot.get(mon, 0) + wt
    return sp.expand(sum(v*k for k, v in tot.items())), xs
