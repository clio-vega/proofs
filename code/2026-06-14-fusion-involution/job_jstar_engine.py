"""
job_jstar_engine.py  --  shared engine for the |J*|-even / SYT-chain study.

M_j(lambda) = <s_lambda, h_1^{2(m-j)} e_2^j>,  lambda |- 2m.

Two independent computations (cross-checked in job_jstar_crosscheck.py):
  (A) build-up Pieri:  e_2^j then h_1^{2(m-j)} from the empty partition, read coeff of lambda.
  (B) SYT-chain model: M_j = sum_mu N^{(j)}_{lambda->mu} f^mu, where N counts length-j chains
      lambda -> mu each step REMOVING a vertical 2-strip (mu |- 2(m-j)).  [adjoint of (A)]

(B) is the load-bearing model for PROVE and is what we scale to m<=16.
"""
from sympy import binomial
from collections import defaultdict
from functools import lru_cache
from math import factorial

# ---------- 2-adic / binomial helpers ----------
def v2(n):
    if n == 0:
        return None
    c = 0
    while n % 2 == 0:
        n //= 2; c += 1
    return c

def C(n, k):
    if k < 0 or k > n or n < 0:
        return 0
    return int(binomial(n, k))

def carries(x, y):
    """v2(C(x+y, x)) = number of carries adding x and y in base 2 (Kummer)."""
    c = 0; cy = 0
    while x > 0 or y > 0 or cy > 0:
        s = (x & 1) + (y & 1) + cy; cy = s >> 1; c += cy; x >>= 1; y >>= 1
    return c

def vC(n, k):
    if k < 0 or k > n:
        return None
    return carries(k, n - k)

# ---------- partitions ----------
@lru_cache(maxsize=None)
def hook_f(lam):
    """Number of SYT of shape lam (hook length formula)."""
    lam = [x for x in lam if x > 0]
    if not lam:
        return 1
    n = sum(lam)
    conj = []
    m = lam[0]
    for c in range(m):
        conj.append(sum(1 for r in lam if r > c))
    prod = 1
    for i, row in enumerate(lam):
        for jcol in range(row):
            arm = row - jcol - 1
            leg = conj[jcol] - i - 1
            prod *= (arm + leg + 1)
    return factorial(n) // prod

def add_single(mu):
    res = []; mu = list(mu); L = len(mu)
    for i in range(L):
        if i == 0 or mu[i-1] > mu[i]:
            nu = mu[:]; nu[i] += 1; res.append(tuple(nu))
    res.append(tuple(mu + [1]))
    return res

def add_vstrip2(mu):
    """All nu with nu/mu a vertical 2-strip (two distinct rows each +1)."""
    mu = list(mu); ext = mu + [0, 0]; n = len(ext); res = []
    for i in range(n):
        for k in range(i+1, n):
            nu = ext[:]; nu[i] += 1; nu[k] += 1
            if all(nu[t] >= nu[t+1] for t in range(n-1)):
                while nu and nu[-1] == 0:
                    nu.pop()
                res.append(tuple(nu))
    return res

def remove_vstrip2(lam):
    """All mu with lam/mu a vertical 2-strip (two distinct rows each -1, mu a partition)."""
    lam = list(lam); n = len(lam); res = []
    for i in range(n):
        for k in range(i+1, n):
            nu = lam[:]; nu[i] -= 1; nu[k] -= 1
            if nu[i] < 0 or nu[k] < 0:
                continue
            if all(nu[t] >= nu[t+1] for t in range(n-1)) and nu[-1] >= 0:
                while nu and nu[-1] == 0:
                    nu.pop()
                res.append(tuple(nu))
    return res

def mult(state, op):
    new = defaultdict(int)
    for mu, c in state.items():
        for nu in op(mu):
            new[nu] += c
    return new

# ---------- M_j ----------
def Mj_buildup(lam, j, m):
    """Method (A): empty -> e_2^j -> h_1^{2(m-j)}, read coeff of lam.  (small m only)"""
    state = {(): 1}
    for _ in range(j):
        state = mult(state, add_vstrip2)
    for _ in range(2*(m-j)):
        state = mult(state, add_single)
    return state.get(tuple(lam), 0)

def Mj_chain_all(lam):
    """Method (B): return [M_0, M_1, ..., M_m] for lam |- 2m via vertical-2-strip removal chains."""
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    assert n % 2 == 0
    m = n // 2
    Ms = []
    state = {lam: 1}
    for level in range(m + 1):
        Ms.append(sum(c * hook_f(mu) for mu, c in state.items()))
        if level < m:
            state = mult(state, remove_vstrip2)
    return Ms

# ---------- partitions of n ----------
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest

# ---------- val / J* analysis ----------
def analyze(lam):
    """Return dict with M, val, Jstar, |Jstar|, V(=min val), parity of Jstar."""
    lam = tuple(x for x in lam if x > 0)
    m = sum(lam) // 2
    Ms = Mj_chain_all(lam)
    vals = {}
    for j in range(m + 1):
        Mj = Ms[j]; cj = C(m, j)
        if Mj == 0 or cj == 0:
            continue
        vals[j] = j + 2 * (v2(cj) + v2(Mj))
    Vmin = min(vals.values())
    Jstar = sorted(j for j, vv in vals.items() if vv == Vmin)
    parities = set(j % 2 for j in Jstar)
    return {
        "lam": lam, "m": m, "M": Ms, "val": vals,
        "Jstar": Jstar, "absJstar": len(Jstar), "V": Vmin,
        "Jparity": (parities.pop() if len(parities) == 1 else "MIXED"),
    }
