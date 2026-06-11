"""
Core machinery for d=4 non-vanishing investigation.

Object: G_lambda(i) = <s_lambda, psi^m>,  psi = h_2 + i e_2 = s_2 + i s_{11}, n=2m.

We work with Gaussian integers exactly: represent a Gaussian integer as a pair (a,b) meaning a + b*i.
"""
from functools import lru_cache
from itertools import product as iproduct

# ---------- partitions ----------

def partitions(n):
    """All partitions of n as tuples (weakly decreasing)."""
    def gen(n, maxpart):
        if n == 0:
            yield ()
            return
        for k in range(min(n, maxpart), 0, -1):
            for rest in gen(n - k, k):
                yield (k,) + rest
    return list(gen(n, n))

def conj(lam):
    """Conjugate (transpose) partition."""
    if not lam:
        return ()
    m = lam[0]
    return tuple(sum(1 for p in lam if p >= j) for j in range(1, m + 1))

def is_partition(lam):
    return all(lam[k] >= lam[k+1] for k in range(len(lam)-1)) and all(p > 0 for p in lam)

# ---------- Pieri: add a horizontal / vertical r-strip ----------
# s_r * s_mu = sum over mu subset nu, nu/mu horizontal r-strip
# s_{1^r} * s_mu = sum over nu/mu vertical r-strip

def pad(lam, length):
    return tuple(lam) + (0,) * (length - len(lam))

def horizontal_strip_adds(mu, r):
    """Return list of partitions nu = mu + horizontal r-strip.
    Horizontal strip: no two added boxes in same column.
    nu/mu horiz strip <=> mu_1>=nu_2>=mu_2>=nu_3>=... (interlacing) and |nu|-|mu|=r.
    """
    L = len(mu) + 1
    mup = pad(mu, L)  # mu_1..mu_L (with mu_{L}=0)
    results = []
    # nu has length <= L. nu_j for j=1..L. Constraint: nu_j integer,
    # nu_1>=nu_2>=...>=0 ; horizontal strip: nu_j >= mu_j and nu_{j} <= mu_{j-1} (mu_0=inf).
    # i.e. mu_j <= nu_j <= mu_{j-1}, and sum(nu_j-mu_j)=r.
    # enumerate
    def rec(j, remaining, cur):
        if j == L:
            if remaining == 0:
                nu = tuple(x for x in cur if x > 0)
                results.append(nu)
            return
        mu_jm1 = mup[j-1] if j-1 >= 0 else 10**9
        if j == 0:
            upper = 10**9
        else:
            upper = mup[j-1]
        mu_j = mup[j]
        # nu_j in [mu_j, upper], add = nu_j-mu_j in [0, remaining]
        lo = mu_j
        hi = min(upper, mu_j + remaining)
        for nu_j in range(lo, hi+1):
            rec(j+1, remaining - (nu_j - mu_j), cur + [nu_j])
    rec(0, r, [])
    # dedupe + validate partitions
    out = []
    seen = set()
    for nu in results:
        if is_partition(nu) and nu not in seen:
            seen.add(nu)
            out.append(nu)
    return out

def vertical_strip_adds(mu, r):
    """nu = mu + vertical r-strip = conjugate of horizontal."""
    muc = conj(mu)
    return [conj(nu) for nu in horizontal_strip_adds(muc, r)]

# ---------- Schur algebra over Gaussian integers ----------
# Represent symmetric function as dict: partition tuple -> (a,b) for a+b*i.

def gadd(x, y):
    return (x[0]+y[0], x[1]+y[1])

def gmul(x, y):
    a,b = x; c,d = y
    return (a*c - b*d, a*d + b*c)

def sf_add(f, g):
    out = dict(f)
    for k,v in g.items():
        out[k] = gadd(out.get(k, (0,0)), v)
    return {k:v for k,v in out.items() if v != (0,0)}

def sf_scale(f, c):
    return {k: gmul(v, c) for k,v in f.items()}

def mult_by_psi(f):
    """Multiply symmetric function f (Schur-basis dict) by psi = s_2 + i*s_{11}.
    Returns new dict."""
    out = {}
    for mu, coeff in f.items():
        # s_2 * s_mu : horizontal 2-strips, weight 1
        for nu in horizontal_strip_adds(mu, 2):
            out[nu] = gadd(out.get(nu,(0,0)), gmul(coeff, (1,0)))
        # i*s_{11} * s_mu : vertical 2-strips, weight i
        for nu in vertical_strip_adds(mu, 2):
            out[nu] = gadd(out.get(nu,(0,0)), gmul(coeff, (0,1)))
    return {k:v for k,v in out.items() if v != (0,0)}

def psi_power_schur(m):
    """psi^m in Schur basis. Returns dict partition->(a,b)."""
    f = {(): (1,0)}
    for _ in range(m):
        f = mult_by_psi(f)
    return f

def G_direct(lam):
    """G_lambda(i) via symbolic psi^m expansion. Returns (a,b)."""
    n = sum(lam)
    assert n % 2 == 0
    m = n // 2
    f = psi_power_schur(m)
    return f.get(tuple(lam), (0,0))

# ---------- Chain enumeration (combinatorial model) ----------

def enumerate_chains(lam):
    """Enumerate all 2-strip chains empty -> lam.
    Each step is ('h', nu) or ('v', nu).
    Returns list of chains; each chain is list of ('h'/'v') step-types in order.
    But careful: a step could be BOTH h and v (disconnected pair diff row diff col).
    We treat h-step and v-step as DISTINCT objects per Pieri coefficient.
    Return: list of v-counts (one per chain), and optionally full chains.
    """
    n = sum(lam)
    m = n // 2
    lam = tuple(lam)
    # BFS/DFS forward from empty. We need chains ending at lam.
    # Represent partial chain by current partition; recurse.
    chains = []  # list of tuples of step letters
    def rec(cur, steps):
        if sum(cur) == n:
            if cur == lam:
                chains.append(tuple(steps))
            return
        if sum(cur) > n:
            return
        # h-steps
        for nu in horizontal_strip_adds(cur, 2):
            # prune: nu must be subset of lam (containment) for chain to reach lam
            if subset(nu, lam):
                rec(nu, steps + ['h'])
        for nu in vertical_strip_adds(cur, 2):
            if subset(nu, lam):
                rec(nu, steps + ['v'])
    rec((), [])
    return chains

def enumerate_chains_full(lam):
    """Like enumerate_chains but returns list of (step_letters, partition_sequence)."""
    n = sum(lam)
    lam = tuple(lam)
    chains = []
    def rec(cur, steps, seq):
        if sum(cur) == n:
            if cur == lam:
                chains.append((tuple(steps), tuple(seq)))
            return
        for nu in horizontal_strip_adds(cur, 2):
            if subset(nu, lam):
                rec(nu, steps + ['h'], seq + [nu])
        for nu in vertical_strip_adds(cur, 2):
            if subset(nu, lam):
                rec(nu, steps + ['v'], seq + [nu])
    rec((), [], [()])
    return chains

def subset(mu, lam):
    """Is mu contained in lam (Young diagram containment)?"""
    if len(mu) > len(lam):
        return False
    return all(mu[k] <= lam[k] for k in range(len(mu)))

def G_chains(lam):
    """G via chain model. Returns (a,b) and the v-distribution dict."""
    chains = enumerate_chains(lam)
    vdist = {}
    for steps in chains:
        v = sum(1 for s in steps if s == 'v')
        vdist[v] = vdist.get(v, 0) + 1
    # G = sum i^v
    a = b = 0
    for v, c in vdist.items():
        r = v % 4
        if r == 0: a += c
        elif r == 1: b += c
        elif r == 2: a -= c
        else: b -= c
    return (a, b), vdist

def residues(vdist):
    n = [0,0,0,0]
    for v,c in vdist.items():
        n[v%4] += c
    return tuple(n)
