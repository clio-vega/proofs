"""Region partial sums G_R for the m=2 cylindric Kostka slice sum.

Everything here is a direct transcription of the PROVED reduction in
proofs/2026-09-30-c1-cylindric-kostka-logconcavity.tex (prop:reform, prop:m2red,
the four-region table).  Free parameters: n, A, B, u, T-, T+  with C=B+1, D=A+n-1.
"""

def conv_interval(a, b, c, d):
    """1_[a,b] * 1_[c,d] as dict s -> count. Empty if a>b or c>d."""
    out = {}
    if a > b or c > d:
        return out
    for i in range(a, b + 1):
        for j in range(c, d + 1):
            out[i + j] = out.get(i + j, 0) + 1
    return out

def region_of(t, tau1, tau2):
    al = 1 if t >= tau1 else 0
    be = 1 if t > tau2 else 0
    return {(0,0):'I', (1,0):'II', (1,1):'III', (0,1):'IV'}[(al,be)]

def endpoints(t, n, A, B, u):
    C, D = B + 1, A + n - 1
    at = max(t, A); bt = min(u - 1 - t, B)
    ct = max(u - t, C); dt = min(t + n - 1, D)
    return at, bt, ct, dt

def region_sums(n, A, B, u, Tm, Tp):
    """Returns dict region -> (dict s->value), plus the t-sets."""
    tau1 = A
    tau2 = u - 1 - B
    G = {'I': {}, 'II': {}, 'III': {}, 'IV': {}}
    tset = {'I': [], 'II': [], 'III': [], 'IV': []}
    for t in range(Tm, Tp + 1):
        at, bt, ct, dt = endpoints(t, n, A, B, u)
        R = region_of(t, tau1, tau2)
        tset[R].append(t)
        for s, v in conv_interval(at, bt, ct, dt).items():
            G[R][s] = G[R].get(s, 0) + v
    return G, tset

def as_seq(g):
    """dict -> (lo, list) with zeros filled; ([],0) if empty."""
    if not g:
        return 0, []
    lo, hi = min(g), max(g)
    return lo, [g.get(s, 0) for s in range(lo, hi + 1)]

def get(g, s):
    return g.get(s, 0)

def is_pf2(g):
    """nonneg, interval support, g(s)^2 >= g(s-1)g(s+1)."""
    if not g:
        return True
    lo, hi = min(g), max(g)
    seq = [g.get(s, 0) for s in range(lo, hi + 1)]
    if any(v < 0 for v in seq):
        return False
    # interval support
    nz = [k for k, v in enumerate(seq) if v > 0]
    if nz and nz[-1] - nz[0] + 1 != len(nz):
        return False
    for s in range(lo - 1, hi + 2):
        if get(g, s) ** 2 < get(g, s - 1) * get(g, s + 1):
            return False
    return True

def star_fails(g1, g2):
    """List of s where 2 g1(s)g2(s) < g1(s-1)g2(s+1) + g1(s+1)g2(s-1)."""
    if not g1 or not g2:
        return []
    lo = min(min(g1), min(g2)) - 2
    hi = max(max(g1), max(g2)) + 2
    bad = []
    for s in range(lo, hi + 1):
        lhs = 2 * get(g1, s) * get(g2, s)
        rhs = get(g1, s - 1) * get(g2, s + 1) + get(g1, s + 1) * get(g2, s - 1)
        if lhs < rhs:
            bad.append((s, lhs, rhs))
    return bad

# ---- cylindric parameters ----
def cyl_params(n, mu, lam, b):
    mu1, mu2 = mu; lam1, lam2 = lam
    A = lam2 - n + 1; B = lam1
    u = mu1 + mu2 + b
    Tm = max(mu1, u - mu1 - n + 1)
    Tp = min(mu2 - 1, u - mu2)
    return n, A, B, u, Tm, Tp

def valid_shape(n, mu, lam):
    mu1, mu2 = mu; lam1, lam2 = lam
    # cylindric shape lambda/mu at m=2, level n: mu_1<=mu_2<=mu_1+n, mu<=lam componentwise,
    # and lam must itself be a cylindric shape: lam1<=lam2<=lam1+n
    return (mu1 <= mu2 <= mu1 + n and lam1 <= lam2 <= lam1 + n
            and mu1 <= lam1 and mu2 <= lam2)
