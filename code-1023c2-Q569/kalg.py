r"""
INDEPENDENT mechanism for c, c~, d, d~ : build the ring directly.

G_lambda = sum_T (-1)^{|T|-|lambda|} x^T  over semistandard SET-VALUED tableaux
(ktiles.tex:711-716, Buch's definition).  Multiply/comultiply as actual
polynomials, then expand in the G-basis by degree-triangular peeling
(G_la = s_la + higher degree).  Nothing here touches a puzzle or any tableau
rule from the paper -- it goes through the ring.
"""
import itertools
from functools import lru_cache

def cells(shape):
    return [(i, j) for i, r in enumerate(shape) for j in range(r)]

def svt(shape, N, maxextra=None):
    """Yield set-valued tableaux as dict cell->frozenset, entries in 1..N.
    rows: max T(i,j) <= min T(i,j+1);  cols: max T(i,j) < min T(i+1,j)."""
    cs = cells(shape)
    subsets = []
    for r in range(1, N + 1):
        for c in itertools.combinations(range(1, N + 1), r):
            subsets.append(frozenset(c))
    sz = len(cells(shape))
    def rec(idx, T, extra):
        if maxextra is not None and extra > maxextra: return
        if idx == len(cs):
            yield dict(T); return
        i, j = cs[idx]
        for S in subsets:
            if j > 0 and not (max(T[(i, j-1)]) <= min(S)): continue
            if i > 0 and not (max(T[(i-1, j)]) < min(S)): continue
            T[(i, j)] = S
            yield from rec(idx + 1, T, extra + len(S) - 1)
            del T[(i, j)]
    yield from rec(0, {}, 0)

def Gpoly(shape, N, maxdeg=None):
    """G_shape(x_1..x_N) as dict: exponent tuple (len N) -> coeff."""
    shape = tuple(x for x in shape if x > 0)
    n0 = sum(shape)
    me = None if maxdeg is None else maxdeg - n0
    out = {}
    for T in svt(shape, N, maxextra=me):
        e = [0] * N
        for S in T.values():
            for v in S: e[v - 1] += 1
        tot = sum(e)
        if maxdeg is not None and tot > maxdeg: continue
        key = tuple(e)
        out[key] = out.get(key, 0) + (-1) ** (tot - n0)
    return {k: v for k, v in out.items() if v}

def Gpoly_split(shape, N, M, maxdeg=None):
    """G_shape(x_1..x_N, y_1..y_M) as dict: (xexp, yexp) -> coeff."""
    shape = tuple(x for x in shape if x > 0)
    n0 = sum(shape)
    me = None if maxdeg is None else maxdeg - n0
    out = {}
    for T in svt(shape, N + M, maxextra=me):
        ex = [0] * N; ey = [0] * M
        for S in T.values():
            for v in S:
                if v <= N: ex[v - 1] += 1
                else: ey[v - N - 1] += 1
        tot = sum(ex) + sum(ey)
        if maxdeg is not None and tot > maxdeg: continue
        key = (tuple(ex), tuple(ey))
        out[key] = out.get(key, 0) + (-1) ** (tot - n0)
    return {k: v for k, v in out.items() if v}

def mul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            k = tuple(x + y for x, y in zip(a, b))
            out[k] = out.get(k, 0) + ca * cb
    return {k: v for k, v in out.items() if v}

@lru_cache(maxsize=None)
def Spoly(shape, N):
    """Schur polynomial via SSYT."""
    shape = tuple(x for x in shape if x > 0)
    if len(shape) > N: return {}
    cs = cells(shape)
    out = {}
    def rec(idx, T):
        if idx == len(cs):
            e = [0] * N
            for v in T.values(): e[v - 1] += 1
            k = tuple(e); out[k] = out.get(k, 0) + 1
            return
        i, j = cs[idx]
        lo = 1
        if j > 0: lo = max(lo, T[(i, j-1)])
        if i > 0: lo = max(lo, T[(i-1, j)] + 1)
        for v in range(lo, N + 1):
            T[(i, j)] = v; rec(idx + 1, T); del T[(i, j)]
    rec(0, {})
    return out

def schur_expand(p, N):
    """Expand symmetric polynomial p in Schur basis by leading-monomial peeling."""
    p = dict(p); res = {}
    while p:
        lead = max(p.keys())                 # lex-largest exponent
        assert all(lead[i] >= lead[i+1] for i in range(N-1)), ("not symmetric?", lead)
        la = tuple(x for x in lead if x > 0)
        c = p[lead]; res[la] = res.get(la, 0) + c
        s = Spoly(la, N)
        for k, v in s.items():
            p[k] = p.get(k, 0) - c * v
            if p[k] == 0: del p[k]
    return {k: v for k, v in res.items() if v}

def G_expand(p, N, maxdeg):
    """Expand p in the G-basis.  G_la = s_la + (strictly higher degree), so
    peel in increasing order of |la|."""
    sc = schur_expand(p, N)
    res = {}
    for d in range(0, maxdeg + 1):
        for la in sorted([l for l in list(sc) if sum(l) == d], reverse=True):
            c = sc.get(la, 0)
            if c == 0: continue
            res[la] = res.get(la, 0) + c
            g = schur_expand(Gpoly(la, N, maxdeg=maxdeg), N)
            for k, v in g.items():
                sc[k] = sc.get(k, 0) - c * v
                if sc[k] == 0: del sc[k]
    leftover = {k: v for k, v in sc.items() if v and sum(k) <= maxdeg}
    assert not leftover, ("G-expansion residue", leftover)
    return {k: v for k, v in res.items() if v}

# ----------------------------------------------------- structure constants
def c_coeffs(la, mu, N, maxdeg):
    """c^nu_{la mu}:  G_la G_mu = sum (-1)^{|nu|-|la|-|mu|} c G_nu."""
    p = mul(Gpoly(la, N, maxdeg=maxdeg), Gpoly(mu, N, maxdeg=maxdeg))
    p = {k: v for k, v in p.items() if sum(k) <= maxdeg}
    A = G_expand(p, N, maxdeg)
    n0 = sum(la) + sum(mu)
    return {nu: (-1) ** (sum(nu) - n0) * v for nu, v in A.items()}

def ct_coeffs(la, mu, N, maxdeg):
    """c~ : Gam_la = G_la (1 - G_1);  Gam_la Gam_mu = sum (-1)^... c~ Gam_nu.
    Equivalently expand (G_la G_mu)(1-G_1) in the G-basis after dividing...
    we instead expand Gam_la*Gam_mu directly in the Gam-basis, which is
    degree-triangular for the same reason (Gam_la = s_la + higher)."""
    G1 = Gpoly([1], N, maxdeg=maxdeg)
    one = {tuple([0]*N): 1}
    fac = {k: v for k, v in ({**one, **{k: -v for k, v in G1.items()}}).items() if v}
    def Gam(sh):
        r = mul(Gpoly(sh, N, maxdeg=maxdeg), fac)
        return {k: v for k, v in r.items() if sum(k) <= maxdeg}
    p = mul(Gam(la), Gam(mu))
    p = {k: v for k, v in p.items() if sum(k) <= maxdeg}
    # peel in Gam basis
    sc = schur_expand(p, N); res = {}
    for d in range(0, maxdeg + 1):
        for nu in sorted([l for l in list(sc) if sum(l) == d], reverse=True):
            cc = sc.get(nu, 0)
            if cc == 0: continue
            res[nu] = res.get(nu, 0) + cc
            g = schur_expand(Gam(nu), N)
            for k, v in g.items():
                sc[k] = sc.get(k, 0) - cc * v
                if sc[k] == 0: del sc[k]
    n0 = sum(la) + sum(mu)
    return {nu: (-1) ** (sum(nu) - n0) * v for nu, v in res.items() if v}

def d_coeffs(nu, N, M, maxdeg):
    """d^nu_{la mu}: Delta(G_nu) = sum (-1)^{|nu|-|la|-|mu|} d G_la (x) G_mu."""
    p = Gpoly_split(nu, N, M, maxdeg=maxdeg)
    # double Schur expansion
    byy = {}
    for (ex, ey), c in p.items():
        byy.setdefault(ey, {})[ex] = byy.setdefault(ey, {}).get(ex, 0) + c
    # expand x-side in Schur, collect
    sxy = {}
    for ey, px in byy.items():
        for la, c in schur_expand(px, N).items():
            sxy.setdefault(la, {})[ey] = sxy.setdefault(la, {}).get(ey, 0) + c
    dbl = {}
    for la, py in sxy.items():
        for mu, c in schur_expand(py, M).items():
            dbl[(la, mu)] = dbl.get((la, mu), 0) + c
    # peel in G (x) G basis by total degree
    res = {}
    for tot in range(0, maxdeg + 1):
        for (la, mu) in sorted([k for k in list(dbl) if sum(k[0]) + sum(k[1]) == tot], reverse=True):
            c = dbl.get((la, mu), 0)
            if c == 0: continue
            res[(la, mu)] = res.get((la, mu), 0) + c
            gx = schur_expand(Gpoly(la, N, maxdeg=maxdeg), N)
            gy = schur_expand(Gpoly(mu, M, maxdeg=maxdeg), M)
            for a, va in gx.items():
                for b, vb in gy.items():
                    if sum(a) + sum(b) > maxdeg: continue
                    dbl[(a, b)] = dbl.get((a, b), 0) - c * va * vb
                    if dbl[(a, b)] == 0: del dbl[(a, b)]
    n0 = sum(nu)
    return {k: (-1) ** (n0 - sum(k[0]) - sum(k[1])) * v for k, v in res.items() if v}
