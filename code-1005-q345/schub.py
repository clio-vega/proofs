"""Independent ground truth for flag-variety structure constants c_{u,v}^w.

Nothing in this module knows about chains, Monk's rule, or Chevalley.
Schubert polynomials are built from the classical definition:
    S_{w0} = x_1^{n-1} x_2^{n-2} ... x_{n-1},     S_{w s_j} = d_j S_w  if l(w s_j) = l(w)-1
and structure constants are extracted by divided differences.
"""
from itertools import permutations
from functools import lru_cache

# ---------- permutations (one-line notation, w = (w(1),...,w(n)) ) ----------

def length(w):
    """Coxeter length = number of inversions."""
    n = len(w)
    return sum(1 for i in range(n) for j in range(i+1, n) if w[i] > w[j])

def rmul_t(w, a, b):
    """w * t_{ab}  (1-indexed a<b): swap the entries in POSITIONS a and b."""
    v = list(w)
    v[a-1], v[b-1] = v[b-1], v[a-1]
    return tuple(v)

def refl_length(w):
    """absolute (reflection) length = n - (number of cycles)."""
    n = len(w); seen = [False]*n; cyc = 0
    for i in range(n):
        if not seen[i]:
            cyc += 1
            j = i
            while not seen[j]:
                seen[j] = True
                j = w[j]-1
    return n - cyc

# ---------- polynomials as dict: exponent tuple -> int coefficient ----------

def pmul(f, g):
    h = {}
    for e1, c1 in f.items():
        for e2, c2 in g.items():
            e = tuple(a+b for a, b in zip(e1, e2))
            h[e] = h.get(e, 0) + c1*c2
    return {e: c for e, c in h.items() if c}

def padd(f, g):
    h = dict(f)
    for e, c in g.items():
        h[e] = h.get(e, 0) + c
    return {e: c for e, c in h.items() if c}

def divided_difference(f, j, N):
    """d_j f = (f - s_j f)/(x_j - x_{j+1}),  j is 1-indexed."""
    i, k = j-1, j
    out = {}
    for e, c in f.items():
        p, q = e[i], e[k]
        if p == q:
            continue
        base = list(e)
        if p > q:
            # x_j^p x_{j+1}^q - x_j^q x_{j+1}^p  over (x_j - x_{j+1})
            #  = sum_{r=0}^{p-q-1} x_j^{p-1-r} x_{j+1}^{q+r}
            sgn, hi, lo = 1, p, q
        else:
            sgn, hi, lo = -1, q, p
        for r in range(hi-lo):
            ee = list(base)
            ee[i] = hi-1-r
            ee[k] = lo+r
            t = tuple(ee)
            out[t] = out.get(t, 0) + sgn*c
    return {e: c for e, c in out.items() if c}

# ---------- Schubert polynomials for S_n, in N >= n variables ----------

def schubert_table(n, N=None):
    """dict: permutation of S_n -> Schubert polynomial in N variables."""
    if N is None:
        N = n
    w0 = tuple(range(n, 0, -1))
    e0 = [0]*N
    for i in range(n-1):
        e0[i] = n-1-i
    tab = {w0: {tuple(e0): 1}}
    # BFS downward in length
    cur = [w0]
    for d in range(length(w0), 0, -1):
        nxt = []
        for w in cur:
            for j in range(1, n):
                v = rmul_t(w, j, j+1)
                if length(v) == length(w)-1 and v not in tab:
                    tab[v] = divided_difference(tab[w], j, N)
                    nxt.append(v)
        cur = nxt
    return tab

def descent_path(w):
    """sequence j_1, j_2, ... of simple reflections with l(w s_{j1} s_{j2} ...)
    dropping by 1 each time, ending at identity."""
    path = []
    cur = w
    while length(cur) > 0:
        for j in range(1, len(w)):
            v = rmul_t(cur, j, j+1)
            if length(v) == length(cur)-1:
                path.append(j); cur = v; break
    return path

def extract_coeff(P, w, N):
    """coefficient of S_w in the Schubert expansion of homogeneous P, deg P = l(w).
    Apply d_{j_1}, d_{j_2}, ... along a descent path of w; result is the constant."""
    f = P
    for j in descent_path(w):
        f = divided_difference(f, j, N)
    if not f:
        return 0
    assert set(f.keys()) == {tuple([0]*N)}, ("not a constant", f)
    return f[tuple([0]*N)]

def structure_constants(n, N=None):
    """c[(u,v,w)] for u,v,w in S_n with l(u)+l(v)=l(w)."""
    if N is None:
        N = n+1
    tab = schubert_table(n, N)
    perms = list(permutations(range(1, n+1)))
    byl = {}
    for w in perms:
        byl.setdefault(length(w), []).append(w)
    C = {}
    prodcache = {}
    for u in perms:
        for v in perms:
            d = length(u)+length(v)
            if d not in byl:
                continue
            P = pmul(tab[u], tab[v])
            for w in byl[d]:
                c = extract_coeff(P, w, N)
                if c:
                    C[(u, v, w)] = c
    return C, tab
