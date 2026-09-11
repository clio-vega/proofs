"""
DERIVER (engine E1): build G(mu) directly from prop:reflect of
proofs/2026-09-10-c2-Q140-local-identity-certificate-family.tex.

VERTEX SET, fixed here in writing (instrument-reports-on-its-referent):
  Z-Maya convention.  B = B(mu) = { mu_i - i : i>=1 }  (mu_i = 0 for i > ell),
  so B = {b_1,...,b_ell} cup {-k : k > ell},  b_i = mu_i - i strictly decreasing.
  A VERTEX is a pair (b, u) with b in B a bead, u notin B a hole, u < b.
  prop:count asserts #V = n.  ASSERTED below.

EDGES, verbatim from prop:reflect:
 (I)  beads c < b, hole b' < c, e := b+c-b' is a hole  ==>  edge {(b,b'),(c,b')},
      relation  w(c,b') = -t^{|B cap (c,b)|} w(b,b').
 (II) bead b, holes b' < e < b, c := b'+e-b is a bead  ==>  edge {(b,b'),(b,e)},
      relation  w(b,e) = -t^{1+|B cap (b',e)|} w(b,b').
Each edge carries label -t^h; we store (u,v,h) meaning w(v) = -t^h w(u).
"""
from itertools import combinations

def maya(mu):
    """Return (ell, beads_list b_1..b_ell, is_bead(v), min_hole = -ell)."""
    ell = len(mu)
    bs = [mu[i] - (i+1) for i in range(ell)]
    bset = set(bs)
    def is_bead(v):
        return v <= -ell-1 or v in bset
    return ell, bs, is_bead

def vertices(mu):
    ell, bs, is_bead = maya(mu)
    V = []
    for b in bs:
        for u in range(-ell, b):          # min hole is -ell; nothing below it is a hole
            if not is_bead(u):
                V.append((b, u))
    return V

def count_beads_open(is_bead, lo, hi, ell):
    """#(B cap (lo,hi)) -- open interval.  lo >= -ell-? ; all args here are > -ell-1."""
    return sum(1 for v in range(lo+1, hi) if is_bead(v))

def edges(mu, use_I=True, use_II=True, strict=True):
    ell, bs, is_bead = maya(mu)
    V = set(vertices(mu))
    E = []
    holes_below = {b: [u for u in range(-ell, b) if not is_bead(u)] for b in bs}
    if use_I:
        for c, b in combinations(sorted(bs), 2):       # c < b
            for bp in holes_below[c]:                  # hole b' < c
                e = b + c - bp
                if not is_bead(e):
                    h = count_beads_open(is_bead, c, b, ell)
                    E.append(((b, bp), (c, bp), h, 'I'))
    if use_II:
        for b in bs:
            for bp, e in combinations(holes_below[b], 2):   # b' < e < b
                c = bp + e - b
                ok = is_bead(c) if strict else (is_bead(c) or True)
                if ok:
                    h = 1 + count_beads_open(is_bead, bp, e, ell)
                    E.append(((b, bp), (b, e), h, 'II'))
    for (x, y, h, ty) in E:
        assert x in V and y in V, (mu, x, y, ty)
    return E

# ---------- graph utilities ----------
def components_and_bipartite(V, E):
    adj = {v: [] for v in V}
    for (x, y, h, ty) in E:
        adj[x].append((y, h, +1)); adj[y].append((x, h, -1))
    colour, comp = {}, {}
    comps = []
    for s in V:
        if s in colour: continue
        cid = len(comps); comps.append([]); colour[s] = 0; comp[s] = cid
        stack = [s]; comps[cid].append(s)
        while stack:
            x = stack.pop()
            for (y, h, sg) in adj[x]:
                if y not in colour:
                    colour[y] = 1 - colour[x]; comp[y] = cid
                    comps[cid].append(y); stack.append(y)
    odd = set()
    for (x, y, h, ty) in E:
        if colour[x] == colour[y]:
            odd.add(comp[x])
    return comps, comp, colour, odd

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,) + rest
