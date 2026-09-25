"""AKO's side: stacked-ribbon decomposition, heights, widths, predicted weight."""
from cyl import *
from weight import geom, w

def intermediates(R):
    "all cylindric sub-diagrams of R with the same inner boundary (as hi'-tuples)"
    x, y, lo, hi = R.x, R.y, R.lo, R.hi
    out = []
    def rec(i, cur):
        if i == y:
            if cur[-1] <= cur[0] + x: out.append(tuple(cur))
            return
        start = max(lo[i]-1, cur[-1] if cur else lo[i]-1)
        for v in range(start, hi[i]+1):
            rec(i+1, cur+[v])
    rec(0, [])
    return out

def region(R, hlo, hhi):
    "the skew region between sub-diagrams with outer boundaries hlo and hhi"
    return Cyl(R.x, R.y, [v+1 for v in hlo], list(hhi))

def size_of(R, h):
    return sum(h[i] - R.lo[i] + 1 for i in range(R.y))

def is_ribbon(S):
    g = geom(S)
    return g['conn'] and g['sq'] == 0

def is_loop_ribbon(S):
    return S.size() == S.x + S.y and is_ribbon(S)

def decompose(R):
    """R = L_1 u ... u L_l u F, loop ribbons peeled from the OUTER boundary.
    returns (l, F) or None if R is not a stacked ribbon."""
    n = R.x + R.y
    r = R.size()
    if r == 0: return None
    if r < n:
        return (0, R) if is_ribbon(R) else None
    for h in intermediates(R):
        if size_of(R, h) != r - n: continue
        S = region(R, h, R.hi)
        if not is_loop_ribbon(S): continue
        if r == n:                      # F itself is the loop ribbon: pure
            return (0, R)
        inner = Cyl(R.x, R.y, R.lo, list(h))
        sub = decompose(inner)
        if sub is not None:
            l, F = sub
            return (l+1, F)
    return None

def ako_weight(R):
    """AKO prop:stackedRibbon / thm:cylindricMNRule, as literally stated:
       ht(R) = l*y + ht(F),  ht(ribbon) = #rows - 1 = #vertical edges,
       wd = x if pure (F is a loop ribbon) else 1."""
    d = decompose(R)
    if d is None: return 0, None
    l, F = d
    pure = is_loop_ribbon(F)
    htF = F.x + F.y if False else geom(F)['vert']      # vertical edges of F
    ht = l*R.y + (R.y if pure else htF)
    wd = R.x if pure else 1
    return (-1)**ht * wd, dict(l=l, F=F, pure=pure, htF=htF, ht=ht, wd=wd)

def ako_weight_corrected(R):
    """same, but ht(R) = sum of constituent heights + l   (one extra vertical edge per
    junction between consecutive constituents)."""
    d = decompose(R)
    if d is None: return 0, None
    l, F = d
    pure = is_loop_ribbon(F)
    htF = R.y if pure else geom(F)['vert']
    ht = l*R.y + htF + l
    wd = R.x if pure else 1
    return (-1)**ht * wd, dict(l=l, F=F, pure=pure, ht=ht, wd=wd)

def mn_rule(D, weight=None):
    """sum over chains lo-1 = h_0 < h_1 < ... < h_l = hi with |step j| = nu_j,
    each step a stacked ribbon; product of step weights.  Returns {nu: value}."""
    if weight is None: weight = lambda R: ako_weight_corrected(R)[0]
    inter = intermediates(D)
    empty = tuple(v-1 for v in D.lo)
    assert empty in inter, (empty, 'empty diagram not an intermediate')
    sz = {h: size_of(D, h) for h in inter}
    # cache step weights
    wt = {}
    def stepw(a, b):
        if (a,b) not in wt:
            try: R = region(D, a, b)
            except AssertionError: wt[(a,b)] = 0; return 0
            wt[(a,b)] = weight(R)
        return wt[(a,b)]
    out = {}
    target = tuple(D.hi)
    def rec(h, acc, hist):
        if h == target:
            nu = tuple(sorted(hist, reverse=True))
            out[nu] = out.get(nu, 0) + acc
            return
        for h2 in inter:
            if h2 == h: continue
            if any(h2[i] < h[i] for i in range(D.y)): continue
            d = sz[h2] - sz[h]
            if d <= 0: continue
            if hist and d > hist[-1]: continue     # parts weakly decreasing: nu a partition, inner first
            sw = stepw(h, h2)
            if sw == 0: continue
            rec(h2, acc*sw, hist+[d])
    rec(empty, 1, [])
    return out

def all_decompositions(R):
    "every way of peeling loop ribbons from the outer boundary -> list of (l, F)"
    n = R.x + R.y; r = R.size(); out = []
    if r == 0: return out
    if r < n:
        return [(0, R)] if is_ribbon(R) else []
    hits = []
    for h in intermediates(R):
        if size_of(R, h) != r - n: continue
        S = region(R, h, R.hi)
        if is_loop_ribbon(S): hits.append(h)
    if r == n:
        return [(0, R)] if hits or is_loop_ribbon(R) else []
    for h in hits:
        for (l, F) in all_decompositions(Cyl(R.x, R.y, R.lo, list(h))):
            out.append((l+1, F))
    return out
