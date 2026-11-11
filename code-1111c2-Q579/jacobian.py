"""Q579: Jacobian of octahedron-recurrence move sequences on stepped surfaces.

State: a height function h on an L x L grid of columns (i,j), with |h(p)-h(q)|=1 for
grid neighbours and h = (i+j) mod 2 initially; boundary columns are frozen.  A move at an
interior column p is allowed when all 4 neighbours have h = h(p)+1; it sets
  x_p <- F(x_N, x_S, x_E, x_W) / x_p,    h(p) <- h(p)+2.
Octahedron:      F = x_E x_W + x_N x_S.
Claim A (geometric): det d(x_final)/d(x_init) = (-1)^N prod_p x_p^final / x_p^init.
Claim B (tropical) : for y_p <- max(yE+yW, yN+yS) - y_p, det = (-1)^N a.e.
Controls: rules NOT of the form x_new x_old = F(others) must FAIL Claim A.
"""
import sys, random, itertools
import sympy as sp
from fractions import Fraction

def moves_to_pyramid(L, rng, max_moves=None):
    """Random admissible move sequence, raising the surface as far as possible."""
    h = {(i, j): (i + j) % 2 for i in range(L) for j in range(L)}
    interior = [(i, j) for i in range(1, L - 1) for j in range(1, L - 1)]
    seq = []
    while True:
        cand = [p for p in interior
                if all(h[(p[0]+a, p[1]+b)] == h[p] + 1 for a, b in ((1,0),(-1,0),(0,1),(0,-1)))
                and h[p] + 2 <= min(p[0], p[1], L-1-p[0], L-1-p[1]) * 2 + (p[0]+p[1]) % 2 + 1]
        if not cand or (max_moves and len(seq) >= max_moves):
            return seq
        p = rng.choice(cand); h[p] += 2; seq.append(p)

def nbrs(p):
    i, j = p
    return (i+1, j), (i-1, j), (i, j+1), (i, j-1)

def run(L, seq, rule, xs):
    x = dict(xs)
    for p in seq:
        E, W, N, S = (x[q] for q in nbrs(p))
        x[p] = rule(x[p], E, W, N, S)
    return x

RULES = {
    'octahedron': lambda o, E, W, N, S: (E*W + N*S) / o,
    'ctrl_plus_old': lambda o, E, W, N, S: (E*W + N*S) / o + o,      # not x_new*x_old = F
    'ctrl_square': lambda o, E, W, N, S: (E*W + N*S) / o**2,         # x_new*x_old^2 = F
    'ctrl_dep_old': lambda o, E, W, N, S: (E*W + N*S + o) / o,       # F depends on x_old
}

def geometric(L, seed, rule_name):
    rng = random.Random(seed)
    seq = moves_to_pyramid(L, rng)
    cols = sorted((i, j) for i in range(L) for j in range(L))
    syms = {p: sp.Symbol(f"x_{p[0]}_{p[1]}", positive=True) for p in cols}
    fin = run(L, seq, RULES[rule_name], syms)
    J = sp.Matrix([[sp.diff(fin[p], syms[q]) for q in cols] for p in cols])
    # evaluate exactly at a random rational point (identity check by exact specialisation
    # at several points; the symbolic det is too slow beyond L=4)
    ok_all = True
    for trial in range(3):
        pt = {syms[p]: sp.Rational(rng.randint(1, 9), rng.randint(1, 9)) for p in cols}
        d = J.subs(pt).det()
        pred = (-1)**len(seq) * sp.prod([fin[p].subs(pt) / syms[p].subs(pt) for p in cols])
        ok_all &= sp.simplify(d - pred) == 0
    return len(seq), ok_all

def tropical(L, seed, coef=1):
    """Exact PL Jacobian: pick a rational point, perturb each coordinate by a tiny exact
    epsilon, finite differences are exact inside a linear region (checked by two eps)."""
    rng = random.Random(seed)
    seq = moves_to_pyramid(L, rng)
    cols = sorted((i, j) for i in range(L) for j in range(L))
    rule = lambda o, E, W, N, S: max(E + W, N + S) - coef * o
    y0 = {p: Fraction(rng.randint(-10**6, 10**6), rng.randint(1, 997)) for p in cols}
    f0 = run(L, seq, rule, y0)
    def col(q, eps):
        y = dict(y0); y[q] += eps
        f = run(L, seq, rule, y)
        return [(f[p] - f0[p]) / eps for p in cols]
    M1 = [col(q, Fraction(1, 10**12)) for q in cols]
    M2 = [col(q, Fraction(1, 10**13)) for q in cols]
    if M1 != M2:
        return len(seq), None   # straddled a wall; skip
    d = sp.Matrix(M1).T.det()
    return len(seq), d

if __name__ == '__main__':
    for L in (3, 4, 5):
        for rn in RULES:
            if L == 5 and rn != 'octahedron':
                seeds = [0]
            else:
                seeds = range(3)
            for s in seeds:
                N, ok = geometric(L, s, rn)
                print(f"geom L={L} rule={rn:14s} seed={s} moves={N:3d} det==(-1)^N*prod(xf/xi): {ok}")
    for L in (3, 4, 5, 6, 7):
        for coef in (1, 2):
            dets = []
            for s in range(20):
                N, d = tropical(L, s, coef)
                dets.append((N, d))
            good = [(N, d) for N, d in dets if d is not None]
            match = sum(1 for N, d in good if d == (-1)**N * coef**N)
            unit = sum(1 for N, d in good if abs(d) == 1)
            print(f"trop L={L} coef={coef} usable={len(good)}/20 det==(-1)^N*coef^N: {match}  |det|==1: {unit}  moves={sorted(set(N for N,_ in good))}")
