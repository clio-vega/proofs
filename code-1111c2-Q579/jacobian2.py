"""Q579 v2: no height cap (Lipschitz + frozen boundary bounds the surface), so columns are
updated REPEATEDLY.  Jacobian computed exactly by forward-mode dual numbers over Fraction
(independent of sympy.diff); determinant by exact sympy on a Fraction matrix."""
import random, sys
from fractions import Fraction
import sympy as sp

def nbrs(p):
    i, j = p
    return (i+1, j), (i-1, j), (i, j+1), (i, j-1)

def full_sequence(L, rng):
    h = {(i, j): (i + j) % 2 for i in range(L) for j in range(L)}
    interior = [(i, j) for i in range(1, L-1) for j in range(1, L-1)]
    seq = []
    while True:
        cand = [p for p in interior if all(h[q] == h[p] + 1 for q in nbrs(p))]
        if not cand: return seq, h
        p = rng.choice(cand); h[p] += 2; seq.append(p)

class D:  # dual number: value + gradient (dict col->Fraction)
    __slots__ = ('v', 'g')
    def __init__(s, v, g): s.v, s.g = v, g
def add(a, b):
    g = dict(a.g)
    for k, w in b.g.items(): g[k] = g.get(k, 0) + w
    return D(a.v + b.v, g)
def mul(a, b):
    g = {k: w * b.v for k, w in a.g.items()}
    for k, w in b.g.items(): g[k] = g.get(k, 0) + w * a.v
    return D(a.v * b.v, g)
def div(a, b):
    g = {k: w / b.v for k, w in a.g.items()}
    for k, w in b.g.items(): g[k] = g.get(k, 0) - w * a.v / b.v**2
    return D(a.v / b.v, g)
def scal(c, a): return D(c * a.v, {k: c * w for k, w in a.g.items()})

RULES = {
  'octahedron': lambda o, E, W, N, S: div(add(mul(E, W), mul(N, S)), o),
  'ctrl_plus_old': lambda o, E, W, N, S: add(div(add(mul(E, W), mul(N, S)), o), o),
  'ctrl_square': lambda o, E, W, N, S: div(add(mul(E, W), mul(N, S)), mul(o, o)),
  'ctrl_weighted': lambda o, E, W, N, S: div(add(mul(E, W), scal(Fraction(2), mul(N, S))), o),
}

def check(L, seed, rn):
    rng = random.Random(seed)
    seq, hfin = full_sequence(L, rng)
    cols = sorted(hfin)
    x0 = {p: Fraction(rng.randint(1, 30), rng.randint(1, 30)) for p in cols}
    x = {p: D(x0[p], {p: Fraction(1)}) for p in cols}
    for p in seq:
        x[p] = RULES[rn](x[p], *(x[q] for q in nbrs(p)))
    J = sp.Matrix([[sp.Rational(x[p].g.get(q, 0).numerator, x[p].g.get(q, 0).denominator)
                    if q in x[p].g else 0 for q in cols] for p in cols])
    d = J.det()
    pred = sp.Integer(-1)**len(seq)
    for p in cols:
        r = x[p].v / x0[p]; pred *= sp.Rational(r.numerator, r.denominator)
    maxupd = max(seq.count(p) for p in set(seq)) if seq else 0
    return len(seq), maxupd, d == pred

if __name__ == '__main__':
    for L in (3, 4, 5, 6, 7):
        for rn in RULES:
            res = [check(L, s, rn) for s in range(4 if L < 7 else 2)]
            print(f"L={L} {rn:14s} moves={res[0][0]:3d} max-updates-per-column={res[0][1]} "
                  f"identity holds: {sum(r[2] for r in res)}/{len(res)}", flush=True)
