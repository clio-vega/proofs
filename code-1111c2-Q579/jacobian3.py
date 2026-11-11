"""Q579 v3: valley initial surface h0 = |i-c|+|j-c| - 2c (c = centre), frozen boundary.
Centre columns are raised many times -> tests the telescoping across repeated updates.
Also: F with arbitrary positive rational coefficients per octahedron (an inhomogeneous,
'weighted' recurrence) must satisfy the SAME identity; non-replacement rules must not."""
import random
from fractions import Fraction
import sympy as sp
from jacobian2 import D, add, mul, div, scal, nbrs

def valley_sequence(L, rng):
    c = (L - 1) // 2
    h = {(i, j): abs(i-c) + abs(j-c) - 2*c for i in range(L) for j in range(L)}
    interior = [(i, j) for i in range(1, L-1) for j in range(1, L-1)]
    seq = []
    while True:
        cand = [p for p in interior if all(h[q] == h[p] + 1 for q in nbrs(p))]
        if not cand: return seq, h
        p = rng.choice(cand); h[p] += 2; seq.append(p)

def check(L, seed, kind):
    rng = random.Random(seed)
    seq, hfin = valley_sequence(L, rng)
    cols = sorted(hfin)
    x0 = {p: Fraction(rng.randint(1, 30), rng.randint(1, 30)) for p in cols}
    x = {p: D(x0[p], {p: Fraction(1)}) for p in cols}
    for p in seq:
        o = x[p]; E, W, N, S = (x[q] for q in nbrs(p))
        a, b = Fraction(rng.randint(1, 9), rng.randint(1, 9)), Fraction(rng.randint(1, 9), rng.randint(1, 9))
        F = add(scal(a, mul(E, W)), scal(b, mul(N, S)))
        if kind == 'octahedron': x[p] = div(add(mul(E, W), mul(N, S)), o)
        elif kind == 'weighted': x[p] = div(F, o)
        elif kind == 'ctrl_plus_old': x[p] = add(div(F, o), o)
        elif kind == 'ctrl_sqrt_like': x[p] = div(F, mul(o, o))
    J = sp.Matrix([[sp.Rational(x[p].g[q].numerator, x[p].g[q].denominator) if q in x[p].g else 0
                    for q in cols] for p in cols])
    pred = sp.Integer(-1)**len(seq)
    for p in cols:
        r = x[p].v / x0[p]; pred *= sp.Rational(r.numerator, r.denominator)
    return len(seq), max(seq.count(p) for p in set(seq)), J.det() == pred

if __name__ == '__main__':
    for L in (3, 5, 7, 9):
        for kind in ('octahedron', 'weighted', 'ctrl_plus_old', 'ctrl_sqrt_like'):
            n = 3 if L < 9 else 1
            res = [check(L, s, kind) for s in range(n)]
            print(f"L={L} {kind:15s} moves={res[0][0]:4d} max-updates/column={res[0][1]:2d} "
                  f"identity: {sum(r[2] for r in res)}/{n}", flush=True)
