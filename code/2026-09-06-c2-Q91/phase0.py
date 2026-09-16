"""Phase 0: two code-disjoint engines must agree on R_e(t), symbolically in t."""
import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
import engine                      # SHAPE side
from bead import R_bead, open_interval_load, trim, t   # BEAD side

assert engine.t is not sp.Symbol   # sanity: both use sympy Symbol('t')
T = sp.Symbol('t')


def shape_R(lam, e):
    return engine.op_R({trim(lam): 1}, e)


def compare(e_list, nmax, **kw):
    agree = 0; total = 0; bad = []
    for e in e_list:
        for n in range(0, nmax + 1):
            for lam in engine.parts_of(n) if n else ((),):
                lhs = {k: sp.expand(v) for k, v in shape_R(lam, e).items()}
                rhs = {k: sp.expand(v) for k, v in R_bead(lam, e, **kw).items()}
                total += 1
                if lhs == rhs:
                    agree += 1
                else:
                    bad.append((e, lam, lhs, rhs))
    return agree, total, bad


if __name__ == "__main__":
    print("=== MAIN CHECK: Conjecture N, e=2,3,4, |lam| <= 7, t symbolic ===")
    a, tot, bad = compare([2, 3, 4], 7)
    print(f"  agreement: {a}/{tot}")
    for b in bad[:5]:
        print("   MISMATCH", b[0], b[1])
        print("     shape:", b[2])
        print("     bead :", b[3])

    print()
    print("=== WITNESS: need #(M cap (b,b+e)) >= 2 somewhere ===")
    best = (-1, None, None)
    for e in (2, 3, 4):
        for n in range(0, 8):
            for lam in engine.parts_of(n) if n else ((),):
                N, b = open_interval_load(lam, e)
                if N > best[0]:
                    best = (N, lam, e, b)
    print("  max open-interval occupancy over the tested range:", best)
    # enumerate all (lam,e,b) with load >= 2 at the smallest |lam|
    hits = []
    for e in (2, 3, 4):
        for n in range(0, 8):
            for lam in engine.parts_of(n) if n else ((),):
                N, b = open_interval_load(lam, e)
                if N >= 2:
                    hits.append((n, lam, e, b, N))
    hits.sort()
    print(f"  #(lam,e) with a bead move of load >= 2: {len(hits)}; smallest:")
    for h in hits[:6]:
        print("   ", h)

    print()
    print("=== does (-t-1)^2 / t^2 actually appear in a coefficient? ===")
    seen2 = []
    for e in (3, 4):
        for n in range(0, 8):
            for lam in engine.parts_of(n) if n else ((),):
                for mu, c in R_bead(lam, e).items():
                    if sp.degree(sp.Poly(c, T)) >= 2:
                        seen2.append((lam, e, mu, c))
    print(f"  coefficients of t-degree >= 2: {len(seen2)}; e.g. {seen2[:3]}")
