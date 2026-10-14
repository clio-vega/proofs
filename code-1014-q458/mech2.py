"""Mechanism 2: no sympy, no symmetric-function library.
Laurent polynomials as dicts {exponent-vector : int}. Builds edot_r by the
elementary-symmetric recursion over the 2n variables x_1,1/x_1,...,x_n,1/x_n,
then reads off [x^0] and the whole k=1 sum.  Shares no code path with core.py."""
from math import comb
from itertools import combinations

def edot_constant_terms(n):
    """Return list ct[r] = [x^0] e_r(x_1,1/x_1,...,x_n,1/x_n), 0<=r<=2n,
    by exact Laurent arithmetic (dict of exponent tuples)."""
    vars_ = []
    for i in range(n):
        v = [0]*n; v[i] = 1;  vars_.append(tuple(v))
        w = [0]*n; w[i] = -1; vars_.append(tuple(w))
    e = [dict() for _ in range(2*n+1)]
    e[0] = {tuple([0]*n): 1}
    for v in vars_:
        for r in range(2*n, 0, -1):
            src = e[r-1]
            if not src: continue
            dst = e[r]
            for mon, c in src.items():
                nm = tuple(a+b for a, b in zip(mon, v))
                dst[nm] = dst.get(nm, 0) + c
                if dst[nm] == 0: del dst[nm]
    zero = tuple([0]*n)
    return [e[r].get(zero, 0) for r in range(2*n+1)]

def edot_constant_terms_bruteforce(n):
    """Third reading: literally enumerate r-subsets of the 2n exponent vectors."""
    ex = []
    for i in range(n):
        ex.append((i, +1)); ex.append((i, -1))
    out = []
    for r in range(2*n+1):
        c = 0
        for S in combinations(range(2*n), r):
            tot = [0]*n
            for s in S: tot[ex[s][0]] += ex[s][1]
            if all(v == 0 for v in tot): c += 1
        out.append(c)
    return out

def c_1n(n, ell):
    """c_{(1^n)}(t) as dict {t-exponent: coeff}, via mechanism 2."""
    M = ell+2; N = 2*M
    ct = edot_constant_terms(n)
    def CT(r): return ct[r] if 0 <= r <= 2*n else 0
    out = {}
    for y in range(-(n+4), n+5):
        d = CT(n-N*y) - CT(n+2-N*y)
        if d:
            out[M*y*y-y] = out.get(M*y*y-y, 0) + d
    return {k: v for k, v in out.items() if v != 0}

if __name__ == '__main__':
    print("=== census check: [x^0]edot_r  vs  binom(n,r/2) ===")
    allok = True
    for n in range(1, 9):
        ct = edot_constant_terms(n)
        pred = [comb(n, r//2) if r % 2 == 0 else 0 for r in range(2*n+1)]
        ok = (ct == pred)
        extra = ""
        if n <= 6:
            bf = edot_constant_terms_bruteforce(n)
            extra = f"  bruteforce_agrees={bf == ct}"
            ok = ok and (bf == ct)
        print(f"  n={n}: dict-Laurent==binom  {ct == pred}{extra}")
        allok = allok and ok
    print("CENSUS ALL OK =", allok)

    print("\n=== c_{(1^n)} mechanism-2 vs closed form ===")
    import sys; sys.path.insert(0, '/home/clio/projects/proofs/code-1014-q458')
    import sympy as sp
    from census import closed_form, predicted_degree
    t = sp.Symbol('t')
    ok = fail = 0
    for ell in (0,1,2,3):
        for n in range(1, 13):
            m2 = c_1n(n, ell)
            m2e = sp.expand(sum(v*t**k for k, v in m2.items())) if m2 else sp.Integer(0)
            cf = closed_form(n, ell)
            same = sp.expand(m2e - cf) == 0
            deg = max(m2) if m2 else None
            pd = predicted_degree(n, ell)[0] if n % 2 == 0 else None
            dok = (deg == pd)
            print(f"  ell={ell} n={n}: match={same}  deg={deg} predicted={pd} degree_law={dok}")
            if same and dok: ok += 1
            else: fail += 1
    print(f"MECH2 ok={ok} fail={fail}")
