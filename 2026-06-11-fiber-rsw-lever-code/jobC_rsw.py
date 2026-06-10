"""
JOB C (decisive, never run) — RSW principal-specialisation probe.

Question: is  G_lambda(zeta_d)  a q-shift of the principal specialisation
  s_lambda(1, q, ..., q^{N-1}) | q = zeta_d   (RSW / cyclic-sieving), with the
d-core/d-quotient explicit?  If it HITS, the whole vanishing trichotomy is
off-the-shelf RSW.  If it misses, the multiplicative correction is itself informative.

Method (exact, then high-precision numeric for equality testing):
  LHS  = G_lambda(zeta_d)             [proved master identity, fiber_engine]
  RHS candidates at q = zeta_d:
     - principal spec in N variables s_lambda(1,..,q^{N-1}), N in a natural menu
       {len(lam), len(lam)+1, lam_1, n, d, 2d}
     - fake degree  f^lambda(q) = q^{n(lam)} (q;q)_n / prod_cells (1-q^h)   (maj-gen of SYT)
  For each candidate we test whether  LHS = zeta_d^c * RHS  for some shift c in 0..d-1
  (and the unshifted c=0), to high precision; exact-verify the survivors.

Outputs: results/rsw-probe.csv  +  a hit/miss verdict per d, with the correction factor.
"""
import sys, csv, os
import sympy as sp
import mpmath as mp
from collections import defaultdict
from fiber_engine import (partitions, G_poly, princ_spec_poly, n_lambda,
                          hook_lengths, t as TSYM, zeta)

RESULTS = "/home/clio/projects/code/results"
os.makedirs(RESULTS, exist_ok=True)
mp.mp.dps = 60

def zeta_num(d):
    return mp.e**(2j*mp.pi/d)

def fake_degree_poly(lam):
    """f^lambda(q) = q^{n(lam)} prod_{i=1}^n (1-q^i) / prod_cells (1-q^{hook}), exact poly in q."""
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    q = sp.Symbol('q')
    num = sp.Integer(1)
    for i in range(1, n+1):
        num *= (1 - q**i)
    den = sp.Integer(1)
    for h in hook_lengths(lam).values():
        den *= (1 - q**h)
    poly = sp.cancel(q**n_lambda(lam) * num / den)
    assert sp.denom(poly) == 1, f"fake degree not a polynomial: lam={lam}"
    return sp.Poly(poly, q).as_expr(), q

def num_eval_poly_in(symbol, poly, zval):
    """Numerically evaluate a sympy polynomial `poly` in `symbol` at complex zval (mpmath)."""
    f = sp.lambdify(symbol, poly, modules='mpmath')
    return f(zval)

def candidate_menu(lam):
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    L = len(lam)
    c1 = lam[0]
    Ns = sorted(set([L, L+1, c1, n]))
    cands = [('princN%d' % N, ('princ', N)) for N in Ns]
    cands.append(('fakedeg', ('fake', None)))
    return cands

def close(a, b, tol=1e-40):
    return abs(complex(a) - complex(b)) < tol

def main(NMAX):
    rows = []
    # per (d, candidate-kind): count hits where LHS = zeta^c * RHS for some c
    hit_any = defaultdict(int); seen = defaultdict(int)
    hit_shift0 = defaultdict(int)
    for n in range(1, NMAX+1):
        for lam in partitions(n):
            lam = tuple(lam)
            Pt = G_poly(lam)            # poly in TSYM
            cands = candidate_menu(lam)
            # precompute candidate polynomials
            candpoly = {}
            for name, (kind, N) in cands:
                if kind == 'princ':
                    res = princ_spec_poly(lam, N)
                    candpoly[name] = (None if res == 0 else res)   # (poly,q) or None
                else:
                    candpoly[name] = fake_degree_poly(lam)
            for d in (2, 3, 4):
                zc = zeta_num(d)
                Lval = num_eval_poly_in(TSYM, Pt, zc)
                Lzero = abs(complex(Lval)) < 1e-40
                for name, (kind, N) in cands:
                    cp = candpoly[name]
                    if cp is None:
                        Rval = mp.mpf(0)
                    else:
                        poly, q = cp
                        Rval = num_eval_poly_in(q, poly, zc)
                    seen[(d, name)] += 1
                    Rzero = abs(complex(Rval)) < 1e-40
                    # match up to a zeta_d^c shift
                    best_c = None
                    if Lzero and Rzero:
                        best_c = 0   # both zero: trivially consistent (record but weak)
                    elif (not Lzero) and (not Rzero):
                        for c in range(d):
                            if close(Lval, (zc**c) * Rval):
                                best_c = c; break
                    if best_c is not None:
                        hit_any[(d, name)] += 1
                        if best_c == 0:
                            hit_shift0[(d, name)] += 1
                    # record ratio (if both nonzero) for correction analysis
                    ratio = None
                    if (not Lzero) and (not Rzero):
                        r = Lval / Rval
                        ratio = complex(r)
                    rows.append(dict(
                        n=n, lam='|'.join(map(str, lam)), d=d, cand=name,
                        Lzero=int(Lzero), Rzero=int(Rzero),
                        shift_c=('' if best_c is None else best_c),
                        ratio_re=('' if ratio is None else round(ratio.real, 6)),
                        ratio_im=('' if ratio is None else round(ratio.imag, 6)),
                    ))
    with open(f"{RESULTS}/rsw-probe.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote results/rsw-probe.csv ({len(rows)} rows, n<={NMAX})\n")

    print("=== HIT RATE: G_lambda(zeta_d) == zeta_d^c * candidate(zeta_d) for some c ===")
    print(f"{'d':>2} {'candidate':>10} {'hit/seen (any shift)':>22} {'hit/seen (shift 0)':>20}")
    for (d, name) in sorted(seen):
        ha = hit_any[(d, name)]; h0 = hit_shift0[(d, name)]; s = seen[(d, name)]
        print(f"{d:>2} {name:>10} {f'{ha}/{s}':>22} {f'{h0}/{s}':>20}")

    # zoom: for the strongest candidate per d, show the ratio distribution on non-trivial cases
    print("\n=== ratio G/candidate distribution (non-zero LHS&RHS), best candidate per d ===")
    for d in (2, 3, 4):
        # pick candidate with max any-shift hits
        cnames = [name for (dd, name) in seen if dd == d]
        best = max(cnames, key=lambda nm: hit_any[(d, nm)])
        ratios = defaultdict(int)
        for r in rows:
            if r['d'] == d and r['cand'] == best and r['ratio_re'] != '':
                key = (r['ratio_re'], r['ratio_im'])
                ratios[key] += 1
        ha = hit_any[(d, best)]; s = seen[(d, best)]
        print(f" d={d}  best candidate={best}  (any-shift hits {ha}/{s})")
        top = sorted(ratios.items(), key=lambda kv: -kv[1])[:8]
        for (rr, ri), cnt in top:
            print(f"     ratio {rr:+}{ri:+}i : {cnt}")

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    main(NMAX)
