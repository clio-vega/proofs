"""
Job B (2026-06-15): is the FIBER value G_lambda(zeta_d) = <s_lambda,
h_1^e prod (h_2 + e_2 zeta^{2k-1})> a zeta_d-shift of a principal
specialisation  s_lambda(1,q,...,q^{k-1})|_{q=zeta_d}?

CODE.md Job B / live-probes #6.  NOTE: this is the same probe run decisively on
2026-06-11 as `jobC_rsw.py` (clean d=2-only verdict).  This run RE-VERIFIES that
verdict from a clean reimplementation AND adds the two specific cuts CODE.md asks
for that the 06-11 run did not isolate:
  (a) the FINITE principal specs at the named variable counts k=2 and k=d (not just
      a blind N=1..15 sweep);
  (b) an explicit integer-shift search  G ?= zeta_d^s * princ_spec  for s in a band,
      and, on a miss, the characterisation of the ratio G/princ_spec.

All exact (sympy), evaluated at exact zeta_d.  Cross-checked against fiber_engine
(which is itself cross-checked vs the Murnaghan-Nakayama G_gaussian and chi(2^m)).
"""
import sympy as sp
from fiber_engine import (G_value, princ_spec_value, princ_spec_poly,
                          zeta, partitions, n_lambda, hook_lengths)

q = sp.Symbol('q')

# ---- fake degree f^lambda(q) = q^{n(lam)} (q;q)_n / prod (1-q^{h})  (= STABLE
#      principal spec s_lam(1,q,q^2,...), the "q-analogue of dimension") ----
def fake_degree_value(lam, d):
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    hooks = hook_lengths(lam)
    num = sp.prod([1 - q**k for k in range(1, n+1)])
    den = sp.prod([1 - q**h for h in hooks.values()])
    expr = q**n_lambda(lam) * num / den
    poly = sp.cancel(expr)
    assert sp.denom(poly) == 1, (lam,)
    return sp.simplify(poly.subs(q, zeta(d)))

def is_unit_shift(G, S, d):
    """Return integer s in [-2n,2n] with G == zeta_d^s * S, else None.
       Also returns the ratio G/S (simplified) when S != 0."""
    if sp.simplify(S) == 0:
        return ('Szero', None)
    ratio = sp.simplify(G / S)
    if sp.simplify(G) == 0:
        return ('Gzero', ratio)
    zt = zeta(d)
    for s in range(-2*d, 2*d+1):
        if sp.simplify(ratio - zt**s) == 0:
            return (s, ratio)
    return (None, ratio)

def run():
    DS = [2, 3, 4]
    NMAX = 14
    # tally: per d, per candidate-k, how many shapes are an exact zeta-shift
    from collections import Counter, defaultdict
    hit = defaultdict(Counter)     # (d, kind) -> Counter{ 'hit', 'Gzero', 'miss' }
    ratios = defaultdict(list)     # (d, kind) -> list of (lam, ratio) for misses
    first_miss = {}                # (d, kind) -> first divergent lam

    for n in range(2, NMAX+1):
        for lam in partitions(n):
            lam = tuple(lam)
            for d in DS:
                G = sp.simplify(G_value(lam, d))
                # candidate specialisations
                cands = {}
                cands['fake(N=inf)'] = fake_degree_value(lam, d)
                for k in sorted(set([2, d, len(lam)])):
                    cands[f'k={k}'] = princ_spec_value(lam, k, d)
                for kind, S in cands.items():
                    s, ratio = is_unit_shift(G, S, d)
                    if s == 'Szero':
                        # principal spec vanished; only a hit if G also 0
                        if sp.simplify(G) == 0:
                            hit[(d, kind)]['hit'] += 1
                        else:
                            hit[(d, kind)]['miss'] += 1
                            if (d, kind) not in first_miss:
                                first_miss[(d, kind)] = (lam, 'S=0,G!=0')
                    elif s == 'Gzero':
                        # G vanished, S didn't: only a "hit" if we count 0 as 0*shift; no.
                        hit[(d, kind)]['Gzero'] += 1
                        if (d, kind) not in first_miss:
                            first_miss[(d, kind)] = (lam, f'G=0,S!=0 ratio={ratio}')
                    elif isinstance(s, int):
                        hit[(d, kind)]['hit'] += 1
                        hit[(d, kind)][f'shift={s}'] += 1
                    else:
                        hit[(d, kind)]['miss'] += 1
                        ratios[(d, kind)].append((lam, ratio))
                        if (d, kind) not in first_miss:
                            first_miss[(d, kind)] = (lam, f'ratio={ratio}')

    print("="*72)
    print(f"Job B: G_lambda(zeta_d) vs zeta_d^s * principal-spec,  all lam |- n<=({NMAX})")
    print("="*72)
    total = sum(1 for n in range(2,NMAX+1) for _ in partitions(n))
    for d in DS:
        nshapes = sum(1 for n in range(2,NMAX+1) for _ in partitions(n))
        print(f"\n--- d = {d}  (zeta_{d}, {nshapes} shapes) ---")
        kinds = sorted({k for (dd,k) in hit if dd==d})
        for kind in kinds:
            c = hit[(d,kind)]
            h = c['hit']; gz = c['Gzero']; ms = c['miss']
            shifts = {k:v for k,v in c.items() if k.startswith('shift=')}
            print(f"  {kind:12s}: hit {h:4d}  Gzero {gz:3d}  miss {ms:4d}   shifts {dict(sorted(shifts.items()))}")
            if (d,kind) in first_miss and ms>0:
                print(f"               first miss: lam={first_miss[(d,kind)][0]}  {first_miss[(d,kind)][1]}")
        # characterise the miss ratios for the fake-degree candidate (the d=2 winner)
        rr = ratios[(d,'fake(N=inf)')]
        if rr:
            vals = sorted(set(str(sp.nsimplify(r)) for _,r in rr[:40]))
            print(f"  ratio spread (fake, sample): {vals[:10]}")

    # ---- the specific d=4 vanisher (2,2) ----
    print("\n" + "="*72)
    print("d=4 vanisher lambda=(2,2):")
    for d in [4]:
        G = sp.simplify(G_value((2,2), d))
        print(f"  G_(2,2)(i) = {G}")
        for k in [2,4]:
            S = sp.simplify(princ_spec_value((2,2), k, d))
            print(f"  s_(2,2)(1,...,q^{k-1})|q=i = {S}")
        print(f"  fake_(2,2)(i) = {sp.simplify(fake_degree_value((2,2), 4))}")

if __name__ == "__main__":
    run()
