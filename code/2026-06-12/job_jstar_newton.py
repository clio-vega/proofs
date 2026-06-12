"""
Job A.2 -- 2-adic Newton polygon of the half-polynomial for every tie shape.

For a tie shape lambda |- 2m with single-parity J* (parity p in {0,1}):
    half-poly coeffs  c_t = C(m, 2t+p) * M_{2t+p}     (t = 0,1,2,...)
Newton-polygon points  P_t = (t, v2(c_t)).
The slope-(-1) edge = the t minimizing  w(t) = t + v2(c_t);  these are exactly J*/.. since
    val(2t+p) = (2t+p) + 2 v2(c_t) = p + 2 w(t).
Edge contributors  E = { t : w(t) = min }.   |E| = |J*|  (must be even when >=2).

Tabulated per shape:
  * the on-edge t's and the j = 2t+p they correspond to
  * reduced edge coefficients (c_t / 2^{v2}) mod 2   -- expect all 1 (odd)
  * edge length = |J*|, evenness
  * t-space offsets of the edge and whether they form a 2-adic box (generators)
Summary: do ALL edges have even length?  Are reduced coeffs always odd?
"""
import json
from collections import Counter
from job_jstar_engine import C, v2

def two_adic_box_gens(offs):
    """A 2-adic box has generators that are distinct powers of 2 (memory: j0+{sum 2^a}).
    Then the offset set is an XOR-closed subspace.  Return the XOR basis (the generators)
    if offs is such a box, else None."""
    offs = sorted(offs)
    if offs[0] != 0:
        return None
    S = set(offs)
    k = len(S)
    if k & (k-1) != 0:
        return None
    # XOR closure
    for a in offs:
        for b in offs:
            if (a ^ b) not in S:
                return None
    # extract a minimal XOR basis (Gaussian elimination over GF(2))
    basis = []
    for x in sorted(offs):
        cur = x
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
    box = {0}
    for d in basis:
        box |= {x ^ d for x in box}
    return sorted(basis) if box == S else None

def main():
    ties = json.load(open("results-jstar-ties.json"))
    all_even = True
    all_odd_coeff = True
    box_ok = 0; box_bad = 0
    edgelen_census = Counter()
    examples = []
    for rec in ties:
        m = rec["m"]; M = rec["M"]; p = rec["Jstar"][0] % 2
        # build half-poly coeffs
        coeffs = {}
        t = 0
        while 2*t + p <= m:
            j = 2*t + p
            cj = C(m, j) * M[j]
            if cj != 0:
                coeffs[t] = cj
            t += 1
        w = {t: t + v2(c) for t, c in coeffs.items()}
        wmin = min(w.values())
        edge = sorted(t for t, ww in w.items() if ww == wmin)
        L = len(edge)
        edgelen_census[L] += 1
        if L % 2 != 0 and L >= 2:
            all_even = False
        # reduced coeffs mod 2
        red = []
        for t in edge:
            c = coeffs[t]
            red.append((c >> v2(c)) & 1)
        if any(r != 1 for r in red):
            all_odd_coeff = False
        # offsets
        offs = [t - edge[0] for t in edge]
        gens = two_adic_box_gens(offs)
        if gens is None:
            box_bad += 1
        else:
            box_ok += 1
        # keep a few illustrative examples (incl. all |J*|>=4)
        if L >= 4 or len(examples) < 8:
            examples.append({
                "lam": tuple(rec["lam"]), "m": m, "parity": p,
                "edge_t": edge, "edge_j": [2*t+p for t in edge],
                "wmin": wmin, "red_coeff_mod2": red, "L": L,
                "t_offsets": offs, "box_gens": gens,
            })

    print("=== Newton-polygon slope(-1) edge, all %d tie shapes ===" % len(ties))
    print("edge-length census:", dict(sorted(edgelen_census.items())))
    print("ALL edge lengths even (when >=2):", all_even)
    print("ALL reduced edge coeffs ODD (==1 mod 2):", all_odd_coeff)
    print("edge t-offsets form 2-adic box: %d ok / %d bad" % (box_ok, box_bad))
    print("\n--- examples (all |J*|>=4 plus first few) ---")
    seen4 = 0
    for e in examples:
        tag = "  *|J*|>=4*" if e["L"] >= 4 else ""
        if e["L"] >= 4:
            seen4 += 1
            if seen4 > 12:
                continue
        if e["L"] < 4 and examples.index(e) >= 8:
            continue
        print(f"lam={e['lam']} m={e['m']} par={e['parity']}: edge j={e['edge_j']} "
              f"redmod2={e['red_coeff_mod2']} L={e['L']} t_off={e['t_offsets']} "
              f"gens={e['box_gens']}{tag}")

if __name__ == "__main__":
    main()
