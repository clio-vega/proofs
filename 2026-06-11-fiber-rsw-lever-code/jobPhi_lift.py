"""
The REAL lever (post step-law-tautology): the engine's exact (1+z)-lift Phi(z).

  Phi(z) = <s_lambda, (p_1^2 + z p_2)^m> = sum_k C(m,k) chi_k z^k         (chi_k = chi^lambda(2^k 1^{2m-2k}))
         = sum_r C(m,r) 2^r R_r (1+z)^{m-r},   R_r = <s_lambda, p_2^{m-r} e_2^r>.

Coarse Newton locus (Prop 3.2):  e = min_k v2(C(m,k) chi_k);  reduce Phi/2^e mod 2.
Claim: on every TIE,  Phi/2^e ≡ z^{j0} (1+z)^g  (mod 2)  with g >= 2  -> coarse box B* even.

This script:
  (1) verifies Phi/2^e ≡ z^{j0}(1+z)^g mod 2 and extracts (j0, g) for ALL lambda |- 2m, m<=MMAX;
      reports g>=2 on every tie (Prop 3.2) and pushes m as far as feasible.
  (2) GAP (i) coarse->sharp: compares coarse box B* (support of Phi/2^e mod2) against the SHARP
      Newton locus J* (argmin_k k + 2 v2(C(m,k) M_k)); tabulates when |B*| != |J*| and whether
      B* is always an even box exactly when J* is.
  (3) GAP (ii) the g>=2 mechanism: on ties chi^lambda(2^m)=R_0 is even, so the r=0 layer (1+z)^m
      dies mod 2; we locate the FIRST layer r* that survives mod 2 in the (1+z)-basis and test
      whether the surviving low-order structure forces g>=2 (the e_2-mod-2 wall).
"""
import sys, math, csv, os
from collections import Counter, defaultdict
import symfunc as sf
from job1_tie_census import partitions, M_vector, v2, chi_b_vector, analyze

RESULTS = "/home/clio/projects/code/results"
os.makedirs(RESULTS, exist_ok=True)

# -------------------------------------------------- F2 polynomial helpers (support as frozenset of exps)
def poly_div_test_shifted_power(supp, m):
    """Given support set `supp` of an F2 polynomial in z, test if it equals z^{j0}(1+z)^g
    for some j0>=0, g>=0.  Returns (j0, g) if so, else None.
    (1+z)^g over F2 has support = {submasks of g} (Lucas).  So z^{j0}(1+z)^g support =
    {j0 + s : s submask of g}."""
    if not supp:
        return None
    j0 = min(supp)
    shifted = sorted(x - j0 for x in supp)
    # shifted must be exactly the submask set of some g; the max element is g (all bits),
    # and submask set has size 2^{popcount(g)}.
    g = max(shifted)
    subs = set()
    s = g
    while True:
        subs.add(s)
        if s == 0:
            break
        s = (s - 1) & g
    return (j0, g) if set(shifted) == subs else None

def phi_coeffs(lam, m):
    """integer coefficients [chi_0*C(m,0), ..., chi_m*C(m,m)] of Phi(z) in the z-monomial basis."""
    chi = chi_b_vector(lam, m)
    return [math.comb(m, k) * chi[k] for k in range(m+1)]

def R_vector(lam, m):
    """R_r = <s_lambda, p_2^{m-r} e_2^r> for r=0..m (exact)."""
    p2 = sf.p_one(2)
    e2 = sf.e_n(2)
    sl = sf.schur_p(lam)
    Rs = []
    for r in range(m+1):
        f = sf.mul(sf.power(p2, m-r), sf.power(e2, r))
        Rs.append(int(sf.inner(sl, f)))
    return Rs

def coarse_box(lam, m):
    """e=min_k v2(c_k); support of Phi/2^e mod 2; (j0,g) if it is a shifted power of (1+z)."""
    c = phi_coeffs(lam, m)
    vals = [v2(ck) for ck in c]          # None where chi_k=0
    fin = [v for v in vals if v is not None]
    if not fin:
        return None
    e = min(fin)
    supp = frozenset(k for k in range(m+1) if vals[k] == e)
    jg = poly_div_test_shifted_power(supp, m)
    return dict(e=e, supp=supp, jg=jg)

def sharp_Jstar(lam, m):
    """J* = argmin_k k + 2 v2(C(m,k) M_k)  (the real pi-adic Newton locus)."""
    Ms = M_vector(lam, m)
    val = []
    for k in range(m+1):
        ck = math.comb(m, k) * Ms[k]
        vv = v2(ck)
        val.append(None if vv is None else k + 2*vv)
    fin = [v for v in val if v is not None]
    mu = min(fin)
    J = frozenset(k for k in range(m+1) if val[k] == mu)
    return J, mu

# -------------------------------------------------- (1)+(2): coarse box + coarse-vs-sharp
def run_box(MMAX):
    print(f"=== (1) coarse box Phi/2^e ≡ z^j0 (1+z)^g, and (2) coarse vs sharp J*  (m<={MMAX}) ===")
    rows = []
    g_on_ties = Counter()
    not_shifted_power = 0
    ties = 0; g_ge2 = 0
    coarse_size_eq_sharp = 0; total = 0
    coarse_box_even_iff_sharp_even = 0
    size_mismatch = []
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            total += 1
            cb = coarse_box(lam, m)
            J, mu = sharp_Jstar(lam, m)
            if cb is None:
                continue
            jg = cb['jg']
            is_tie = (len(J) >= 2)
            if jg is None:
                not_shifted_power += 1
            else:
                j0, g = jg
            Bsize = len(cb['supp'])
            Jsize = len(J)
            if Bsize == Jsize:
                coarse_size_eq_sharp += 1
            else:
                size_mismatch.append((m, tuple(lam), Bsize, Jsize, jg))
            # even-box agreement
            cb_even = (Bsize % 2 == 0)
            sh_even = (Jsize % 2 == 0)
            if cb_even == sh_even:
                coarse_box_even_iff_sharp_even += 1
            if is_tie:
                ties += 1
                if jg is not None:
                    g_on_ties[jg[1]] += 1
                    if jg[1] >= 2:
                        g_ge2 += 1
            rows.append(dict(m=m, lam='|'.join(map(str, lam)),
                             e=cb['e'], coarse_supp='|'.join(map(str, sorted(cb['supp']))),
                             g=(jg[1] if jg else -1), j0=(jg[0] if jg else -1),
                             Bsize=Bsize, Jstar='|'.join(map(str, sorted(J))), Jsize=Jsize,
                             is_tie=int(is_tie)))
    print(f"  total shapes: {total};  Phi/2^e NOT a shifted power of (1+z): {not_shifted_power}")
    print(f"  TIES (|J*|>=2): {ties};  with coarse g>=2: {g_ge2}  ({'ALL' if g_ge2==ties else 'NOT ALL'})")
    print(f"  coarse-g distribution on ties: {dict(sorted(g_on_ties.items()))}")
    print(f"  |B*| == |J*| (coarse size = sharp size): {coarse_size_eq_sharp}/{total}")
    print(f"  coarse box even <=> sharp box even:      {coarse_box_even_iff_sharp_even}/{total}")
    if size_mismatch:
        print(f"  ** {len(size_mismatch)} shapes with |B*| != |J*| (coarse coarser/finer than sharp):")
        for r in size_mismatch[:15]:
            print(f"       m={r[0]} lam={r[1]}  |B*|={r[2]} |J*|={r[3]}  (j0,g)={r[4]}")
    with open(f"{RESULTS}/phi-coarse-box.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"  wrote results/phi-coarse-box.csv ({len(rows)} rows)")
    return size_mismatch

# -------------------------------------------------- (3): the g>=2 mechanism (e_2 mod-2 layer)
def run_layers(MMAX):
    print(f"\n=== (3) g>=2 mechanism: the (1+z)-layer that survives mod 2 on ties  (m<={MMAX}) ===")
    # Phi = sum_r C(m,r) 2^r R_r (1+z)^{m-r}.  Coeff of (1+z)^{m-r} is C(m,r) 2^r R_r, v2 >= r.
    # On a tie, R_0 = chi^lambda(2^m) is even, so the (1+z)^m layer dies mod 2^? — find the
    # smallest r with v2(C(m,r) 2^r R_r) minimal, i.e. the layer that sets the leading 2-adic content.
    first_layer = Counter()
    g_vs_firstlayer = Counter()
    R0_parity_on_ties = Counter()
    tie_examples = []
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            J, mu = sharp_Jstar(lam, m)
            if len(J) < 2:
                continue
            Rs = R_vector(lam, m)
            cb = coarse_box(lam, m)
            if cb is None or cb['jg'] is None:
                continue
            g = cb['jg'][1]
            # 2-adic content of each (1+z)-layer coefficient
            layer_v = []
            for r in range(m+1):
                coef = math.comb(m, r) * (2**r) * Rs[r]
                layer_v.append(v2(coef))   # None if Rs[r]==0
            fin = [(r, vv) for r, vv in enumerate(layer_v) if vv is not None]
            emin = min(vv for _, vv in fin)
            rstar = min(r for r, vv in fin if vv == emin)   # first minimal layer
            first_layer[rstar] += 1
            g_vs_firstlayer[(rstar, g)] += 1
            R0_parity_on_ties[Rs[0] % 2] += 1
            if len(tie_examples) < 6:
                tie_examples.append((m, tuple(lam), g, rstar, [Rs[r] for r in range(m+1)], layer_v))
    print(f"  R_0=chi^lambda(2^m) parity on ties (0=even):  {dict(R0_parity_on_ties)}")
    print(f"  first minimal (1+z)-layer r* on ties:         {dict(sorted(first_layer.items()))}")
    print(f"  (r*, coarse g) joint distribution on ties:    {dict(sorted(g_vs_firstlayer.items()))}")
    print(f"  worked tie examples (m, lam, g, r*, R_vec, layer_v2):")
    for ex in tie_examples:
        print(f"     m={ex[0]} lam={ex[1]} g={ex[2]} r*={ex[3]} R={ex[4]} layerv2={ex[5]}")

# -------------------------------------------------- cross-check R_r vs chi (binomial transform)
def crosscheck_R(MMAX=6):
    print("=== cross-check: Phi via chi_k == Phi via (1+z)-lift R_r ===")
    import sympy as sp
    z = sp.Symbol('z')
    bad = 0; tot = 0
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            chi = chi_b_vector(lam, m)
            Rs = R_vector(lam, m)
            lhs = sum(math.comb(m, k) * chi[k] * z**k for k in range(m+1))
            rhs = sum(math.comb(m, r) * 2**r * Rs[r] * (1+z)**(m-r) for r in range(m+1))
            tot += 1
            if sp.expand(lhs - rhs) != 0:
                bad += 1
                if bad <= 3:
                    print(f"  MISMATCH lam={lam}")
    print(f"  {tot-bad}/{tot} match")

if __name__ == "__main__":
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    crosscheck_R(5)
    run_box(MMAX)
    run_layers(MMAX)
