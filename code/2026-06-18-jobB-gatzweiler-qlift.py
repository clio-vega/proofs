"""
2026-06-18 Job B (probe 15): Gatzweiler-Krattenthaler q-lift of the v2(prod C) walls.

Thesis: my boundary/NL_c/Lemma-F2 walls are v2(product-of-binomials) inequalities.
v2 lives at the cyclotomic factor Phi_2(q)=q+1 of the q-analogue, via the classical
dictionary
    v2( C(n,k) )  =  sum_{j>=1} mult_{Phi_{2^j}} ( [n choose k]_q ),
and  mult_{Phi_d}([n;k]_q) = #carries when adding k+(n-k) in base d  (Kummer, q-form).
So each v2 wall is the q=-1 shadow of a q-binomial cyclotomic-multiplicity inequality
-- exactly the kind of object G-K (2502.06032) prove positive.

This script (no browsing; arithmetic only):
 (A) verifies the dictionary v2(C) = sum_j mult_{Phi_{2^j}}[n;k]_q numerically;
 (B) lifts Lemma F2 (Q>=6 consecutive ints, drop 2, product stays even) to the
     q-Pochhammer and checks (q+1) | remaining (the Phi_2 shadow) AND whether the
     remaining product is cyclotomic-nonneg (stronger G-K-type positivity);
 (C) lifts NL_2 to q-binomials and checks the cyclotomic-multiplicity inequality.
"""
import sympy as sp

q = sp.symbols('q')


def v2(n):
    n = abs(int(n))
    if n == 0:
        return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def qint(n):
    return sp.expand(sum(q**i for i in range(n))) if n > 0 else sp.Integer(1)


def qbinom(n, k):
    if k < 0 or k > n:
        return sp.Integer(0)
    num = sp.Integer(1)
    den = sp.Integer(1)
    for i in range(k):
        num *= qint(n - i)
        den *= qint(i + 1)
    return sp.cancel(num / den)


def phi_mult(poly, d):
    """multiplicity of cyclotomic Phi_d(q) in poly (a polynomial in q)."""
    poly = sp.Poly(sp.expand(poly), q)
    phi = sp.Poly(sp.cyclotomic_poly(d, q), q)
    m = 0
    while True:
        quo, rem = sp.div(poly, phi)
        if rem == 0 and poly.degree() >= phi.degree():
            poly = quo
            m += 1
        else:
            break
    return m


def carries_base(k, nk, d):
    """#carries adding k+nk in base d."""
    c = 0
    carry = 0
    while k or nk or carry:
        s = (k % d) + (nk % d) + carry
        carry = 1 if s >= d else 0
        c += carry
        k //= d
        nk //= d
    return c


# ---------- (A) dictionary v2(C) = sum_j mult_{Phi_{2^j}} ----------
def check_dictionary(N=22):
    bad = 0
    tested = 0
    for n in range(1, N):
        for k in range(0, n + 1):
            qb = qbinom(n, k)
            lhs = v2(sp.binomial(n, k))                       # 2-adic valuation
            jmax = n.bit_length() + 1
            rhs = sum(phi_mult(qb, 2**j) for j in range(1, jmax + 1))  # sum_j Phi_{2^j} mult
            kum = carries_base(k, n - k, 2)                   # Kummer: carries in BASE 2
            tested += 1
            if not (lhs == rhs == kum):
                bad += 1
                if bad <= 5:
                    print("   MISMATCH n=%d k=%d v2=%d sumPhi=%d kummer=%d"
                          % (n, k, lhs, rhs, kum))
    return bad, tested


# ---------- (B) Lemma F2 q-lift ----------
def qpoch_drop(m, Q, drop):
    """q-product of [m+1]..[m+Q] with factors at positions in `drop` removed."""
    prod = sp.Integer(1)
    for s in range(1, Q + 1):
        if s in drop:
            continue
        prod *= qint(m + s)
    return sp.expand(prod)


def cyclotomic_nonneg(poly):
    """True if poly factors into cyclotomics with nonneg mult only (it always does
    for products of q-integers; we report the coeff-nonnegativity of poly itself)."""
    p = sp.Poly(sp.expand(poly), q)
    return all(c >= 0 for c in p.all_coeffs())


def check_F2():
    from itertools import combinations
    allphi2 = True
    allnonneg = True
    worst = None
    for Q in (6, 7, 8):
        for m in range(0, 8):
            for drop in combinations(range(1, Q + 1), 2):
                rem = qpoch_drop(m, Q, set(drop))
                m2 = phi_mult(rem, 2)            # Phi_2 = q+1 -> v2>=1 shadow
                # integer v2 of the remaining product of consecutive ints
                iv = sum(v2(m + s) for s in range(1, Q + 1) if s not in drop)
                if m2 < 1:
                    allphi2 = False
                if m2 != iv:
                    # Phi_2 mult counts #even factors; integer v2 counts total 2s.
                    # they differ (v2>=Phi_2 mult); we only need Phi_2 mult>=1.
                    pass
                if not cyclotomic_nonneg(rem):
                    allnonneg = False
                if worst is None or m2 < worst[0]:
                    worst = (m2, Q, m, drop, iv)
    return allphi2, allnonneg, worst


# ---------- (C) NL_2 q-lift: v2 C(F+3,j) + v2(j(j-1)(j-2)(j-3)) >= v2(F)+1 ----------
def check_NL2_qlift(Fmax=18, jmax=18):
    """integer NL_2 and its Phi_2-multiplicity q-lift."""
    bad_int = 0
    bad_q = 0
    tested = 0
    for F in range(4, Fmax):
        for j in range(4, min(jmax, F + 3)):
            lhs_int = v2(sp.binomial(F + 3, j)) + v2(j * (j - 1) * (j - 2) * (j - 3))
            rhs_int = v2(F) + 1
            if lhs_int < rhs_int:
                bad_int += 1
            # q-lift: Phi_2 mult of [F+3 choose j]_q + Phi_2 mult of q-falling [j]_4
            qb = qbinom(F + 3, j)
            qfall = qint(j) * qint(j - 1) * qint(j - 2) * qint(j - 3)
            lhs_q = phi_mult(qb, 2) + phi_mult(sp.expand(qfall), 2)
            rhs_q = phi_mult(qint(F), 2) + 1
            if lhs_q < rhs_q:
                bad_q += 1
            tested += 1
    return bad_int, bad_q, tested


if __name__ == "__main__":
    print("=" * 76)
    print("(A) dictionary  v2(C(n,k)) = sum_j mult_{Phi_{2^j}}[n;k]_q = Kummer")
    bad, tested = check_dictionary(22)
    print("    tested=%d  bad=%d  -> %s" % (tested, bad,
          "DICTIONARY HOLDS" if bad == 0 else "FAILS"))

    print("\n(B) Lemma F2 q-lift: Q>=6 consec, drop 2, remaining q-Pochhammer")
    allphi2, allnonneg, worst = check_F2()
    print("    (q+1)=Phi_2 divides remaining for ALL drops: %s" % allphi2)
    print("    remaining q-product has nonneg coeffs (cyclotomic-positive): %s" % allnonneg)
    print("    worst Phi_2 mult = %d at Q=%d m=%d drop=%s (integer v2=%d)"
          % (worst[0], worst[1], worst[2], worst[3], worst[4]))

    print("\n(C) NL_2 q-lift: Phi_2 mult inequality vs integer v2 inequality")
    bi, bq, t = check_NL2_qlift()
    print("    integer NL_2 violations: %d/%d" % (bi, t))
    print("    Phi_2-multiplicity q-lift violations: %d/%d" % (bq, t))
    print("    -> %s" % ("q-LIFT HOLDS (NL_2 is the Phi_2 shadow of a q-binomial ineq)"
                         if bq == 0 else "q-lift FAILS"))
