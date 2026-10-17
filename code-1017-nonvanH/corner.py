"""gap:nonvanH via a bidegree corner.

If P = sum_{i,j} c_{ij} w_a^i w_b^j  and  P|_H (w_b -> w_a^{-1}) = sum_d (sum_{i-j=d} c_{ij}) w_a^d,
then i-j is maximised uniquely at (i,j) = (A, B) := (max i, min j) PROVIDED that corner is
occupied; and then coeff of w_a^{A-B} in P|_H is exactly c_{A,B} /= 0, so P|_H /= 0.

MEASURE: over all (a<b, M >= I, {a,b} cap M = empty), is the corner (max_i, min_j) occupied?
And if so, what does c_{A,B} factor as?
"""
import sympy as sp, itertools
import q484_engine as E
from q484_engine import q, wvars

def bidegree_support(P, wa, wb):
    """P Laurent in wa,wb over the rest. Return dict (i,j)->coeff, cleared of wa,wb denominators."""
    P = sp.cancel(sp.together(P))
    num, den = sp.fraction(P)
    # den must be free of wa, wb for this to be a clean Laurent in wa,wb; check.
    assert sp.degree(sp.Poly(den, wa), wa) == 0 or den.has(wa) is False, f"den has wa: {den}"
    # clear negative powers by multiplying by wa^p wb^p
    num = sp.expand(num)
    pa = sp.Poly(sp.expand(num * wa**50 * wb**50), wa, wb)
    out = {}
    for (i, j), c in zip(pa.monoms(), pa.coeffs()):
        if c != 0:
            out[(i - 50, j - 50)] = sp.cancel(c / den)
    return out

print("=== is the bidegree corner (max_i, min_j) occupied? ===")
for N in (2, 3, 4):
    w = wvars(N)
    nocc = occ = 0
    examples = []
    for a, b in itertools.combinations(range(1, N + 1), 2):
        wa, wb = w[a - 1], w[b - 1]
        for M in E.subsets(N):
            if a in M or b in M: continue
            for I in E.subsets(N):
                if not E.geq(M, I): continue
                A = E.Atilde(N, w, M, I, X='C')
                if A == 0:
                    print(f"  !! Atilde_{{{sorted(M)},{sorted(I)}}} == 0 with M>=I  (converse of tri FAILS)")
                    continue
                supp = bidegree_support(A, wa, wb)
                if not supp: continue
                imax = max(i for (i, j) in supp)
                jmin = min(j for (i, j) in supp)
                if (imax, jmin) in supp:
                    occ += 1
                    if len(examples) < 4:
                        examples.append((a,b,sorted(M),sorted(I),imax,jmin,sp.factor(supp[(imax,jmin)])))
                else:
                    nocc += 1
                    print(f"  CORNER EMPTY: N={N} (a,b)=({a},{b}) M={sorted(M)} I={sorted(I)} "
                          f"imax={imax} jmin={jmin}")
    print(f"  N={N}: corner occupied {occ}, empty {nocc}")
    for ex in examples:
        print(f"     e.g. (a,b)=({ex[0]},{ex[1]}) M={ex[2]} I={ex[3]}  corner=({ex[4]},{ex[5]})")
        print(f"          c = {ex[6]}")
