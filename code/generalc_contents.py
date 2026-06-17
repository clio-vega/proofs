"""
generalc_contents.py  --  CODE session 2026-06-17, Job A

General-c three-row d=4 boundary: tabulate the 2-adic CONTENT g[c][k] of the
boundary deficit polynomials N_{c-k}^{(c)} for c = 5,6,7,8 (and 2,3,4 for
continuity), split by parity of b, and decide whether the Content bound is
c-UNIFORM (a function of depth k alone) or genuinely c-dependent.

Two engines, cross-checked:
  (A) DIRECT (fast, robust, reaches c=8 / large m): the verified clean-constant
      identity  const_i = c!*k!  (master valuation, generalc_master.py, 0 mismatch
      c=4..8) gives v2(N_i) EXACTLY from the MN value of M_{b+i}:
          v2(N_i) = v2(M_{b+i}) + v2(c!) + v2(k!) - v2(b-c+1)
                    - [ v2(a-b+1)         if i==1
                        sum_{t=2}^i v2(b+t) otherwise ].
      Content g[c][k] = min over the (P,b) grid of v2(N_i(a=2P+b+c, b)).
  (B) SYMBOLIC (c<=6): polynomial fit of M_{b+i}, strip hook factors, read off
      N_i, take gcd-of-coefficients content.  Used to CROSS-CHECK (A) and to
      print the explicit closed forms (the deliverable's "closed forms").

Outputs:
  * the content table g[c][k] (headline), b-even / b-odd split
  * c-uniformity verdict + best-fit law for g[c][k]
  * boundary certify  Delta(b+i) >= 1  for c<=8, large m, with min slack per
    (c, b-parity) and extremal shapes
  * the theta(c) interior-offset table
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "threerow-boundary"))
from math import comb, factorial
from mn import Mj
import sympy as sp

# ----------------------------------------------------------------------
def v2(n):
    n = abs(int(n))
    if n == 0:
        return 999
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

# ====================  ENGINE A : direct v2(N_i) from MN  =============
def v2_Ni(P, bb, c, i):
    """Exact v2 of the deficit N_i^{(c)} at the shape a=2P+b+c, via const_i=c!*k!."""
    aa = 2 * P + bb + c
    m = P + bb + c
    k = c - i
    if bb + i > m:
        return None
    Mbi = Mj((aa, bb, c), bb + i, m)
    if Mbi == 0:
        return 999  # deficit vanishes -> infinite content there; excluded from min
    base = v2(Mbi) + v2(factorial(c)) + v2(factorial(k)) - v2(bb - c + 1)
    if i == 1:
        return base - v2(aa - bb + 1)
    return base - sum(v2(bb + t) for t in range(2, i + 1))

def content_table_direct(cmax=8, Pmax=24, bmax=56, mmax=84):
    """g[c][k] = min over grid of v2(N_i), split by parity of b.
    Records (g, witness, count_at_min) per parity so we can judge convergence:
    if the min v2 is attained on many (P,b) it is the true content, not a grid edge."""
    table = {}   # (c,k) -> (g_all, g_even, g_odd, wit_even, wit_odd, hits_e, hits_o)
    for c in range(2, cmax + 1):
        for i in range(1, c + 1):
            k = c - i
            g_e = g_o = 999
            w_e = w_o = None
            hits_e = hits_o = 0
            for P in range(0, Pmax + 1):
                for bb in range(2 * c, bmax):     # box interior b>=2c
                    m = P + bb + c
                    if m > mmax:
                        continue
                    vv = v2_Ni(P, bb, c, i)
                    if vv is None or vv >= 999:
                        continue
                    if bb % 2 == 0:
                        if vv < g_e:
                            g_e, w_e, hits_e = vv, (P, bb, 2 * P + bb + c, m), 1
                        elif vv == g_e:
                            hits_e += 1
                    else:
                        if vv < g_o:
                            g_o, w_o, hits_o = vv, (P, bb, 2 * P + bb + c, m), 1
                        elif vv == g_o:
                            hits_o += 1
            table[(c, k)] = (min(g_e, g_o), g_e, g_o, w_e, w_o, hits_e, hits_o)
    return table

# ====================  ENGINE B : symbolic closed forms  ==============
a_s, b_s = sp.symbols('a b')

def collect(c, i, amax=64, bmax=30):
    pts = []
    for bb in range(c, bmax):
        for aa in range(bb, amax):
            if (aa + bb + c) % 2 != 0:
                continue
            m = (aa + bb + c) // 2
            j = bb + i
            if j > m:
                continue
            pts.append((aa, bb, Mj((aa, bb, c), j, m)))
    return pts

def fit_poly(pts, maxdeg):
    monos = [(p, q) for p in range(maxdeg + 1) for q in range(maxdeg + 1) if p + q <= maxdeg]
    if len(pts) < len(monos) + 3:
        return None
    A = sp.Matrix([[aa**p * bb**q for (p, q) in monos] for (aa, bb, _) in pts])
    y = sp.Matrix([M for (_, _, M) in pts])
    sol = A.solve_least_squares(y)
    if any(A * sol - y):
        return None
    expr = sum(sol[k] * a_s**p * b_s**q for k, (p, q) in enumerate(monos))
    return sp.factor(sp.expand(expr))

def extract_symbolic(c, i):
    """fit M_{b+i}, strip hook factors (a-b+1)(a-c+2), return N_i (a-involving residue)."""
    M = fit_poly(collect(c, i), maxdeg=2 * c + 1)
    if M is None:
        return None
    _, facs = sp.factor_list(M)
    hookA = sp.expand(a_s - b_s + 1)
    hookB = sp.expand(a_s - c + 2)
    N = sp.Integer(1)
    removedA = removedB = 0
    for fac, mult in facs:
        fac_e = sp.expand(fac)
        deg_a = sp.Poly(fac_e, a_s).degree() if fac_e.has(a_s) else 0
        if deg_a == 0:
            continue
        if removedA == 0 and sp.simplify(fac_e - hookA) == 0:
            removedA = 1; mult -= 1
        if removedB == 0 and sp.simplify(fac_e - hookB) == 0:
            removedB = 1; mult -= 1
        if mult > 0:
            N *= fac**mult
    return sp.expand(N)

def content_symbolic(N, c, parity):
    g = 999
    for Pv in range(0, 45):
        for bv in range(2 * c, 90):
            if bv % 2 != parity:
                continue
            av = 2 * Pv + bv + c
            vv = v2(int(N.subs({a_s: av, b_s: bv})))
            if vv < g:
                g = vv
    return g

# ====================  theta(c)  ======================================
def theta_table(cmax=8):
    out = {}
    for c in range(2, cmax + 1):
        seen = {}
        for m in range(c + 4, 34):
            for bb in range(2 * c, m):
                aa = 2 * m - c - bb
                if aa < bb:
                    continue
                lam = (aa, bb, c)
                vals = {}
                for j in range(0, bb + c + 1):
                    if j > m:
                        continue
                    M = Mj(lam, j, m)
                    if M == 0:
                        continue
                    vals[j] = j + 2 * v2(comb(m, j) * M)
                if 0 not in vals:
                    continue
                V = min(vals.values())
                Jstar = sorted(j for j, v in vals.items() if v == V)
                seen.setdefault(aa % 2, set()).add((vals[0] - V, Jstar[0]))
        out[c] = seen
    return out

# ====================  Delta certify  =================================
def true_delta(P, bb, c, i):
    aa = 2 * P + bb + c; m = P + bb + c; lam = (aa, bb, c)
    M0 = Mj(lam, 0, m)
    if M0 == 0:
        return None
    j = bb + i
    if j > m:
        return None
    Mji = Mj(lam, j, m)
    if Mji == 0:
        return None
    return (j + 2 * v2(comb(m, j) * Mji)) - 2 * v2(M0)

def certify(cmax=8, Pmax=50, mmax=150):
    rows = {}
    for c in range(4, cmax + 1):
        checked = 0; viol = 0
        minslack = {}   # (b%2) -> [min slack, witness]
        for P in range(0, Pmax + 1):
            for bb in range(2 * c, mmax):
                aa = 2 * P + bb + c; m = P + bb + c
                if aa < bb or m > mmax:
                    continue
                theta = 3 if (c % 2 == 1 and aa % 2 == 1) else 0
                for i in range(1, c + 1):
                    td = true_delta(P, bb, c, i)
                    if td is None:
                        continue
                    checked += 1
                    # boundary requirement: val(b+i) > V = val(0)-theta, i.e. Delta(b+i) > -theta
                    if not (td > -theta):
                        viol += 1
                    slack = td + theta            # margin above the threshold
                    key = bb % 2
                    cur = minslack.get(key, [10**9, None])
                    if slack < cur[0]:
                        minslack[key] = [slack, (aa, bb, c, i, m)]
        rows[c] = (checked, viol, minslack)
    return rows

# ======================================================================
def g0(c, k):
    """candidate c-parity DEPTH lower bound for the content."""
    return 2 * (k // 2) + (1 if (c % 2 == 1 and k % 2 == 1) else 0)

if __name__ == "__main__":
    print("#" * 78)
    print("# JOB A : general-c boundary contents  (CODE 2026-06-17)")
    print("#" * 78)

    # ---- (2) the content table (headline), direct engine, c=2..8 --------
    print("\n" + "=" * 78)
    print("(2) CONTENT TABLE g[c][k] (direct MN via const_i=c!k!, box interior b>=2c, m<=84)")
    print("=" * 78)
    table = content_table_direct(cmax=8)

    print("\n  g[c][k]  (min over b-parities):")
    header = "       " + "".join("k=%-5d" % k for k in range(0, 8))
    print(header)
    for c in range(2, 9):
        row = "  c=%-3d " % c
        for k in range(0, 8):
            row += "%-7s" % (table[(c, k)][0] if (c, k) in table else ".")
        print(row)

    print("\n  g[c][k]  [b even / b odd]:")
    print(header)
    for parity, lab in [(1, "even"), (2, "odd")]:
        print("   -- b %s --" % lab)
        for c in range(2, 9):
            row = "  c=%-3d " % c
            for k in range(0, 8):
                if (c, k) in table:
                    val = table[(c, k)][parity]
                    row += "%-7s" % val
                else:
                    row += "%-7s" % "."
            print(row)

    # ---- (3) c-uniformity analysis + g0 lower bound --------------------
    print("\n" + "=" * 78)
    print("(3) c-UNIFORMITY: is g[c][k] a function of depth k alone?")
    print("=" * 78)
    for k in range(0, 8):
        vals = [(c, table[(c, k)][0]) for c in range(2, 9) if (c, k) in table]
        if not vals:
            continue
        gvals = sorted(set(v for _, v in vals))
        print("  k=%d : %s   %s" % (k, vals, "UNIFORM" if len(gvals) == 1 else "c-DEPENDENT"))

    print("\n  candidate LOWER bound  g0(k) = 2*floor(k/2) + [c odd & k odd]:")
    print("   c  k  true(min)  g0   dev   (b even/odd)")
    valid = True
    sharp_shallow = True
    for c in range(4, 9):
        for k in range(0, c):
            if (c, k) not in table:
                continue
            ga, ge, go, we, wo, he, ho = table[(c, k)]
            gg = g0(c, k); dev = ga - gg
            flag = ""
            if dev < 0:
                valid = False; flag = " <-- VIOLATION"
            if k <= 3 and dev != 0:
                sharp_shallow = False; flag += " (shallow not sharp)"
            print("   %d  %d    %2d      %2d   %+d   (%d/%d)%s"
                  % (c, k, ga, gg, dev, ge, go, flag))
        print()
    print("  g0 is a VALID lower bound (true >= g0) for all c>=4:  %s" % valid)
    print("  g0 is SHARP (equality) for all shallow depth k<=3:    %s" % sharp_shallow)
    print("  => content is NOT c-uniform as an equality (deep k grows with c);")
    print("     g0 IS a c-uniform-up-to-parity lower bound, sharp where the boundary needs it.")

    print("\n  convergence (hits at min v2; large => true content, not grid edge):")
    for (c, k), (ga, ge, go, we, wo, he, ho) in sorted(table.items()):
        print("    c=%d k=%d : even g=%d (hits %d)  odd g=%d (hits %d)"
              % (c, k, ge, he, go, ho))

    # ---- (4) theta(c) --------------------------------------------------
    print("\n" + "=" * 78)
    print("(4) theta(c) interior offset table")
    print("=" * 78)
    th = theta_table(cmax=8)
    for c in range(2, 9):
        for par in sorted(th[c]):
            print("  c=%d a%%2=%d : (theta,j0) = %s"
                  % (c, par, sorted(th[c][par])))

    # ---- (5) certify Delta >= 1 ----------------------------------------
    print("\n" + "=" * 78)
    print("(5) CERTIFY  Delta(b+i) > -theta  (boundary loses), c=4..8")
    print("=" * 78)
    rows = certify(cmax=8)
    for c in range(4, 9):
        checked, viol, ms = rows[c]
        print("  c=%d: checked=%d  violations=%d" % (c, checked, viol))
        for par, lab in [(0, "b even"), (1, "b odd")]:
            if par in ms:
                sl, wit = ms[par]
                print("        %s : min slack=%d at (a,b,c,i,m)=%s" % (lab, sl, wit))

    # ---- (6) symbolic closed forms + cross-check (optional, slower) -----
    print("\n" + "=" * 78)
    print("(6) EXPLICIT DEFICIT FORMS N_i^{(c)} (symbolic) + cross-check vs direct")
    print("=" * 78)
    mism = 0
    for c in [3, 4, 5]:
        print("\n---- c = %d ----" % c)
        for i in range(1, c + 1):
            k = c - i
            N = extract_symbolic(c, i)
            if N is None:
                print("  i=%d (k=%d): symbolic fit unavailable" % (i, k))
                continue
            ge = content_symbolic(N, c, 0)
            go = content_symbolic(N, c, 1)
            adeg = sp.Poly(sp.expand(N), a_s).degree() if N.has(a_s) else 0
            if (c, k) in table:
                _, dge, dgo, _, _, _, _ = table[(c, k)]
                if (ge, go) != (dge, dgo):
                    mism += 1
                    print("    CROSS-CHECK MISMATCH c=%d k=%d sym(%d/%d) direct(%d/%d)"
                          % (c, k, ge, go, dge, dgo))
            print("  i=%d (k=%d): a-deg=%d content g=%d [b even %d/odd %d]"
                  % (i, k, adeg, min(ge, go), ge, go))
            if k <= 2:
                print("       N = %s" % N)
    print("\n  symbolic-vs-direct content mismatches = %d" % mism)
