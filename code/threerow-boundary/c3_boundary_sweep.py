"""
c3_boundary_sweep.py  -- Job A (2026-06-16 code session)

Pin the sharp 2-adic facts behind the c=3 boundary lemma for lambda=(a,b,3),
a+b=2m-3.  We confirm what is TRUE and document what is FALSE, with equality
loci, to large m.  The canonical object is the DIRECT Delta:

    G_j = <s_lambda, e2^j h1^{2(m-j)}> = C(m,j) M_j,
    val(j) = j + 2 v2(G_j) = j + 2( v2 C(m,j) + v2 M_j ),
    Delta(j) := val(j) - val(0)   (parity = j mod 2).

Boundary indices are j = b+1, b+2, b+3.  The lemma needs Delta(b+i) > -theta
(theta = 0 for a even, 3 for a odd) for every box-interior shape.

The two a-linear factors that appear in the closed forms M_{b+i}:
    N2 = a(b+3) - (b^2+2b+3)
    N1 = a b^2 + 5 a b + 6 a - b^3 - b^2 - 4 b - 6
and the abscissa P = m - b - 3.

CLAIM (A): is  v2(N2) >= v2(P+1) - 1  ?   (CODE.md hoped yes)
CLAIM (B): fit the sharp RHS for v2(N1); test the COMBINED direct-Delta ineq.
CLAIM (C): verify the two-factor Lemma F2 premise directly.

Direct Delta formulas (cross-checked vs Murnaghan-Nakayama, 0 mismatch):
    Delta(b+2) = (b+2) - 2 s2(b+3) + 2 K(2(P+1), b+3) + 2 v2(N2) - 2 v2(a-1) - 2 v2(P+2)
    Delta(b+1) = (b-1) - 2 s2(b+1) + 2 K(2(P+2), b+1) + 2 v2(N1) - 2 v2(a-1)
    Delta(b+3) = (b+5) - 2 s2(b+5) + 2 K(2P,     b+5)             - 2 v2(a-1) - 2 v2(P+2)
where K(x,y) = v2 C(x+y, x) = #carries adding x+y in base 2 (Kummer).
"""
from mn import Mj
from math import comb

def v2(n):
    n = abs(n)
    if n == 0:
        return 99
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

def s2(n):
    return bin(n).count('1')

def K(x, y):            # = v2 C(x+y, x), Kummer carry count
    return v2(comb(x + y, x))

def val(lam, j, m):
    M = Mj(lam, j, m)
    if M == 0:
        return None
    return j + 2 * v2(comb(m, j) * M)

def N2_of(a, b):
    return a * (b + 3) - (b**2 + 2*b + 3)

def N1_of(a, b):
    return a*b**2 + 5*a*b + 6*a - b**3 - b**2 - 4*b - 6

def shapes(MMAX):
    """box-interior c=3 shapes: a even b>=6, a odd b>=9."""
    for m in range(8, MMAX + 1):
        for b in range(6, m):
            a = 2*m - 3 - b
            if a < b:
                continue
            P = m - b - 3
            if P < 0:
                continue
            if a % 2 == 0 and b < 6:
                continue
            if a % 2 == 1 and b < 9:
                continue
            yield (a, b, m, P)

# --------------------------------------------------------------------------
# CLAIM (A): v2(N2) vs v2(P+1)
# --------------------------------------------------------------------------
def claim_A(MMAX):
    print("=" * 72)
    print("CLAIM (A): test  v2(N2) >= v2(P+1) - 1   (box-interior, m<=%d)" % MMAX)
    print("=" * 72)
    viol = []        # v2(N2) < v2(P+1)-1
    eqloc = []       # v2(N2) == v2(P+1)-1   (equality of the hoped bound)
    # also compare v2(N2) to v2(P+1) directly, and split a even/odd
    aeven_vN2 = {}
    n = 0
    for (a, b, m, P) in shapes(MMAX):
        n += 1
        vN2 = v2(N2_of(a, b))
        vP1 = v2(P + 1)
        if a % 2 == 0:
            aeven_vN2[vN2] = aeven_vN2.get(vN2, 0) + 1
        if vN2 < vP1 - 1:
            viol.append((a, b, m, vN2, vP1))
        if vN2 == vP1 - 1:
            eqloc.append((a, b, m, vN2, vP1, a % 4, b % 4, P % 4))
    print("shapes scanned: %d" % n)
    print("VIOLATIONS of v2(N2) >= v2(P+1)-1 : %d" % len(viol))
    for x in viol[:12]:
        print("    a=%d b=%d m=%d  v2(N2)=%d v2(P+1)=%d" % x)
    print("a-even distribution of v2(N2):", dict(sorted(aeven_vN2.items())))
    print("[expect a-even => v2(N2)=1 exactly: %s]" %
          (set(aeven_vN2) == {1}))
    # equality-of-bound locus signature
    sig = {}
    for *_, am, bm, pm in eqloc:
        sig[(am, bm, pm)] = sig.get((am, bm, pm), 0) + 1
    print("equality v2(N2)=v2(P+1)-1 occurs %d times; (a%%4,b%%4,P%%4) signature:"
          % len(eqloc))
    for k in sorted(sig):
        print("    a%%4=%d b%%4=%d P%%4=%d : %d" % (k[0], k[1], k[2], sig[k]))
    return len(viol)

# --------------------------------------------------------------------------
# CLAIM (B): sharp RHS for v2(N1); the COMBINED direct-Delta inequality
# --------------------------------------------------------------------------
def claim_B(MMAX):
    print()
    print("=" * 72)
    print("CLAIM (B): v2(N1) bound + the combined Delta(b+1) inequality (m<=%d)" % MMAX)
    print("=" * 72)
    cand = {
        "v2(N1) >= v2(P+2)-1":            lambda vN1, P, a, b: vN1 >= v2(P+2) - 1,
        "v2(N1) >= v2(P+2)-1-v2(a-b+1)":  lambda vN1, P, a, b: vN1 >= v2(P+2) - 1 - v2(a - b + 1),
        "v2(N1) >= 1":                    lambda vN1, P, a, b: vN1 >= 1,
    }
    fails = {k: [] for k in cand}
    # what is the actual gap v2(N1) - v2(P+2)?  tabulate
    gapdist = {}
    n = 0
    for (a, b, m, P) in shapes(MMAX):
        n += 1
        vN1 = v2(N1_of(a, b))
        g = vN1 - v2(P + 2)
        gapdist[g] = gapdist.get(g, 0) + 1
        for k, f in cand.items():
            if not f(vN1, P, a, b):
                fails[k].append((a, b, m, vN1, v2(P+2), v2(a-b+1)))
    print("shapes scanned: %d" % n)
    for k in cand:
        print("  candidate  %-34s failures=%d  e.g. %s"
              % (k, len(fails[k]), fails[k][:3]))
    print("distribution of  v2(N1)-v2(P+2):", dict(sorted(gapdist.items())))

    # the inequality the proof actually needs: Delta(b+1) > -theta.
    # Use the direct Delta(b+1) formula and check > -theta, record equality slack.
    print("\n  direct Delta(b+1) > -theta  check + minimum slack:")
    minslack = 99; minshape = None; bad = 0; cnt = 0
    slackdist = {}
    for (a, b, m, P) in shapes(MMAX):
        j = b + 1
        if j > m:
            continue
        theta = 0 if a % 2 == 0 else 3
        D1 = ((b - 1) - 2*s2(b+1) + 2*K(2*(P+2), b+1)
              + 2*v2(N1_of(a, b)) - 2*v2(a - 1))
        slack = D1 + theta            # need slack > 0
        slackdist[slack] = slackdist.get(slack, 0) + 1
        cnt += 1
        if slack <= 0:
            bad += 1
        if slack < minslack:
            minslack = slack; minshape = (a, b, m)
    print("    indices checked=%d  violations(slack<=0)=%d  min slack=%d at (a,b,m)=%s"
          % (cnt, bad, minslack, minshape))
    print("    slack distribution:", dict(sorted(slackdist.items())))
    return bad

# --------------------------------------------------------------------------
# CLAIM (C): two-factor Lemma F2 premise, verified directly on the shapes
# --------------------------------------------------------------------------
def claim_C(MMAX):
    print()
    print("=" * 72)
    print("CLAIM (C): two-factor Lemma F2 premise on a-odd tops (m<=%d)" % MMAX)
    print("=" * 72)
    # a odd, b=2*beta.  The carry expansion builds the consecutive product
    #   U = prod_{t=1}^{Q}(P+t),  Q = beta+2  (= b//2 + 2).
    # The TWO deficits in Delta are 2 v2(a-1) and 2 v2(a-b+1); since a is odd,
    #   a-1   = 2(P+beta+1),   a-b+1 = 2(P+2),
    # so the OPERATIVE members of the consecutive product are the two integers
    #   f1 = P+2   (t=2)     and    f2 = P+beta+1   (t=beta+1).
    # The stray factor of 2 in each of a-1, a-b+1 is absorbed separately (gives
    # the +2 that the proof banks).  F2 PREMISE (the real one, matches
    # aodd_check2.py):   v2(U) - v2(f1) - v2(f2) >= 1   (Q>=6 guarantees it).
    def vprod(start, cnt):
        return sum(v2(x) for x in range(start, start + cnt))
    minslack = 99; minshape = None; bad = 0; cnt = 0
    slackdist = {}; smallQ = 0
    for (a, b, m, P) in shapes(MMAX):
        if a % 2 != 1:
            continue
        beta = b // 2
        Q = beta + 2
        if Q < 6:                 # threshold: box-interior a-odd has b>=10 => Q>=7
            smallQ += 1
            continue
        U = vprod(P + 1, Q)                # prod (P+1)...(P+Q)
        f1 = v2(P + 2)                     # designated member t=2
        f2 = v2(P + beta + 1)             # designated member t=beta+1
        slack = U - f1 - f2 - 1           # need slack >= 0  (premise >= 1)
        slackdist[slack] = slackdist.get(slack, 0) + 1
        cnt += 1
        if slack < 0:
            bad += 1
        if slack < minslack:
            minslack = slack; minshape = (a, b, m, P)
    print("a-odd shapes checked=%d (skipped Q<6: %d)" % (cnt, smallQ))
    print("F2 premise  v2(U)-v2(P+2)-v2(P+beta+1) >= 1   failures=%d" % bad)
    print("min surplus over the >=1 bound = %d  at (a,b,m,P)=%s" % (minslack, minshape))
    print("surplus distribution (premise value minus 1):", dict(sorted(slackdist.items())))
    return bad

# --------------------------------------------------------------------------
# cross-check the three DIRECT Delta formulas against Murnaghan-Nakayama
# --------------------------------------------------------------------------
def crosscheck_MN(MMAX):
    print()
    print("=" * 72)
    print("CROSS-CHECK: direct Delta(b+i) formulas vs Murnaghan-Nakayama (m<=%d)" % MMAX)
    print("=" * 72)
    bad = []; cnt = 0
    for (a, b, m, P) in shapes(MMAX):
        lam = (a, b, 3)
        v0 = val(lam, 0, m)
        if v0 is None:
            continue
        N2 = N2_of(a, b); N1 = N1_of(a, b)
        forms = {
            2: (b+2) - 2*s2(b+3) + 2*K(2*(P+1), b+3) + 2*v2(N2) - 2*v2(a-1) - 2*v2(P+2),
            1: (b-1) - 2*s2(b+1) + 2*K(2*(P+2), b+1) + 2*v2(N1) - 2*v2(a-1),
            3: (b+5) - 2*s2(b+5) + 2*K(2*P,     b+5)             - 2*v2(a-1) - 2*v2(P+2),
        }
        for i, Gi in forms.items():
            j = b + i
            if j > m:
                continue
            vj = val(lam, j, m)
            if vj is None:
                continue
            cnt += 1
            if vj - v0 != Gi:
                bad.append((lam, m, i, vj - v0, Gi))
    print("indices cross-checked=%d  mismatches=%d  %s"
          % (cnt, len(bad), bad[:8]))
    return len(bad)

if __name__ == "__main__":
    MMAX = 400
    vA = claim_A(MMAX)
    vB = claim_B(MMAX)
    vC = claim_C(MMAX)
    vX = crosscheck_MN(min(MMAX, 60))   # MN is expensive; cap at 60
    print()
    print("=" * 72)
    print("SUMMARY (m<=%d): A-bound-violations=%d  B-Delta-violations=%d  "
          "C-F2-premise-fails=%d  MN-mismatches=%d"
          % (MMAX, vA, vB, vC, vX))
    print("=" * 72)
