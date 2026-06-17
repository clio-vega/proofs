"""
ayyer_kumari_valuebound.py  --  CODE session 2026-06-17, Job B

Probe (06-29 connection): does the Ayyer-Kumari {0,+-1,+-2} universal-character
value bound at roots of unity give a STRUCTURAL cause for "|J*| even" at d=4,
independent of my per-family closed-form arithmetic?  And does it pin (2,2)?

DOMAIN CHECK FIRST (this regime burned me once -- the fiber-splitting lever was DEAD).

Object:  G_lambda(i) = <s_lambda, (p2 + i e2)^m>  for lambda |- 2m  (exact Gaussian int).
  term_j := C(m,j) i^j (1+i)^j M_j,    G = sum_j term_j.
  v_pi(term_j) = j + 2 v2(C(m,j) M_j) = val(j)   (pi = 1+i, |pi|^2 = 2).
  mu := min_j val(j),   J* := {j : val(j)=mu},   |J*|.

Two DIFFERENT vanishing phenomena (the heart of the domain check):
  (a) LEADING-LAYER cancellation:  v_pi(G) > mu  <=>  the |J*| leading units sum
      to 0 mod pi  <=>  |J*| even.        [this is the "tie", the even-|J*| program]
  (b) FULL vanishing: G = 0 entirely  <=>  lambda = (2,2)  [trichotomy unique vanisher].

The probe asks whether a {0,+-1,+-2} bound on some normalised G forces (a).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import math
from job1_tie_census import partitions, M_vector, G_gaussian, vpi, v2, analyze

def v2_or(n):
    return v2(n)

def norm(re, im):
    return re * re + im * im

def pi_unit(re, im):
    """Divide Gaussian integer by pi=(1+i) until odd-norm; return (unit-part re,im, v_pi)."""
    if re == 0 and im == 0:
        return (0, 0, None)
    c = 0
    while (re + im) % 2 == 0:
        nre = (re + im) // 2
        nim = (im - re) // 2
        re, im = nre, nim
        c += 1
    return (re, im, c)

# ----------------------------------------------------------------------
def census(mmax):
    rows = []
    for m in range(1, mmax + 1):
        for lam in partitions(2 * m):
            r = analyze(lam, m)
            re, im = r['re'], r['im']
            ure, uim, vp = pi_unit(re, im)
            rows.append(dict(lam=lam, m=m, Jstar=r['Jstar'], absJ=len(r['Jstar']),
                             mu=r['mu'], re=re, im=im, vpi=vp,
                             norm=norm(re, im), ure=ure, uim=uim))
    return rows

# ======================================================================
if __name__ == "__main__":
    MMAX = 13
    print("#" * 78)
    print("# JOB B : Ayyer-Kumari {0,+-1,+-2} value bound as cause of |J*| even?")
    print("#" * 78)
    rows = census(MMAX)
    print("\nshapes computed: %d  (m<=%d)" % (len(rows), MMAX))

    # ---- (0) domain check: leading-layer cancel vs |J*| parity --------
    print("\n" + "=" * 78)
    print("(0) DOMAIN CHECK: is  v_pi(G) > mu  <=>  |J*| even ?")
    print("=" * 78)
    bad = 0
    for r in rows:
        cancels = (r['vpi'] is None) or (r['vpi'] > r['mu'])   # G=0 counts as cancel
        even = (r['absJ'] % 2 == 0)
        if cancels != even:
            bad += 1
            if bad <= 8:
                print("  MISMATCH lam=%s m=%d |J*|=%d mu=%d vpi=%s"
                      % (r['lam'], r['m'], r['absJ'], r['mu'], r['vpi']))
    print("  mismatches (leading-cancel xor |J*| even) = %d / %d" % (bad, len(rows)))

    # ---- (1) is G_lambda(i) itself bounded?  (value bound on the RAW object) ----
    print("\n" + "=" * 78)
    print("(1) IS G_lambda(i) BOUNDED?  (the literal AK value-bound question)")
    print("=" * 78)
    norms = sorted(set(r['norm'] for r in rows))
    big = max(rows, key=lambda r: r['norm'])
    print("  distinct |G|^2 values (first 20): %s ..." % norms[:20])
    print("  MAX |G|^2 = %d at lam=%s (m=%d)  -> G grows; NOT a bounded set"
          % (big['norm'], big['lam'], big['m']))
    # growth with m
    print("  max |G|^2 per m:")
    for m in range(1, MMAX + 1):
        mx = max((r['norm'] for r in rows if r['m'] == m), default=0)
        print("     m=%2d : max|G|^2=%d" % (m, mx))

    # ---- (2) normalised objects: do ANY land in {0,+-1,+-2}-type small set? ----
    print("\n" + "=" * 78)
    print("(2) NORMALISED G: pi-unit part  G/pi^{v_pi}  (odd-norm Gaussian int)")
    print("=" * 78)
    unitnorms = sorted(set(r['ure']**2 + r['uim']**2 for r in rows if r['vpi'] is not None))
    print("  distinct |G/pi^{vpi}|^2 (first 25): %s" % unitnorms[:25])
    bigu = max((r for r in rows if r['vpi'] is not None),
               key=lambda r: r['ure']**2 + r['uim']**2)
    print("  MAX unit-part |.|^2 = %d at lam=%s (m=%d) -> unit part ALSO unbounded"
          % (bigu['ure']**2 + bigu['uim']**2, bigu['lam'], bigu['m']))

    # ---- (3) THE MAKE-OR-BREAK: does any value-bound DECIDE |J*| parity? ----
    print("\n" + "=" * 78)
    print("(3) MAKE-OR-BREAK: does a value bound on G force |J*| parity?")
    print("=" * 78)
    # Test: among shapes, is |J*| parity a function of (G mod small)?  i.e. is there
    # a value-level invariant of G that determines |J*| even/odd WITHOUT M_j data?
    # The leading-layer cancellation is BY DEFINITION |J*| mod 2; the question is
    # whether AK gives this for free.  AK bounds the VALUE; |J*| is a SUPPORT count.
    # Concrete falsifier: find two shapes with the SAME normalised value G/pi^{vpi}
    # (same unit) but DIFFERENT |J*| parity -> value can't decide parity.
    from collections import defaultdict
    byunit = defaultdict(list)
    for r in rows:
        if r['vpi'] is not None:
            byunit[(abs(r['ure']), abs(r['uim']))].append(r)
    agnostic = []
    for key, lst in byunit.items():
        pars = set(x['absJ'] % 2 for x in lst)
        if len(pars) > 1:
            agnostic.append((key, lst))
    print("  unit-value classes carrying BOTH |J*| parities: %d" % len(agnostic))
    if agnostic:
        key, lst = agnostic[0]
        print("  e.g. |unit|=(%d,%d) realised by:" % key)
        for x in lst[:6]:
            print("       lam=%s m=%d |J*|=%d (%s)" %
                  (x['lam'], x['m'], x['absJ'], "even" if x['absJ'] % 2 == 0 else "odd"))

    # ---- (4) pin (2,2): is G=0 ONLY at (2,2)? --------------------------
    print("\n" + "=" * 78)
    print("(4) PIN (2,2): full vanishers G_lambda(i)=0")
    print("=" * 78)
    vanish = [r['lam'] for r in rows if r['re'] == 0 and r['im'] == 0]
    print("  G=0 for: %s" % vanish)
    print("  (trichotomy predicts ONLY (2,2))")

    # ---- (5) |J*| distribution: powers of 2 ----------------------------
    print("\n" + "=" * 78)
    print("(5) |J*| distribution (the even-|J*| program: |J*|=2^|S|)")
    print("=" * 78)
    from collections import Counter
    cnt = Counter(r['absJ'] for r in rows)
    print("  |J*| -> count: %s" % dict(sorted(cnt.items())))
    print("  even-|J*| (ties): %d ; odd-|J*|: %d"
          % (sum(v for k, v in cnt.items() if k % 2 == 0),
             sum(v for k, v in cnt.items() if k % 2 == 1)))
