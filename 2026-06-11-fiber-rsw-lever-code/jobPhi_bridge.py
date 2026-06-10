"""
GAP (i) coarse->sharp, done right: the parity-restricted bridge.

Finding (jobPhi_lift): the coarse box B* (support of Phi/2^e mod 2) is even on ~everything,
INCLUDING unique-min |J*|=1 shapes, because the (1+z)-lift discards the parity tilt
(val(j+1)-val(j) is ALWAYS odd, so J* sits in a single parity class).  So "B* even" does NOT
imply "|J*| even".  The honest bridge must intersect B* with the parity class of J*.

Hypotheses tested (all lambda |- 2m, m<=MMAX):
  (H1) J* sits in a single parity class p (no two J* indices of opposite parity).
  (H2) J* is an affine 2-adic box  j0 + submask(S),  |J*| = 2^|S| in {1,2,4}.
  (H3) BRIDGE:  J* = B* ∩ {k ≡ p (mod 2)}  for p = parity(min J*).
  (H4) If H3 holds: B* = j0 + submask(g); intersecting with a fixed lowest bit gives a sub-box.
       Then |J*| even  <=>  g has a SECOND set bit beyond bit 0, i.e. the sharp box survives.
       Test the refined statement:  |J*| = 2^{(#set bits of g in positions >=1)}  when the
       parity slice is nonempty -- i.e. the bit-0 of g is "spent" by the parity restriction.
"""
import sys, math, csv, os
from collections import Counter
from job1_tie_census import partitions, M_vector, v2, chi_b_vector
from jobPhi_lift import coarse_box, sharp_Jstar

RESULTS = "/home/clio/projects/code/results"

def popcount(x): return bin(x).count('1')

def box_struct(J):
    """If J is an affine 2-adic box j0+submask(S), return (j0, S_exponents, mask). Else None."""
    J = sorted(J)
    j0 = J[0]
    diffs = [j - j0 for j in J]
    mask = 0
    for d in diffs:
        mask |= d
    # submask set of `mask`:
    subs = set()
    s = mask
    while True:
        subs.add(s)
        if s == 0: break
        s = (s-1) & mask
    if set(diffs) == subs:
        exps = [i for i in range(mask.bit_length()) if (mask>>i)&1]
        return (j0, exps, mask)
    return None

def main(MMAX):
    H1 = H1tot = 0
    H2 = H2tot = 0
    H3 = H3tot = 0
    H4 = H4tot = 0
    Jsize_dist = Counter()
    h3_fail = []
    h4_fail = []
    g0bit_vs_parity = Counter()
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            cb = coarse_box(lam, m)
            if cb is None or cb['jg'] is None:
                continue
            j0c, g = cb['jg']
            B = set(cb['supp'])
            J, mu = sharp_Jstar(lam, m)
            Jset = set(J)
            Jsize_dist[len(J)] += 1

            # H1: single parity class
            H1tot += 1
            parities = set(k % 2 for k in J)
            single_parity = (len(parities) == 1)
            if single_parity:
                H1 += 1
            p = min(J) % 2

            # H2: J* is an affine 2-adic box
            H2tot += 1
            bs = box_struct(J)
            if bs is not None and len(J) in (1, 2, 4):
                H2 += 1

            # H3: J* = B* ∩ parity-class(p)
            H3tot += 1
            B_par = set(k for k in B if k % 2 == p)
            if B_par == Jset:
                H3 += 1
            else:
                if len(h3_fail) < 20:
                    h3_fail.append((m, tuple(lam), sorted(J), sorted(B), p))

            # H4: |J*| = 2^{popcount of g restricted to bits >=1}  when bit0 of g matches parity shift
            # the coarse box mask is g (since j0c + submask(g)); parity slice fixes bit 0.
            H4tot += 1
            # bits of g above position 0:
            g_hi = g & ~1
            pred = 2 ** popcount(g_hi)
            # but only the slice with the correct lowest bit is nonempty; J* size should be pred
            if len(J) == pred:
                H4 += 1
            else:
                if len(h4_fail) < 20:
                    h4_fail.append((m, tuple(lam), len(J), g, pred))
            g0bit_vs_parity[(g & 1, p, len(J))] += 1

    print(f"=== coarse->sharp parity bridge  (m<={MMAX}) ===")
    print(f"  |J*| distribution: {dict(sorted(Jsize_dist.items()))}")
    print(f"  H1  J* single parity class:                 {H1}/{H1tot}")
    print(f"  H2  J* affine 2-adic box, |J*| in {{1,2,4}}:  {H2}/{H2tot}")
    print(f"  H3  J* = B* ∩ parity-class:                  {H3}/{H3tot}")
    if h3_fail:
        print(f"      H3 failures (first {len(h3_fail)}):")
        for r in h3_fail:
            print(f"        m={r[0]} lam={r[1]} J*={r[2]} B*={r[3]} p={r[4]}")
    print(f"  H4  |J*| = 2^popcount(g & ~1):               {H4}/{H4tot}")
    if h4_fail:
        print(f"      H4 failures (first {len(h4_fail)}):")
        for r in h4_fail:
            print(f"        m={r[0]} lam={r[1]} |J*|={r[2]} g={r[3]} pred={r[4]}")

if __name__ == "__main__":
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    main(MMAX)
