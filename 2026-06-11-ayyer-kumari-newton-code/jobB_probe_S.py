"""
JOB B probe -- what controls the box generator set S (offsets in {2,4}) ?

H1 forces offsets to be EVEN (single parity).  Observed generators: 2 (=2^1),
4 (=2^2), or both.  So S subset {1,2} (as exponents) at m<=8.  We hunt for a
clean determinant of S, prioritising Tool A:

   D_j := 2^j M_j = sum_{b<=j} C(j,b)(-1)^b chi_b   (signed binomial transform),
   v2(M_j) = v2(D_j) - j,   val(j) = 2 v2(C(m,j)) + 2 v2(D_j) - j.

Candidate determinants for which a in {1,2} is a generator:
  * binary digits / v2 of m,
  * the 2-adic valuation profile of D_j near j0,
  * whether the toggle j0 -> j0+2^a stays valuation-minimal (definitional).

Output: for each tie, j0, S, m, v2(m), and the local val profile, to eyeball the rule.
"""
import math
from collections import Counter, defaultdict

from job1_tie_census import partitions, v2, M_vector
from job_box_detail import is_affine_box


def data(MMAX=8):
    rows = []
    for m in range(1, MMAX + 1):
        for lam in partitions(2 * m):
            Ms = M_vector(lam, m)
            val = {j: j + 2 * v2(math.comb(m, j) * Ms[j]) for j in range(m + 1) if Ms[j] != 0}
            mu = min(val.values())
            J = sorted(k for k in val if val[k] == mu)
            if len(J) < 2:
                continue
            ok, gens = is_affine_box(J)
            S = tuple(int(g).bit_length() - 1 for g in gens) if ok else None  # exponents
            rows.append(dict(m=m, lam=lam, j0=J[0], S=S, J=J, val=val, Ms=Ms))
    return rows


def main():
    rows = data(8)
    print(f"ties m<=8: {len(rows)}")

    # 1) S vs (m parity, v2(m))
    print("\n[1] generator exponent multiset S vs (m mod 4, v2(m)):")
    by = defaultdict(Counter)
    for r in rows:
        by[(r['m'] % 4, v2(r['m']))][r['S']] += 1
    for key in sorted(by):
        print(f"   m%4={key[0]} v2(m)={key[1]:>2}: {dict(by[key])}")

    # 2) Is membership of exponent a in S equivalent to: val(j0+2^a)==mu  (definitional
    #    -- yes) but ALSO val(j0+2^{a'})>mu for the OTHER exponent?  i.e. independence.
    #    Test: does the pair (val(j0+2)-mu, val(j0+4)-mu) determine S?  (trivially yes)
    #    More useful: is there a *gap law* val(j0+2^a)-mu as a clean function?
    print("\n[2] offset deltas delta_a = val(j0+2^a) - mu  for a=1,2 (when in range):")
    delta_dist = Counter()
    for r in rows:
        j0, val, mu = r['j0'], r['val'], min(r['val'].values())
        d1 = val.get(j0 + 2, None)
        d2 = val.get(j0 + 4, None)
        d1 = (d1 - mu) if d1 is not None else None
        d2 = (d2 - mu) if d2 is not None else None
        delta_dist[(d1, d2, r['S'])] += 1
    for k, c in sorted(delta_dist.items(), key=lambda kv: -kv[1]):
        print(f"   (delta@+2, delta@+4)={k[:2]!s:16s} S={k[2]!s:8s}: {c}")

    # 3) the val curve is piecewise: val(j)-val(j0) for j near j0.  Tabulate the
    #    shape of the local Newton vertices that produce a width-4 box.
    print("\n[3] width-4 boxes (S=(1,2)): local val profile j0..j0+4")
    for r in rows:
        if r['S'] == (1, 2):
            j0, val, mu = r['j0'], r['val'], min(r['val'].values())
            prof = [(j, val.get(j, None) and val[j] - mu) for j in range(j0, j0 + 5)]
            print(f"   lam={r['lam']!s:18s} m={r['m']} j0={j0} J*={r['J']} "
                  f"prof(j-mu)={[p[1] for p in prof]}")

    # 4) PROVE-HANDOFF candidate: does the generator a appear  <=>  v2(M_{j0+2^a})
    #    equals a specific value tied to v2(M_{j0})?  Tabulate v2(M_{j0+2^a})-v2(M_{j0}).
    print("\n[4] v2(M_{j0+2^a}) - v2(M_{j0}) for generators a present:")
    jump = Counter()
    for r in rows:
        if r['S'] is None:
            continue
        Ms, j0 = r['Ms'], r['j0']
        for a in r['S']:
            off = 1 << a
            if j0 + off < len(Ms) and Ms[j0] and Ms[j0 + off]:
                jump[(a, v2(Ms[j0 + off]) - v2(Ms[j0]))] += 1
    for k, c in sorted(jump.items()):
        a, dv = k
        # val-constancy: 0 = 2^a + 2(v2 C(m,j0+2^a)-v2 C(m,j0)) + 2 dv  =>
        # dv = -2^{a-1} - (v2C jump).  Record the binomial part too.
        print(f"   gen a={a} (offset {1<<a}): v2(M jump)={dv:>3} : {c} pairs")


if __name__ == "__main__":
    main()
