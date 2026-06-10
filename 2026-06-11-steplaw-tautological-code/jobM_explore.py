"""
Explore the M-half step law via the finite-difference reformulation.

  D_j := sum_{b=0}^j C(j,b)(-1)^b chi_b     (so M_j = D_j / 2^j, integer)
  v2(M_j) = v2(D_j) - j.
  chi_b = chi^lambda(2^b 1^{2m-2b}),  all chi_b == f^lambda (mod 2).

Claim (M*) on pure-M toggle pairs {j, j+8} of J*:  v2(M_{j+8})-v2(M_j) = -4,
  i.e. v2(D_{j+8}) - v2(D_j) = +4.

Goal: find the MECHANISM. Print chi_b 2-adic data, D_j data, for tie shapes.
"""
import math, itertools
from job1_tie_census import chi_b_vector, M_vector, analyze, partitions, v2

def D_vector(lam, m):
    chi = chi_b_vector(lam, m)
    D = []
    for j in range(m+1):
        s = sum(math.comb(j,b)*((-1)**b)*chi[b] for b in range(j+1))
        D.append(s)
    return D, chi

def pure_M_pairs(m_max=12):
    """Yield (lam,m,j,a) for every fpf toggle pair {j,j+2^a} of J* with a=3 (pure-M)."""
    out = []
    for m in range(1, m_max+1):
        for lam in partitions(2*m):
            r = analyze(lam, m)
            J = r['Jstar']
            if len(J) < 2:
                continue
            Jset = set(J)
            for j in J:
                for a in (1,2,3):
                    if (j + (1<<a)) in Jset and (j ^ (1<<a)) == j + (1<<a):
                        # {j, j+2^a} a toggle pair (j has bit a = 0)
                        out.append((tuple(lam), m, j, a))
    return out

if __name__ == "__main__":
    pairs = pure_M_pairs(12)
    # focus a=3
    a3 = [p for p in pairs if p[3]==3]
    print(f"total toggle pairs found: {len(pairs)};  a=3 pairs: {len(a3)}")
    # verify the M* claim and dump mechanism for first several a=3 pairs
    print("\n=== a=3 (pure-M) pairs: v2(D) jump check ===")
    ok = 0
    for (lam,m,j,a) in a3:
        D,chi = D_vector(lam,m)
        Ms = M_vector(lam,m)
        vj, vj8 = v2(Ms[j]), v2(Ms[j+8])
        dM = vj8 - vj
        vDj, vDj8 = v2(D[j]), v2(D[j+8])
        ok += (dM == -4)
        if (lam,m,j) in [a3[k][:3] for k in range(min(8,len(a3)))]:
            print(f"lam={lam} m={m} j={j}->{j+8}: v2(M_j)={vj} v2(M_j+8)={vj8} dM={dM} | v2(D_j)={vDj} v2(D_j+8)={vDj8} dD={vDj8-vDj}")
    print(f"\nM* holds (dM=-4): {ok}/{len(a3)} a=3 pairs")
