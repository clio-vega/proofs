"""
g>=1 test for the e2-mod-2 layer of the exact lift Phi.

Phi_lambda(z) = <s_lambda, (p1^2 + z p2)^m>.  With w = 1+z:
  hatPhi(w) = <s_lambda, (w p2 + 2 e2)^m> = sum_r C(m,r) 2^r R_r w^{m-r},
  R_r = <s_lambda, p2^{m-r} e2^r>.
Coefficient of w^{m-r} is c_{m-r} = C(m,r) 2^r R_r.
  L_r := v2(c_{m-r}) = v2(C(m,r)) + r + v2(R_r).
  e = min_r L_r,   r* = max{r : L_r = e},   g = m - r*.
g>=1  <=>  L_m > e  <=>  exists r<m with L_r < L_m = m + v2(R_m).

"Tie shape": chi^lambda(2^m) = R_0 even.

We compute, for all lambda |- 2m, the L_r vector, e, g, and flag any TIE with g=0.
Also record diagnostic quantities.
"""
import sys
from sympy import Integer
import symfunc as sf

def v2(n):
    n = int(n)
    if n == 0:
        return None  # +infinity
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k

def binom(n, k):
    from math import comb
    return comb(n, k)

def analyze(lam, m):
    slam = sf.schur_p(lam)
    p2 = {(2,): Integer(1)}
    e2 = sf.e_n(2)
    # R_r = <s_lam, p2^{m-r} e2^r>
    R = []
    for r in range(m+1):
        f = sf.mul(sf.power(p2, m-r), sf.power(e2, r))
        R.append(int(sf.inner(slam, f)))
    # L_r
    L = []
    for r in range(m+1):
        vr = v2(R[r])
        if vr is None:
            L.append(None)  # R_r = 0  -> term absent
        else:
            L.append(v2(binom(m, r)) + r + vr)
    # e = min over r with L_r not None
    finite = [(r, L[r]) for r in range(m+1) if L[r] is not None]
    e = min(l for _, l in finite)
    rstar = max(r for r, l in finite if l == e)
    g = m - rstar
    R0 = R[0]  # = chi^lambda(2^m)
    is_tie = (R0 % 2 == 0)
    return dict(lam=lam, m=m, R=R, L=L, e=e, rstar=rstar, g=g, R0=R0, is_tie=is_tie)

def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    bad = []
    n_tie = 0
    e_dist = {}
    for m in range(1, mmax+1):
        for lam in sf.partitions(2*m):
            d = analyze(lam, m)
            if d['is_tie']:
                n_tie += 1
                e_dist[d['e']] = e_dist.get(d['e'], 0) + 1
                if d['g'] < 1:
                    bad.append(d)
    print(f"ties (m<= {mmax}): {n_tie}")
    print(f"e distribution among ties: {dict(sorted(e_dist.items()))}")
    print(f"TIES WITH g=0 (would refute g>=1): {len(bad)}")
    for d in bad[:20]:
        print("  BAD", d['lam'], "m=", d['m'], "e=", d['e'], "g=", d['g'])

if __name__ == "__main__":
    main()
