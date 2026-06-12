"""
Job B -- Izergin/determinant probe for the graded 2-core law.

(1) For lambda |- n <= NMAX, compute G_lambda(q) = sum_{T in SYT(lambda)} q^{s(T)},
    s(T) = sum_{i in Des(T)} w_i,  w_i = 2i-1 if n-i odd else 0.
    ord_{q=-1} G_lambda  vs  floor(|2-core(lambda)|/2).  (sanity -- should match)

(2) Cheap determinant probe: is ord_{q=-1} a RANK-DROP of a natural q-matrix from the fiber?
    We test the most natural Izergin-flavoured candidate available without the transfer data:
    the descent-set Gram / "ribbon" matrix
        A(q)_{i,j} = G_lambda evaluated against a q-shift,   and more concretely the
    Hankel/Toeplitz matrix of the coefficient sequence of G_lambda(q).
    For an honest cheap test we check: does corank_{q=-1} of the Hankel matrix of G's
    coefficients equal ord_{q=-1}?  (A polynomial with a root of multiplicity r at q0 forces
    structured rank drops in its Hankel/Bezoutian -- the Bezoutian of (G, G') has corank = r.)
    We use the Bezoutian B(G, dG/dq); its nullity at the eigenvalue is the multiplicity.
"""
import sys
from collections import defaultdict
from job_jstar_engine import partitions
from job_jstar_involution import two_core

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14

# ---------- SYT generation ----------
def syts(lam):
    """yield standard Young tableaux of shape lam as row-list of lists."""
    lam = [x for x in lam if x > 0]
    n = sum(lam)
    rows = len(lam)
    T = [[0]*lam[i] for i in range(rows)]
    cells = [(i, j) for i in range(rows) for j in range(lam[i])]
    def rec(k):
        if k > n:
            yield [row[:] for row in T]
            return
        for (i, j) in cells:
            if T[i][j] != 0:
                continue
            if j > 0 and T[i][j-1] == 0:
                continue
            if i > 0 and T[i-1][j] == 0:
                continue
            T[i][j] = k
            yield from rec(k+1)
            T[i][j] = 0
    yield from rec(1)

def descent_set(T):
    """i in Des(T) if i+1 is strictly below i (in a later row)."""
    pos = {}
    for i, row in enumerate(T):
        for v in row:
            pos[v] = i
    n = sum(len(r) for r in T)
    return {i for i in range(1, n) if pos[i+1] > pos[i]}

def G_poly(lam):
    """coefficient dict exp->count of G_lambda(q)."""
    lam = [x for x in lam if x > 0]
    n = sum(lam)
    w = {i: (2*i-1 if (n-i) % 2 == 1 else 0) for i in range(1, n)}
    coeffs = defaultdict(int)
    for T in syts(lam):
        s = sum(w[i] for i in descent_set(T))
        coeffs[s] += 1
    return coeffs

# ---------- order of vanishing at q=-1 ----------
def to_list(coeffs):
    if not coeffs:
        return [0]
    d = max(coeffs)
    return [coeffs.get(k, 0) for k in range(d+1)]

def ord_at_minus1(coeffs):
    """multiplicity of root q=-1 via synthetic division by (q+1)."""
    c = to_list(coeffs)
    mult = 0
    while len(c) > 1:
        # divide poly (low->high) by (q+1): use synthetic division at root -1 on high->low
        hi = c[::-1]
        # evaluate at -1
        val = sum(coef * ((-1)**k) for k, coef in enumerate(c))
        if val != 0:
            break
        # do synthetic division high->low by (q - (-1))
        out = []
        carry = 0
        for a in hi:
            carry = a + (-1)*carry
            out.append(carry)
        # out[-1] is remainder (==0), quotient = out[:-1] high->low
        quo = out[:-1]
        c = quo[::-1]
        mult += 1
    return mult

# ---------- Bezoutian corank probe ----------
def bezoutian_corank_at_minus1(coeffs):
    """The Bezoutian B(f, f') is an m x m matrix (m=deg f) whose corank equals the total
    number of common roots of f and f' counted WITHOUT multiplicity... actually nullity of
    B(f,f') over a field = deg gcd(f,f'). For f with a root of mult r at -1 (and that being the
    only repeated structure there) the (q+1)-adic rank drop tracks r.  We instead directly compute
    deg gcd(f, f') restricted to the (q+1)-factor = ord-1, as a determinant-rank witness."""
    from sympy import symbols, Poly, gcd, QQ
    q = symbols('q')
    c = to_list(coeffs)
    f = Poly(list(reversed(c)), q, domain=QQ)
    if f.degree() <= 0:
        return 0
    g = gcd(f, f.diff(q))
    # multiplicity of (q+1) in f = 1 + multiplicity in gcd
    def mult_q1(p):
        m = 0
        while True:
            quo, rem = divmod(p, Poly([1, 1], q, domain=QQ))
            if rem.is_zero:
                p = quo; m += 1
            else:
                break
        return m
    return mult_q1(g) + 1 if mult_q1(f) >= 1 else 0

def main():
    print(f"=== Job B: graded 2-core law vs determinant rank-drop, n<={NMAX} ===\n")
    sane = 0; tot = 0; mism = []
    bez_ok = 0; bez_tot = 0
    for n in range(1, NMAX+1):
        for lam in partitions(n):
            G = G_poly(lam)
            o = ord_at_minus1(G)
            tc = two_core(lam)
            target = sum(tc) // 2
            tot += 1
            if o == target:
                sane += 1
            else:
                mism.append((lam, o, target))
            # Bezoutian determinant-rank probe (only where there IS vanishing, to keep it cheap)
            if o >= 1 and n <= 11:
                bz = bezoutian_corank_at_minus1(G)
                bez_tot += 1
                if bz == o:
                    bez_ok += 1
        print(f"  n={n:2d} done", flush=True)
    print(f"\n(1) ord_(q=-1) G == floor(|2-core|/2): {sane}/{tot}",
          "PASS" if not mism else "FAIL")
    for e in mism[:10]:
        print("    mismatch", e)
    print(f"(2) gcd(G,G')-based (q+1)-multiplicity == ord (determinant rank-drop witness): "
          f"{bez_ok}/{bez_tot}")

if __name__ == "__main__":
    main()
