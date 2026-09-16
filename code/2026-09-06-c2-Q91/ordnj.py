"""C2: ord(n_j) = infinity.  Verify the proof's explicit witness numerically.

Claim proved on paper:  (ad_{M_{p1}}^n n_j)_{mu,lam} = sum_{paths} Delta^n nu(path),
and for a single-ROW (or single-COLUMN) path the path is unique, giving +-1.
Here we recompute ad^n n_j honestly as an operator on Lambda (Schur basis,
Pieri for M_{p1}) and report its norm, for n up to 6.
"""
import itertools
from bead import maya, trim

def nu(lam, j):
    lo = -(len(trim(lam)) + 4 + abs(j))
    return 1 if j in maya(lam, lo) else 0

def addbox(lam):
    lam = list(trim(lam)) + [0]
    out = []
    for i in range(len(lam)):
        if i == 0 or lam[i] < lam[i - 1]:
            m = lam[:]; m[i] += 1
            out.append(trim(tuple(m)))
    return out

def parts(n):
    if n == 0: return [()]
    out = []
    def rec(rem, cap, cur):
        if rem == 0: out.append(tuple(cur)); return
        for p in range(min(rem, cap), 0, -1): rec(rem - p, p, cur + [p])
    rec(n, n, [])
    return out

from math import comb
def ad_n(j, n, lam):
    """(ad_{M_p1}^n n_j) applied to s_lam -> dict mu:coeff."""
    # sum over paths lam=l_0 -> ... -> l_n, weight Delta^n nu
    out = {}
    frontier = {(lam,): 1}
    for _ in range(n):
        nf = {}
        for path, w in frontier.items():
            for m in addbox(path[-1]):
                nf[path + (m,)] = nf.get(path + (m,), 0) + w
        frontier = nf
    for path, w in frontier.items():
        d = sum((-1) ** (n - i) * comb(n, i) * nu(path[i], j) for i in range(n + 1))
        if d:
            out[path[-1]] = out.get(path[-1], 0) + w * d
    return {k: v for k, v in out.items() if v != 0}

print("nonvanishing of ad_{M_p1}^n (n_j)  --  entries at the predicted witness")
for j in (-3, -2, -1, 0, 1, 2):
    row = []
    for n in range(1, 7):
        # scan all lam with |lam| <= 6 for a nonzero image
        hit = None
        for k in range(0, 7):
            for lam in parts(k):
                r = ad_n(j, n, lam)
                if r: hit = (lam, r); break
            if hit: break
        row.append("nonzero" if hit else "ZERO")
    print(f"  j={j:3d}:  " + "  ".join(f"n={n}:{v}" for n, v in zip(range(1, 7), row)))

print()
print("the PROOF's explicit witness (single row / single column, unique path):")
for j in (0, 1, 2):
    for n in (1, 2, 3, 4, 5):
        lam = (j + 1,); mu = (j + 1 + n,)
        r = ad_n(j, n, lam)
        print(f"  j={j} n={n}:  lam={lam} -> coeff at mu={mu} is {r.get(mu, 0)}"
              f"   (predicted {(-1)**n})")
    break
for j in (-1,):
    for n in (1, 2, 3, 4, 5):
        r = ad_n(j, n, ())
        print(f"  j={j} n={n}:  lam=() -> coeff at mu={(n,)} is {r.get((n,), 0)}"
              f"   (predicted {(-1)**n})")
for j in (-2, -3):
    for n in (1, 2, 3, 4):
        lam = tuple([1] * (-j)); mu = tuple([1] * (-j + n))
        r = ad_n(j, n, lam)
        print(f"  j={j} n={n}:  lam={lam} -> coeff at mu={mu} is {r.get(mu, 0)}"
              f"   (predicted {-(-1)**n})")
