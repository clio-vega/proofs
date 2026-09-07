"""Engine A: R_e(t) from the COMBINATORIAL definition on partitions.
Adds a connected e-ribbon (border strip) with weight t^ht.  No abacus, no fermions.
"""
from itertools import product
import sympy as sp

t = sp.Symbol('t')

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,) + rest

def cells(lam):
    return {(i, j) for i, r in enumerate(lam) for j in range(r)}

def is_border_strip(sk):
    """sk: set of cells forming mu/lam.  Connected (edgewise) and no 2x2 block."""
    if not sk: return False
    for (i, j) in sk:
        if (i, j+1) in sk and (i+1, j) in sk and (i+1, j+1) in sk:
            return False
    # edge-connectivity
    start = next(iter(sk)); seen = {start}; stack = [start]
    while stack:
        (i, j) = stack.pop()
        for nb in ((i+1,j),(i-1,j),(i,j+1),(i,j-1)):
            if nb in sk and nb not in seen:
                seen.add(nb); stack.append(nb)
    return len(seen) == len(sk)

def contains(mu, lam):
    if len(mu) < len(lam): return False
    return all(mu[i] >= lam[i] for i in range(len(lam)))

def R_e(lam, e, N=None):
    """returns dict mu -> t^ht(mu/lam), over connected e-ribbons mu/lam."""
    out = {}
    n = sum(lam) + e
    for mu in partitions(n):
        if not contains(mu, lam): continue
        sk = cells(mu) - cells(lam)
        if len(sk) != e: continue
        if not is_border_strip(sk): continue
        rows = len({i for (i, j) in sk})
        out[mu] = t**(rows - 1)
    return out

def apply_op(vec, e):
    """vec: dict lam -> coeff.  Apply R_e."""
    out = {}
    for lam, c in vec.items():
        for mu, w in R_e(lam, e).items():
            out[mu] = sp.expand(out.get(mu, 0) + c*w)
    return {k: v for k, v in out.items() if v != 0}

def commutator(lam, e, f):
    """[R_e,R_f] s_lam as dict."""
    a = apply_op(apply_op({lam: sp.Integer(1)}, f), e)   # R_e R_f
    b = apply_op(apply_op({lam: sp.Integer(1)}, e), f)   # R_f R_e
    out = {}
    for k in set(a) | set(b):
        v = sp.expand(sp.simplify(a.get(k, 0) - b.get(k, 0)))
        if v != 0: out[k] = v
    return out

if __name__ == "__main__":
    print("R_1 s_() =", R_e((), 1))
    print("R_2 s_(1) =", R_e((1,), 2))
    print("[R_1,R_2] s_() =", commutator((), 1, 2))
    print("[R_1,R_2] s_(1) =", commutator((1,), 1, 2))
    print("[R_1,R_3] s_() =", commutator((), 1, 3))
    print("[R_2,R_3] s_() =", commutator((), 2, 3))
