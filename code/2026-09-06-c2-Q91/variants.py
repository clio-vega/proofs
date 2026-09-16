"""Control variants, with a general weight callback so that DISTINCTNESS of the
controls can itself be tested.  (witnesses-must-be-checked-for-distinctness)"""
import sympy as sp
from bead import maya, from_maya, trim, bilinear, t


def R_variant(lam, e, weight, use_sign=True):
    """weight(sign, M, b, e) -> coefficient."""
    lam = trim(lam)
    lo = -(len(lam) + 2 * e + 6)
    M = maya(lam, lo)
    hi = (lam[0] if lam else 0) + 2 * e + 6
    out = {}
    for b in range(lo, hi + 1):
        if b not in M or b + e in M:
            continue
        res = bilinear(b + e, b, M)
        if res is None:
            continue
        sign, Mp = res
        c = weight(sign if use_sign else 1, M, b, e)
        mu = from_maya(Mp, lo)
        out[mu] = sp.expand(out.get(mu, 0) + c)
    return {k: v for k, v in out.items() if sp.expand(v) != 0}


def occ(M, a, bb):
    """#(M cap [a,bb])"""
    return sum(1 for j in range(a, bb + 1) if j in M)


# ---- the conjecture, and the variants ------------------------------------
def W_conj(s, M, b, e):    return s * (-t) ** occ(M, b + 1, b + e - 1)
def W_halfopen(s, M, b, e): return s * (-t) ** occ(M, b + 1, b + e)      # (b,b+e]
def W_closed(s, M, b, e):  return s * (-t) ** occ(M, b, b + e)           # [b,b+e]
def W_shift1(s, M, b, e):  return s * (-t) ** (occ(M, b + 1, b + e - 1) + 1)
def W_tsign(s, M, b, e):   return s * (t) ** occ(M, b + 1, b + e - 1)
def W_window(s, M, b, e):  return s * (-t) ** occ(M, b - e + 1, b - 1)    # wrong window
def W_holes(s, M, b, e):   return s * (-t) ** (e - 1 - occ(M, b + 1, b + e - 1))
def W_nosign_neg(s, M, b, e): return (-t) ** occ(M, b + 1, b + e - 1)     # sign dropped
def W_nosign_pos(s, M, b, e): return (t) ** occ(M, b + 1, b + e - 1)      # sign dropped, t^N
def W_nodress(s, M, b, e): return s                                       # dressing dropped

VARIANTS = [
    ("V0  conjecture:       sign * (-t)^#(M cap (b,b+e))", W_conj, True),
    ("V1  half-open (b,b+e]                             ", W_halfopen, True),
    ("V2  closed   [b,b+e]                              ", W_closed, True),
    ("V3  exponent N+1                                  ", W_shift1, True),
    ("V4  t^N instead of (-t)^N                         ", W_tsign, True),
    ("V5  wrong window (b-e,b)                          ", W_window, True),
    ("V6  holes: exponent (e-1)-N                       ", W_holes, True),
    ("V7  fermionic sign DROPPED, (-t)^N                ", W_nosign_neg, False),
    ("V8  fermionic sign DROPPED, t^N                   ", W_nosign_pos, False),
    ("V9  dressing DROPPED (bare bilinear = MN)         ", W_nodress, True),
]
