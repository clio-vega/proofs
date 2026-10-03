"""y-coordinates, width map, nest.  All identities from A-half-width-is-l1-distance.

   a_i   = lam_{i-1}+1,  y_i = nu_i - a_i,  g_i = lam_i - lam_{i-1} - 1
   w_i   = g_i + 1 - (y_i)_+ - (-y_{i+1})_+                              (eq:w)
   Lam   = (G - ||y||_1)/2,   k(y) = sum_i (-y_i)_+ = (||y||_1 - sigma)/2
   nest  Sigma^{<=j} = { y in Y : k(y) <= j }
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen
from gen import lam_prev, nu_next, conv_intervals, is_pf2
from itertools import product

pos = lambda t: t if t > 0 else 0

def gaps(lam, n, m):
    return [lam[i] - lam_prev(lam, i, n, m) - 1 for i in range(m)]

def avec(lam, n, m):
    return [lam_prev(lam, i, n, m) + 1 for i in range(m)]

def to_y(nu, lam, n, m):
    a = avec(lam, n, m)
    return tuple(nu[i] - a[i] for i in range(m))

def wmap(y, g):
    """eq:w -- widths from y and the gap vector g (cyclic in i)."""
    m = len(y)
    return tuple(g[i] + 1 - pos(y[i]) - pos(-y[(i+1) % m]) for i in range(m))

def defect(y):
    return sum(pos(-t) for t in y)

def trap(w):
    """convolution of 1_[0,w_i-1], CENTRED: returns coeff list (symmetric)."""
    return conv_intervals(list(w))

def centred_sum(ws):
    """sum of Trap_w for w in ws, all centred at a common centre.
       half-width Lam = (sum w_i - m)/2 ; returns coeff list centred."""
    if not ws: return []
    m = len(ws[0])
    half = [ (sum(w) - m) for w in ws ]        # = 2*Lam, must all have same parity
    H = max(half)
    L = H // 2 if H % 2 == 0 else None
    # work on a grid of half-integer offsets: index 2*x
    acc = {}
    for w in ws:
        co = trap(w)
        if not co: continue
        # co has length 2*Lam+1 ... actually len(co)-1 = sum(w_i)-m = 2*Lam
        tw = len(co) - 1          # = 2*Lam
        for idx, c in enumerate(co):
            key = 2*idx - tw      # in units of 1/2, centred at 0
            acc[key] = acc.get(key, 0) + c
    if not acc: return []
    # GUARD: all summands must sit on ONE parity class of the half-integer grid,
    # i.e. all sum_i w_i must share a parity.  Otherwise the summands are not
    # concentric on a common integer lattice and the question is ill-posed.
    # (Without this guard range(lo,hi+1,2) silently DROPPED one parity class.)
    if len({k % 2 for k in acc}) != 1:
        return None
    lo, hi = min(acc), max(acc)
    return [acc.get(k, 0) for k in range(lo, hi+1, 2)]
