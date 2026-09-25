"""Korff-Palazzo side (arXiv:1804.05647), independently implemented from their definitions.

Dictionary (proved/calibrated below):
  KP cylinder C_{k,n-k}  <->  AKO C_{x,y} with x = n-k, y = k.
  boxed partition lam-bar (k parts, each <= n-k), lam_{i+k} = lam_i - (n-k).
  lam/d/mu  =  boxes in KP row i, columns mu_i+1 .. lam_{i-d}+d.
  drawn row r = k-1-... : we set KP row i = k - r so that lo,hi are weakly increasing.
KP height of a cylindric ribbon = number of rows its (extended) ribbon occupies;
  each circular ribbon contributes k, consecutive pieces SHARE a row, so
  ht = rows(tail) + l*(k-1)   (non-pure),   ht = m*(k-1) + 1  (m circular ribbons, no tail).
Per-ribbon weight in chi is (-1)^{ht-1}.
"""
from cyl import *
from ako import decompose, is_loop_ribbon, intermediates, region, size_of
from weight import geom

def boxed(lam, k, nk):
    "extend a boxed partition to all of Z:  lam_{i+k} = lam_i - nk"
    lam = list(lam) + [0]*(k-len(lam))
    assert all(nk >= lam[0] for _ in [0]) and all(lam[i] >= lam[i+1] for i in range(k-1)) and lam[-1] >= 0
    def f(i):
        q, r = divmod(i-1, k)
        return lam[r] - q*nk
    return f

def build(n, k, lam, mu, d):
    "the cylindric skew shape lam/d/mu as a Cyl on C_{x=n-k, y=k}"
    nk = n-k
    L, M = boxed(lam, k, nk), boxed(mu, k, nk)
    lo = [M(k-r)+1 for r in range(k)]
    hi = [L(k-r-d)+d for r in range(k)]
    return Cyl(nk, k, lo, hi)

def kp_weight(R):
    "(-1)^{ht_KP - 1} for a cylindric ribbon R; 0 if R is not one"
    dd = decompose(R)
    if dd is None: return 0
    l, F = dd
    y = R.y
    if is_loop_ribbon(F):
        m = l+1
        ht = m*(y-1) + 1
    else:
        ht = (geom(F)['vert'] + 1) + l*(y-1)
    return (-1)**(ht-1)

def conj_partition(lam, k, nk):
    "conjugate inside the k x nk box -> a partition with <= nk parts each <= k"
    lam = list(lam) + [0]*(k-len(lam))
    return tuple(sum(1 for a in lam if a >= j) for j in range(1, nk+1))

def canon(D):
    "canonical form of a cylindric diagram up to row rotation + column translation"
    best = None
    lo, hi, x, y = list(D.lo), list(D.hi), D.x, D.y
    for t in range(y):
        nlo = [D.LO(t+r) for r in range(y)]
        nhi = [D.HI(t+r) for r in range(y)]
        s = nlo[0]
        key = (x, y, tuple(v-s for v in nlo), tuple(v-s for v in nhi))
        if best is None or key < best: best = key
    return best

def kp_chi(D):
    "chi_{lam/d/mu}(nu) for all nu: chains of cylindric ribbons, weight prod (-1)^{ht-1}"
    from ako import mn_rule
    return mn_rule(D, weight=kp_weight)

def kp_weight_fixed(R, mult):
    """KP's sign, but giving a fully-wound ('all circular ribbons') step the multiplicity
    that the matrix element <v^lam, P*_r v_mu> actually carries."""
    dd = decompose(R)
    if dd is None: return 0
    l, F = dd
    if is_loop_ribbon(F):
        return (-1)**((l+1)*(R.y-1)) * mult
    return (-1)**(geom(F)['vert'] + l*(R.y-1))
