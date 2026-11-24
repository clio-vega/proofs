"""Planar geometry for tiles: phi from localenum (= engine Geo(SIGN=-1)); interior overlap by separating axes."""
import cmath, math
ZETA = cmath.exp(1j*math.pi/3); z0 = cmath.exp(1j*math.pi/12); z1 = cmath.exp(-1j*math.pi/12)
def phi(P): return z0*(float(P[0]) + float(P[1])*ZETA) + z1*(float(P[2]) + float(P[3])*ZETA)
def poly(t):
    zs = [phi(P) for P in t[1]]; c = sum(zs)/len(zs)
    return sorted(zs, key=lambda z: cmath.phase(z - c))
def overlap(s, t, eps=1e-7):
    """True iff the open interiors of convex tiles s, t intersect."""
    if s == t: return True
    A, B = poly(s), poly(t)
    for P in (A, B):
        for i in range(len(P)):
            e = P[(i+1) % len(P)] - P[i]; n = complex(-e.imag, e.real)
            pa = [(z*n.conjugate()).real for z in A]; pb = [(z*n.conjugate()).real for z in B]
            if max(pa) <= min(pb) + eps or max(pb) <= min(pa) + eps: return False
    return True
def conflict(S, T):
    """pairs (s,t), s in S, t in T, distinct tiles with overlapping interiors (cannot coexist in one tiling)"""
    return [(s, t) for s in S for t in T if s != t and overlap(s, t)]
