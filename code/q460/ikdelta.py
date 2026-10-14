"""
Q460 -- Izergin-Korepin 19-vertex scattering factors Delta_{U,V}(x,y;z).

Transcribed EXACTLY from Garbali-Guo-Wheeler arXiv:2412.18085 v1,
eq:IK-delta1 (= 3.8.37), eq:IK-delta2, eq:IK-delta3, with the symmetry
eq:delta-sym  Delta_{U,V}(x,y;z) = Delta_{V,U}(y,x;z).

A weight-2 row U of a 2-permutation matrix is encoded as its sorted
endpoint pair (u1,u2):
    "2 at column k"      ->  (k,k)        [doubled endpoint]
    "1 at k, 1 at l<k"   ->  (k,l), k<l   [an arc]
Column sums are at most 2, so a collision (shared column) can only occur
between two arcs, never involving a doubled endpoint.  This is exactly why
eq:IK-delta2/3 only ever display U=11, V=11 diagrams.
"""
import sympy as sp

q = sp.Symbol('q')

def Z(k):
    return sp.Symbol('z%d' % k)

# ---------------------------------------------------------------- alphabet
# eq:IK-delta1 -- tie-free.  The three diagram groups are, in the paper's
# order: (A) all of U strictly before all of V; (B) V strictly nested
# inside U; (C) U and V cross.
def dA(x, y):            # 4 diagrams: 2/2, 2/11, 11/2, 11/11  (disjoint)
    return (y - q**2*x)*(y - q**3*x) / ((y - x)*(y - q*x))

def dB(x, y):            # 2 diagrams: 1_1/2, 1__1/_11_        (V nested in U)
    return (y - q**2*x)*(q**2*y - x) / (y - x)**2

def dC(x, y):            # 1 diagram: 1_1_/_1_1               (crossing)
    return (y - q**2*x)**2 * (x - q*y) / ((y - x)**2 * (y - q*x))

# eq:IK-delta2 -- exactly one shared column, at absolute position k.
def d2i(x, y, k):        # u1=v1=k, u2<v2
    zk = Z(k)
    return (q**2*x - y)*(x*y + q**2*x*y - q**3*y*zk - q**3*x*zk) \
           / (q*(1-q)*(y-x)*(y-q*x)*zk)

def d2ii(x, y, k):       # u2=v1=k  (U closes where V opens)
    zk = Z(k)
    return (y - q**2*x)*(x*y + q**2*x*y - q**2*y*zk - q**3*x*zk) \
           / (q*(1-q)*(y-x)**2*zk)

def d2iii(x, y, k):      # u2=v2=k, u1<v1
    zk = Z(k)
    return (q**2*x - y)*(x*y + q**2*x*y - q**2*y*zk - q**2*x*zk) \
           / ((1-q)*(y-x)*(y-q*x)*zk)

# eq:IK-delta3 -- both columns shared: u1=v1=k, u2=v2=l.
def d3(x, y, k, l):
    zk, zl = Z(k), Z(l)
    return (x*y + q**2*x*y - q**3*y*zk - q**3*x*zk) \
           * (x*y + q**2*x*y - q**2*y*zl - q**2*x*zl) \
           / (q*(1-q)**2*(y-q*x)*(x-q*y)*zk*zl)

LISTED = ('A', 'B', 'C', '2i', '2ii', '2iii', '3')   # 7 rational functions

# ------------------------------------------------------------- classifier
def pattern(u, v):
    """Return (name, args) for the LISTED case of (U,V), or None if (U,V)
    is not listed and must be obtained from (V,U) via eq:delta-sym."""
    u1, u2 = u
    v1, v2 = v
    shared = sorted(set([u1, u2]) & set([v1, v2])) if u1 != u2 and v1 != v2 else []
    if not shared:
        if u2 < v1:                        return ('A', ())
        if u1 < v1 and v2 < u2:            return ('B', ())
        if u1 < v1 < u2 < v2:              return ('C', ())
        return None                        # A^T, B^T, C^T
    if len(shared) == 2:
        return ('3', (u1, u2))             # u1=v1=k, u2=v2=l
    k = shared[0]
    if u1 == v1 and u2 < v2:               return ('2i', (k,))
    if u2 == v1:                           return ('2ii', (k,))
    if u2 == v2 and u1 < v1:               return ('2iii', (k,))
    return None                            # the transposes

FN = {'A': dA, 'B': dB, 'C': dC, '2i': d2i, '2ii': d2ii, '2iii': d2iii, '3': d3}

def Delta(u, v, x, y, used_sym=None):
    """Delta_{U,V}(x,y;z), U=u=(u1,u2), V=v=(v1,v2)."""
    p = pattern(u, v)
    if p is not None:
        name, args = p
        if used_sym is not None: used_sym.add(name)
        return FN[name](x, y, *args)
    p = pattern(v, u)
    if p is None:
        raise ValueError("unclassified pair %r %r" % (u, v))
    name, args = p
    if used_sym is not None: used_sym.add(name + '^T')
    return FN[name](y, x, *args)          # eq:delta-sym

def row_to_endpoints(row):
    """A row of a 2-permutation matrix (list of 0/1/2) -> (u1,u2), 1-indexed."""
    pts = []
    for j, e in enumerate(row, start=1):
        pts += [j]*e
    assert len(pts) == 2, row
    return (pts[0], pts[1])
