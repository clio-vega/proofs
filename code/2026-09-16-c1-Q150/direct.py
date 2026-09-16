"""
Q150: independent NUMERIC commutator checker, written from Convention conv:def of
the Q149 paper.  No reuse of solver.py's algebra.

Maya set = (all sites < FLOOR are beads) u (beads in [FLOOR, TOP)).
Initial states: FLOOR..0 all beads, T an arbitrary subset of [0, Lw).
FLOOR = -2g is deep enough that no bead below it ever acquires a legal move
(a hole can first appear only at a site >= -f, and a second move can originate
only at a site >= -g).  Landing sites up to TOP = Lw+g are tracked exactly, so
no legal move is ever dropped.
"""
from itertools import product

def words(k):
    return [tuple(w) for w in product((0,1), repeat=k)]

def apply_R(vec, e, W, FLOOR, TOP):
    out = {}
    for S, c in vec.items():
        for b in range(FLOOR, TOP):
            if b not in S:      continue
            if (b+e) in S:      continue
            if b+e >= TOP:      raise RuntimeError("window too small: move escapes")
            u = tuple(1 if (b+i) in S else 0 for i in range(1, e))
            w = W[u]
            if w == 0:          continue
            S2 = S - {b} | {b+e}
            out[S2] = out.get(S2, 0) + c*w
    return out

def commutator_nonzero(e, f, W, Wb, Lw=14, cap=40, mod=None):
    """All nonzero matrix elements of [R_e^W, R_f^Wb] over the window states."""
    g = e+f
    FLOOR, TOP = -2*g, Lw + 2*g
    deep = frozenset(range(FLOOR, 0))
    bad = []
    for mask in range(1 << Lw):
        S0 = deep | frozenset(i for i in range(Lw) if (mask >> i) & 1)
        v = {S0: 1}
        A = apply_R(apply_R(v, f, Wb, FLOOR, TOP), e, W,  FLOOR, TOP)   # R_e R_f
        B = apply_R(apply_R(v, e, W,  FLOOR, TOP), f, Wb, FLOOR, TOP)   # R_f R_e
        for k in set(A) | set(B):
            val = A.get(k, 0) - B.get(k, 0)
            if mod: val %= mod
            if val != 0:
                bad.append((sorted(x for x in S0 if x >= FLOOR+g),
                            sorted(x for x in k if x >= FLOOR+g), val))
                if cap and len(bad) >= cap: return bad
    return bad
