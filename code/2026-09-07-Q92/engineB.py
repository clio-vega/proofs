"""Engine B: the CLOSED FORM for [R_e,R_f] in the Maya/Clifford picture.
Code-disjoint from engineA: no partitions, no border strips, only bead moves
and explicit psi/psi* signs.  Partitions appear only at the two interfaces.
"""
import sympy as sp
t = sp.Symbol('t')

W = 24   # window: all integers < -W are beads, never touched

def maya(lam):
    """explicit part M cap [-W,oo) of the Maya set of lam."""
    lam = list(lam) + [0]*(W - len(lam))
    return frozenset(lam[j-1] - j for j in range(1, W+1))

def unmaya(M):
    a = sorted(M, reverse=True)
    lam = [a[j-1] + j for j in range(1, len(a)+1)]
    while lam and lam[-1] == 0: lam.pop()
    return tuple(lam)

def above(M, a):     # #{x in M : x > a}
    return sum(1 for x in M if x > a)

def ann(state, a):
    """psi*_a on (sign, M) -> (sign, M) or None"""
    s, M = state
    if a not in M: return None
    return (s * (-1)**above(M, a), M - {a})

def cre(state, a):
    """psi_a on (sign, M) -> (sign, M) or None"""
    s, M = state
    if a in M: return None
    return (s * (-1)**above(M, a), M | {a})

def occ(M, a):  return 1 if a in M else 0
def Ncount(M, lo, hi):   # #(M cap (lo,hi)) open interval
    return sum(1 for j in range(lo+1, hi) if j in M)

def closed_form(lam, e, f):
    """[R_e,R_f] s_lam  via  -(1+t)/t * T1  -  (t^2-1)/t * T2."""
    M = maya(lam)
    out = {}
    def add(mu, c):
        out[mu] = sp.expand(out.get(mu, 0) + c)

    # ---- T1 : resonant, an (e+f)-move dressed by (n_{b+f} - n_{b+e}) ----
    pref1 = -(1 + t)/t
    for b in range(-W+1, W):
        if b + e + f >= W: continue
        d = occ(M, b+f) - occ(M, b+e)
        if d == 0: continue
        st = ann((sp.Integer(1), M), b)
        if st is None: continue
        st = cre(st, b+e+f)
        if st is None: continue
        s, Mp = st
        coeff = pref1 * s * (-t)**Ncount(M, b, b+e+f) * d
        add(unmaya(Mp), coeff)

    # ---- T2 : quartic, a genuine two-bead move ----
    pref2 = -(t**2 - 1)/t
    for b in range(-W+1, W):
        for c in range(-W+1, W):
            if b == c: continue
            if b+e >= W or c+f >= W: continue
            k = (1 if c < b+e < c+f else 0) - (1 if c < b < c+f else 0)
            if k == 0: continue
            st = ann((sp.Integer(1), M), c)          # psi*_c
            if st is None: continue
            st = ann(st, b)                          # psi*_b
            if st is None: continue
            st = cre(st, c+f)                        # psi_{c+f}
            if st is None: continue
            st = cre(st, b+e)                        # psi_{b+e}
            if st is None: continue
            s, Mp = st
            coeff = (pref2 * k * s
                     * (-t)**Ncount(M, b, b+e) * (-t)**Ncount(M, c, c+f))
            add(unmaya(Mp), coeff)

    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

if __name__ == "__main__":
    for (lam, e, f) in [((), 1, 2), ((1,), 1, 2), ((), 1, 3), ((), 2, 3)]:
        print(lam, e, f, closed_form(lam, e, f))
