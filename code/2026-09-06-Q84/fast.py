"""fast.py -- integer-coefficient operators at t=-1, for the order computation.

All three operators act on states = dict{partition: int}.

  op_p(.,a)   : multiplication by p_a.  MN rule, abacus.
  op_N(.,e)   : N_e = "add a size-e broken strip with EXACTLY 2 components",
                weight (-1)^ht.       [ = d/dt E_e at t=-1 ]
  op_Psi(.,e) : Psi_{e,1} = "add a CONNECTED e-ribbon", weight ht*(-1)^ht.
                [ = -d/dt R_e at t=-1 ]

Since g_e = R_e + E_e has ord 0 (multiplication by P_{(e)}(X;-t)), and
d/dt g_e|_{-1} is again a multiplication operator, we have
        N_e = M[ d/dt P_(e)|_{-1} ] + Psi_{e,1},
so N_e and Psi_{e,1} differ by an order-0 operator and
        [N_e, X] = [Psi_{e,1}, X]   for every X commuting with multiplication.
op_N and op_Psi are built from disjoint code paths (broken-shape enumeration vs
abacus bead moves); their first brackets agreeing is a real cross-check.
"""
from functools import lru_cache

L = 40   # runners; must exceed len(lam) for every lam that occurs

def trim(lam):
    lam = tuple(lam)
    while lam and lam[-1] == 0: lam = lam[:-1]
    return lam

def beta(lam):
    lam = trim(lam) + (0,) * (L - len(trim(lam)))
    return tuple(sorted(lam[i] + L - 1 - i for i in range(L)))

def unbeta(B):
    B = sorted(B, reverse=True)
    return trim(tuple(B[i] - (L - 1 - i) for i in range(L)))

@lru_cache(maxsize=None)
def ribbons(lam, e):
    """connected e-ribbon additions: tuple of (mu, ht)."""
    B = set(beta(lam)); out = []
    for b in sorted(B):
        if b + e in B: continue
        h = sum(1 for c in B if b < c < b + e)
        out.append((unbeta((B - {b}) | {b + e}), h))
    return tuple(out)

@lru_cache(maxsize=None)
def broken2(lam, e):
    """size-e skew shapes mu/lam, no 2x2, EXACTLY 2 edge-components: (mu, ht)."""
    lam = trim(lam)
    nr = len(lam) + e
    lamp = lam + (0,) * (nr - len(lam))
    res = []
    def rec(i, prev, rem, mu):
        if rem == 0:
            res.append(tuple(mu) + tuple(lamp[i:]))
            return
        if i == nr: return
        hi = min(prev, lamp[i] + rem)
        for v in range(lamp[i], hi + 1):
            rec(i + 1, v, rem - (v - lamp[i]), mu + [v])
    rec(0, 10**9, e, [])
    out = []
    for mu in res:
        mu = trim(mu)
        cells = set((i + 1, j + 1) for i in range(len(mu))
                    for j in range(lamp[i] if i < len(lamp) else 0, mu[i]))
        if any((i, j) in cells and (i+1, j) in cells and (i, j+1) in cells
               and (i+1, j+1) in cells for (i, j) in cells):
            continue
        seen = set(); m = 0; ht = 0
        for c in cells:
            if c in seen: continue
            m += 1
            if m > 2: break
            st = [c]; seen.add(c); rows = set()
            while st:
                (a, b) = st.pop(); rows.add(a)
                for nb in ((a+1, b), (a-1, b), (a, b+1), (a, b-1)):
                    if nb in cells and nb not in seen:
                        seen.add(nb); st.append(nb)
            ht += len(rows) - 1
        if m == 2:
            out.append((mu, ht))
    return tuple(out)

def _acc(out, mu, c):
    out[mu] = out.get(mu, 0) + c

def op_p(state, a):
    out = {}
    for lam, w in state.items():
        for (mu, h) in ribbons(lam, a): _acc(out, mu, w * (-1)**h)
    return {k: v for k, v in out.items() if v}

def op_N(state, e):
    out = {}
    for lam, w in state.items():
        for (mu, h) in broken2(lam, e): _acc(out, mu, w * (-1)**h)
    return {k: v for k, v in out.items() if v}

def op_Psi(state, e):
    out = {}
    for lam, w in state.items():
        for (mu, h) in ribbons(lam, e): _acc(out, mu, w * h * (-1)**h)
    return {k: v for k, v in out.items() if v}
