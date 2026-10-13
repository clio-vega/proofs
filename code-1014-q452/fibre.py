"""Q452 engine: fibres T(k,ell,(1^n)) in the CHIP formulation.

Reformulation (proved in the paper): a single box at cyclic slot i is legal from x
iff m_{i-1}(x) >= 1, and it reads exactly that multiplicity, decrementing it by 1 and
incrementing m_i by 1.  So the walk is a CHIP-HOPPING on the cycle Z_h with w chips:
   state  m = (m_1,...,m_h),  sum m_i = w,  start m^(0) = (0,...,0,w)
   move at SOURCE j (needs m_j>=1): read the pair (j, m_j), then m_j -=1, m_{j+1} +=1,
   and lambda_{j+1} += 1   (indices in 1..h, h+1 == 1).
The weight read is g[(j, m_j_before)] -- j is Korff's cyclic index J(e_i)={i-1}.
"""
import itertools, sympy as sp
from functools import lru_cache

def sign_hkko(lam, h, w):
    """HKKO c^-_{h,w}(lambda). lam is a tuple of length h, weakly decreasing, >=0."""
    k = h // 2
    if all(lam[2*i] == lam[2*i+1] for i in range(k)):      # lam_{2i-1}=lam_{2i}, 1-indexed
        return 1
    if lam[0] - lam[h-1] == w and all(lam[2*i+1] == lam[2*i+2] for i in range(k-1)):
        return -1
    return 0

def enumerate_fibre(h, w, n):
    """All single-box paths of length n from the origin.
    Returns list of (reads, lam, sign) where reads is a tuple of (j, m) pairs."""
    start = tuple([0]*(h-1) + [w])
    states = {(start, tuple([0]*h)): [()]}      # (m, lam) -> list of read-words
    for _ in range(n):
        nxt = {}
        for (m, lam), words in states.items():
            for j in range(1, h+1):
                if m[j-1] == 0: continue
                i = 1 if j == h else j+1                 # destination slot
                m2 = list(m); m2[j-1] -= 1; m2[i-1] += 1
                l2 = list(lam); l2[i-1] += 1
                key = (tuple(m2), tuple(l2))
                r = (j, m[j-1])
                nxt.setdefault(key, [])
                for wd in words: nxt[key].append(wd + (r,))
        states = nxt
    out = []
    for (m, lam), words in states.items():
        s = sign_hkko(lam, h, w)
        for wd in words:
            out.append((wd, lam, s))
    return out

def sigma_symbolic(h, w, n, gsym):
    """sum_T c^-(T) prod g[(j,m)] as a sympy expression in the gsym dict."""
    tot = sp.Integer(0)
    for reads, lam, s in enumerate_fibre(h, w, n):
        if s == 0: continue
        term = sp.Integer(s)
        for r in reads: term *= gsym[r]
        tot += term
    return sp.expand(tot)

def gvars(h, w):
    return {(j, m): sp.Symbol(f'g_{j}_{m}') for j in range(1, h+1) for m in range(1, w+1)}
