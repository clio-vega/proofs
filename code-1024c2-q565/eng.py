"""Independent engine for Q565 c2.  Pure Fraction arithmetic.

psi_pair : Psi^C_{(1^r),theta}(M) via the alternant expansion with B(T) in the
           PAIR FORM (t1_pairform.py, verified this morning):
               B(T) = N(T)/D(T),
               N(T) = prod_{a<=b, {a,b} cap T != {}} (w_a w_b - 1)
               D(T) = e^r * prod_{a in T} w_a * prod_{a in T, t notin T}(w_a - w_t)
           A(T) = det(x_i^{m_j}) / prod_{i<k}(x_k - x_i),  m_j = theta_j + j - 1.

Cross-checked against psi_headsum (eq:Psidef, sum over ordered r-tuples) from
code-1023-q559/psi.py -- a different mechanism entirely.
"""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from itertools import combinations, permutations


def perm_sign(p):
    s, seen = 1, [False]*len(p)
    for i in range(len(p)):
        if seen[i]: continue
        j, L = i, 0
        while not seen[j]:
            seen[j] = True; j = p[j]; L += 1
        if L % 2 == 0: s = -s
    return s


def alternant(x, m):
    """det(x_i^{m_j})_{i,j} / prod_{i<k}(x_k-x_i).  m_j may be negative."""
    r = len(m)
    num = F(0)
    for p in permutations(range(r)):
        t = F(perm_sign(p))
        for i in range(r):
            t *= x[i]**m[p[i]]
        num += t
    den = F(1)
    for i in range(r):
        for k in range(i+1, r):
            den *= (x[k]-x[i])
    return num/den


def N_pair(w, T):
    n = len(w); Ts = set(T); out = F(1)
    for a in range(n):
        for b in range(a, n):
            if a in Ts or b in Ts:
                out *= (w[a]*w[b]-1)
    return out


def D_pair(w, T):
    n = len(w); r = len(T); Ts = set(T)
    e = F(1)
    for x in w: e *= x
    out = e**r
    for a in T: out *= w[a]
    for a in T:
        for t in range(n):
            if t not in Ts: out *= (w[a]-w[t])
    return out


def terms(w, theta):
    """yield (T, B(T), A(T)) over all r-subsets T, exactly."""
    n = len(w); r = len(theta)
    m = [theta[j]+j for j in range(r)]   # m_j = theta_j + j - 1 with j from 1
    m = [theta[j]+(j+1)-1 for j in range(r)]
    for T in combinations(range(n), r):
        d = D_pair(w, T)
        if d == 0:
            raise ZeroDivisionError(f"degenerate D at T={T}")
        B = N_pair(w, T)/d
        A = alternant([w[i] for i in T], m)
        yield T, B, A


def psi_pair(w, theta):
    return sum((B*A for _, B, A in terms(w, theta)), F(0))


def psi_head(w, theta):
    from psi import psi_headsum
    return psi_headsum(w, theta)


def admissible(w, require_recip_free=False):
    """distinct, nonzero; optionally no w_a w_b = 1 (a<=b, so also w_a^2 != 1)."""
    n = len(w)
    if any(x == 0 for x in w): return False
    if len(set(w)) != n: return False
    if require_recip_free:
        for a in range(n):
            for b in range(a, n):
                if w[a]*w[b] == 1: return False
    return True
