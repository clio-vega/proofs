"""Q484 engine: Motegi's type B/C trigonometric weight functions, from the
symmetrization formula eq (Motherfunction) of 2610.15344 (L473-490).

Pure SymPy -- no combinatorics library.  Every formula transcribed from the
e-print source /tmp/motegi/typebcweightfunctions_arxivversionfinal-3.tex,
line numbers in the docstrings.

Conventions.  Subsets of [N]={1..N} are python tuples of 1-indexed ints,
strictly increasing.  I <-> strict partition lambda(I) = (i_r > ... > i_1).
"""
import sympy as sp
from itertools import permutations, combinations
from math import factorial

q = sp.symbols('q')


def wvars(N):
    return list(sp.symbols(f'w1:{N+1}'))


# ---------------------------------------------------------------- Motherfunction
def F(N, n, r, u, w, I, beta, gamma):
    """F_{N,n,r}(u|w|I), eq (Motherfunction), L474-487.

    u : list of n expressions (the spectral parameters, in order)
    w : list of N expressions
    I : increasing tuple (i_1 < ... < i_r), 1-indexed
    """
    assert len(u) == n and len(w) == N and len(I) == r
    total = sp.Integer(0)
    for sigma in permutations(range(n)):
        us = [u[sigma[j]] for j in range(n)]          # sigma acts on the u's
        t = sp.Integer(1)
        # line 1 of the bracket:  prod_{j=1}^r {(1-q)(u_j^2-1) prod_{k=i_j+1}^N (u_j-q w_k) prod_{k=1}^{i_j-1}(u_j-w_k)}
        for j in range(r):
            ij = I[j]
            t *= (1 - q) * (us[j]**2 - 1)
            for k in range(ij + 1, N + 1):
                t *= (us[j] - q * w[k - 1])
            for k in range(1, ij):
                t *= (us[j] - w[k - 1])
        # line 2a: prod_{1<=j<k<=r} (q u_j - u_k)/(u_j - u_k) (u_j u_k - 1)
        for j in range(r):
            for k in range(j + 1, r):
                t *= (q * us[j] - us[k]) / (us[j] - us[k]) * (us[j] * us[k] - 1)
        # line 2b: prod_{j=r+1}^n prod_{k=1}^N (u_j - q w_k)
        for j in range(r, n):
            for k in range(N):
                t *= (us[j] - q * w[k])
        # line 3a: prod_{j=r+1}^n prod_{k=1}^r (q u_j - u_k)/(u_j - u_k) (u_j u_k - 1)
        for j in range(r, n):
            for k in range(r):
                t *= (q * us[j] - us[k]) / (us[j] - us[k]) * (us[j] * us[k] - 1)
        # line 3b: prod_{r+1<=j<k<=n} (u_j u_k - q)
        for j in range(r, n):
            for k in range(j + 1, n):
                t *= (us[j] * us[k] - q)
        # line 3c: prod_{j=r+1}^n (beta u_j - gamma)
        for j in range(r, n):
            t *= (beta * us[j] - gamma)
        total += t
    pref = sp.Integer(1)
    for k in range(n):
        pref *= u[k]**(-(N + n))
    return sp.together(pref * total / factorial(n - r))


def EC(N, n, r, u, w):
    """Normalization factor E^C_{N,n,r}(u,w), eq (typeCnormalizationfactor), L1057-1062."""
    num = sp.Integer(1)
    for j in range(n):
        for k in range(n):
            num *= (1 - q * u[k] / u[j])
    for j in range(n):
        for k in range(N):
            num *= (1 - q * u[j] * w[k])
    den = sp.Integer(1)
    for j in range(n):
        for k in range(j + 1, n):
            den *= (1 - q * u[j] * u[k])
    for j in range(N):
        for k in range(j, N):
            den *= (1 - q * w[j] * w[k])
    return num / den / (1 - q)**(n - r)


def WX(N, n, r, u, w, I, X):
    """WC (X='C') resp. WB (X='B'), L494-505."""
    if X == 'C':
        return F(N, n, r, u, w, I, sp.Integer(1), sp.Integer(0))
    elif X == 'B':
        pre = sp.Integer(1)
        for k in range(n):
            pre *= 1 / (1 + 1 / u[k])
        return pre * F(N, n, r, u, w, I, sp.Integer(1), sp.Integer(-1))
    raise ValueError(X)


def Atilde(N, w, M, I, X='C', simplify=True):
    """Atilde_{M,I} = widetilde{WX}_{N,n,|I|}(w_M | w | I) with n >= N.

    Computed via the stability relation eq:typeC-stability-s-ge-r (L2822-2828):
    for s=|M| >= r=|I| the value equals widetilde{WX}_{N,s,r}(w_{m_1},..,w_{m_s}|w|I),
    i.e. NO trailing 1's -- which is what makes the (q u_j - u_k)/(u_j - u_k)
    factors non-singular at generic w.
    For s < r the value is 0 (Prop L2843, eq:typeC-lower-cardinality-vanishing).
    """
    s, r = len(M), len(I)
    if s < r:
        return sp.Integer(0)
    n = s
    if n == 0:
        # F_{N,0,0}=WC_{N,0,0}=WB_{N,0,0}=1 (L406-409) and "the normalized
        # functions for n=0 are defined by the same normalization formulas"
        # (L409-411), so the value is 1/E^C_{N,0,0} = prod_{j<=k}(1-q w_j w_k).
        return sp.cancel(sp.Integer(1) / EC(N, 0, 0, [], w))
    u = [w[m - 1] for m in M]
    val = WX(N, n, r, u, w, I, X) / EC(N, n, r, u, w)
    if X == 'B':
        # the WB normalization factor is folded into WX already; E^C is the
        # common normalizer, cf. typeBnormalizeddefinition L1358-1366
        pass
    return sp.cancel(sp.together(val)) if simplify else val


# ------------------------------------------------- independent diagonal formula
def diag_rootproduct(N, w, I):
    """prod_{alpha in Inv(lambda)}(1-e^{-alpha}) prod_{alpha in T(lambda)}(1-q e^{alpha}),
    eq (normalizationWC) L1078, written out via L1017-1022 and L1041-1046.

    This is an INDEPENDENT mechanism: a product over positive roots of C_N,
    not a sum over S_n.  Used as the control on the engine.
    """
    r = len(I)
    Ic = [j for j in range(1, N + 1) if j not in I]
    P = sp.Integer(1)
    # prod_{1<=j<=k<=r} (1 - w_{i_j}^{-1} w_{i_k}^{-1})
    for a in range(r):
        for b in range(a, r):
            P *= (1 - 1 / (w[I[a] - 1] * w[I[b] - 1]))
    # prod_{j in I^c} prod_{l : i_l > j} (1 - w_{i_l}^{-1} w_j)
    for j in Ic:
        for il in I:
            if il > j:
                P *= (1 - w[j - 1] / w[il - 1])
    # prod_{j<=k, j,k in I^c} (1 - q w_j w_k)
    for a in range(len(Ic)):
        for b in range(a, len(Ic)):
            P *= (1 - q * w[Ic[a] - 1] * w[Ic[b] - 1])
    # prod_{j in I^c} prod_{l : i_l < j} (1 - q w_{i_l}^{-1} w_j)
    for j in Ic:
        for il in I:
            if il < j:
                P *= (1 - q * w[j - 1] / w[il - 1])
    return P


# ------------------------------------------------------------------- the order
def geq(M, I):
    """Motegi's partial order M >= I, L640-650:  s>=r and m_{s-r+j} >= i_j for all j."""
    s, r = len(M), len(I)
    if s < r:
        return False
    return all(M[s - r + j] >= I[j] for j in range(r))


def subsets(N):
    out = []
    for k in range(N + 1):
        out.extend(combinations(range(1, N + 1), k))
    return out


def lam_size(I):
    """|lambda(I)| = sum of the parts = sum of I."""
    return sum(I)
