"""Q442 -- sector matching.  PROVE 2026-10-13 c2.

Three engines with NO shared code path beyond poly.py's arithmetic:

  A  determinant engine : det_k(f_{i-j} - f_{i+j})                      (symmetric functions)
  B  Warnaar engine     : [t^0] sum_{y in Z^k} det_k(t^{My^2-jy}(edot.. - edot..)) * (x1..xn)^k
                          (Laurent polynomials in x_i^{+-1}, independent of A)
  C  path/DP engine     : sum_lambda (weight) s_{lambda'}  via vertical-strip lattice paths
                          (no determinants, no elementary symmetric functions)

Verdicts printed with their domain (k, ell, n, w) attached, every time.
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1013c2-q442')
from poly import *


# =============================================================== engine A
def engineA(k, n):
    """det_{1<=i,j<=k}( f_{i-j} - f_{i+j} ).   HKKO eq:even-sp RHS (h=k)."""
    E = elem(n)
    M = [[sub(f_r(E, i - j, n), f_r(E, i + j, n)) for j in range(1, k + 1)]
         for i in range(1, k + 1)]
    return det(M, n)

def engineA_transposed(k, n):
    """det_{1<=i,j<=k}( f_{-i+j} - f_{i+j} ).  Must agree with engineA (transpose)."""
    E = elem(n)
    M = [[sub(f_r(E, -i + j, n), f_r(E, i + j, n)) for j in range(1, k + 1)]
         for i in range(1, k + 1)]
    return det(M, n)


# =============================================================== engine B
def edot_table(n):
    """edot_r = e_r(x_1,..,x_n,x_1^{-1},..,x_n^{-1}) for 0<=r<=2n, as Laurent polys."""
    # product_{i} (1 + x_i z)(1 + x_i^{-1} z); collect by z-degree
    cur = [const(1, n)]                     # cur[d] = coeff of z^d
    for i in range(n):
        ep = [0] * n; ep[i] = 1;  xp = {tuple(ep): 1}
        em = [0] * n; em[i] = -1; xm = {tuple(em): 1}
        nxt = [{} for _ in range(len(cur) + 2)]
        for d, c in enumerate(cur):
            if not c: continue
            nxt[d]     = add(nxt[d], c)
            nxt[d + 1] = add(nxt[d + 1], add(mul(c, xp), mul(c, xm)))
            nxt[d + 2] = add(nxt[d + 2], mul(mul(c, xp), xm))
        cur = nxt
    return cur[:2 * n + 1]

def edot(tab, r, n):
    if r < 0 or r > 2 * n: return {}
    return tab[r]

def engineB(k, ell, n, ymax=None):
    """[t^0] of Warnaar Eq_generalise, times (x1..xn)^k.

    Computes G with t symbolic (coefficients in Z[t]) and extracts t^0 afterwards,
    so it does NOT presuppose the t=0 theorem.  Also returns G at t=1.
    """
    M = k + ell + 1
    N = 2 * M
    tab = edot_table(n)
    if ymax is None:
        ymax = (n + 2 * k) // N + 2
    total = {}
    for y in itertools.product(range(-ymax, ymax + 1), repeat=k):
        rows = []
        dead = False
        for i in range(1, k + 1):
            yi = y[i - 1]
            row = []
            for j in range(1, k + 1):
                ent = sub(edot(tab, n - i + j - N * yi, n), edot(tab, n + i + j - N * yi, n))
                ex = M * yi * yi - j * yi
                row.append(Tscal_t(Tlift(ent), ex) if ent else {})
            rows.append(row)
        d = Tdet(rows, n)
        if d: total = Tadd(total, d)
    pref = {tuple([k] * n): 1}
    t0 = mul(Tcoeff_t(total, 0), pref)
    # G at t=1
    g1 = {}
    for m, cc in total.items():
        v = sum(cc.values())
        if v: g1[m] = v
    g1 = mul(g1, pref)
    tdegs = sorted({d for cc in total.values() for d in cc})
    return t0, g1, tdegs


# =============================================================== engine C
def strip_successors(lam, m):
    """all mu with mu_i - lam_i in {0,1} and mu weakly decreasing (vertical strips)."""
    out = []
    for S in range(1 << m):
        mu = tuple(lam[i] + ((S >> i) & 1) for i in range(m))
        ok = all(mu[i] >= mu[i + 1] for i in range(m - 1))
        if ok: out.append((mu, bin(S).count('1')))
    return out

def c_minus(lam, m, w):
    """HKKO eq:wt_C  c^-_{m,w}(lambda),  m = 2h."""
    h = m // 2
    if all(lam[2 * i] == lam[2 * i + 1] for i in range(h)):
        return 1
    if lam[0] - lam[m - 1] == w and all(lam[2 * i + 1] == lam[2 * i + 2] for i in range(h - 1)):
        return -1
    return 0

def engineC(k, n, w, bounded):
    """DP over vertical strips.  Returns dict: endpoint lambda -> polynomial sum_T x^T.

    bounded=True  : every intermediate lambda obeys lam_1 - lam_m <= w  (cylindric tableaux)
    bounded=False : no width constraint                                 (ordinary SSYT)
    """
    m = 2 * k
    start = tuple([0] * m)
    dp = {start: const(1, n)}
    for v in range(n):
        xv = [0] * n; xv[v] = 1
        nxt = {}
        for lam, p in dp.items():
            for mu, a in strip_successors(lam, m):
                if mu[0] > n: continue
                if bounded and mu[0] - mu[m - 1] > w: continue
                mono = [0] * n; mono[v] = a
                term = mul(p, {tuple(mono): 1})
                nxt[mu] = add(nxt.get(mu, {}), term)
        dp = nxt
    return dp


def engineC_sums(k, n, w):
    m = 2 * k
    dpU = engineC(k, n, w, bounded=False)
    dpB = engineC(k, n, w, bounded=True)
    Ainf, Aplus_bd, Bbd, kappa0, overwide = {}, {}, {}, {}, {}
    for lam, p in dpU.items():
        cm = c_minus(lam, m, w)
        if cm == 1:
            Ainf = add(Ainf, p)
            if lam[0] - lam[m - 1] > w: overwide = add(overwide, p)
        if lam[0] - lam[m - 1] <= w and cm:
            kappa0 = add(kappa0, scal(cm, p))
    Binf = {}
    for lam, p in dpU.items():
        if c_minus(lam, m, w) == -1: Binf = add(Binf, p)
    for lam, p in dpB.items():
        cm = c_minus(lam, m, w)
        if cm == 1: Aplus_bd = add(Aplus_bd, p)
        if cm == -1: Bbd = add(Bbd, p)
    return dict(Ainf=Ainf, A=Aplus_bd, B=Bbd, kappa0=kappa0,
                overwide=overwide, Binf=Binf)
