"""Engine 2 -- from the DEFINITION, no tableaux, no charge.
   Hall-Littlewood P_lam is the unique family with
       P_lam = m_lam + (lower in dominance),   <P_lam, P_mu>_t = 0 for lam != mu,
   where  <p_lam, p_mu>_t = delta_{lam,mu} z_lam prod_i (1-t^{lam_i})^{-1}  (Macdonald III).
   Everything is done in the p-basis.  Green polynomials come from p_rho = sum_mu Y^mu_rho P_mu.
"""
import sympy as sp
from functools import lru_cache
t = sp.Symbol('t')

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,)+rest

def zlam(lam):
    from collections import Counter
    c = Counter(lam); z = 1
    for part, m in c.items():
        z *= part**m * sp.factorial(m)
    return sp.Integer(z)

def dominates(lam, mu):
    """lam >= mu in dominance (same size)."""
    sl = sm = 0
    for i in range(max(len(lam), len(mu))):
        sl += lam[i] if i < len(lam) else 0
        sm += mu[i]  if i < len(mu)  else 0
        if sl < sm: return False
    return True

# ---- p-basis vectors: dict  partition(sorted desc tuple) -> coeff ----
def pmul(A, B):
    out = {}
    for a, ca in A.items():
        for b, cb in B.items():
            k = tuple(sorted(a+b, reverse=True))
            out[k] = out.get(k, 0) + ca*cb
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

def h_in_p(r):
    """h_r = sum_{nu |- r} p_nu / z_nu."""
    return {nu: sp.Rational(1,1)/zlam(nu) for nu in partitions(r)}

def hlam_in_p(mu):
    out = {(): sp.Integer(1)}
    for r in mu: out = pmul(out, h_in_p(r))
    return out

def kostka_number(lam, mu):
    """#SSYT(lam, mu) by direct count."""
    cells = [(i,j) for i in range(len(lam)) for j in range(lam[i])]
    k = len(mu); rows=[[None]*lam[i] for i in range(len(lam))]; rem=list(mu); cnt=0
    def rec(idx):
        nonlocal cnt
        if idx == len(cells): cnt += 1; return
        i,j = cells[idx]; lo = 1
        if j>0: lo = max(lo, rows[i][j-1])
        if i>0: lo = max(lo, rows[i-1][j]+1)
        for v in range(lo, k+1):
            if rem[v-1]==0: continue
            rows[i][j]=v; rem[v-1]-=1; rec(idx+1); rem[v-1]+=1; rows[i][j]=None
    rec(0); return cnt

@lru_cache(maxsize=None)
def basis_data(n):
    """Return (parts, m_in_p dict per partition).  m = K^{-1} h with K the Kostka matrix."""
    parts = list(partitions(n))
    K = sp.Matrix(len(parts), len(parts),
                  lambda i, j: kostka_number(parts[i], parts[j]))   # h_mu = sum_lam K_{lam,mu} m_lam
    # CORRECT relations (two convention errors measured 2026-10-09 before this line was right):
    #   s_lam = sum_mu K_{lam,mu} m_mu          =>  s = K m
    #   h_mu  = sum_lam K_{lam,mu} s_lam        =>  h = K^T s
    # hence h = K^T K m  and  m = (K^T K)^{-1} h.
    Kinv = (K.T*K).inv()
    H = [hlam_in_p(mu) for mu in parts]
    M = []
    for i in range(len(parts)):
        acc = {}
        for j in range(len(parts)):
            c = Kinv[i, j]
            if c == 0: continue
            for k, v in H[j].items():
                acc[k] = acc.get(k, 0) + c*v
        M.append({k: sp.together(sp.expand(v)) for k, v in acc.items() if sp.expand(v) != 0})
    return parts, M

def ip(A, B):
    """<A,B>_t in the p-basis."""
    tot = sp.Integer(0)
    for lam, ca in A.items():
        cb = B.get(lam)
        if cb is None: continue
        w = zlam(lam)
        for r in lam: w = w/(1-t**r)
        tot += ca*cb*w
    return sp.cancel(sp.together(tot))

@lru_cache(maxsize=None)
def HL_P(n):
    """Gram-Schmidt in a linear extension of dominance (increasing).  Returns (parts, P_in_p)."""
    parts, M = basis_data(n)
    order = sorted(range(len(parts)), key=lambda i: (sum(1 for j in range(len(parts))
                                                        if dominates(parts[i], parts[j])),))
    Ps = {}
    for i in order:
        cur = dict(M[i])
        for j in order:
            if j == i: break
            Pj = Ps[j]
            c = sp.cancel(ip(M[i], Pj)/ip(Pj, Pj))
            if c == 0: continue
            for k, v in Pj.items():
                cur[k] = sp.cancel(cur.get(k, 0) - c*v)
        Ps[i] = {k: v for k, v in cur.items() if sp.cancel(v) != 0}
    return parts, Ps

@lru_cache(maxsize=None)
def Y_def(mu, rho):
    """Y^mu_rho(t) = [P_mu] p_rho, from the definition-based P."""
    n = sum(mu)
    parts, Ps = HL_P(n)
    # expand p_rho in the P basis: solve linear system in the p-basis
    idx = {p: i for i, p in enumerate(parts)}
    A = sp.zeros(len(parts), len(parts))
    for j, pj in enumerate(parts):
        for k, v in Ps[j].items():
            A[idx[k], j] = v
    b = sp.zeros(len(parts), 1)
    b[idx[tuple(sorted(rho, reverse=True))], 0] = 1
    sol = A.solve(b)
    return sp.cancel(sp.expand(sol[idx[mu], 0]))
