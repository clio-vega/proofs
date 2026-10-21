"""G3: ord_S Atilde^X_{I,I} = [m in I] on S = {w_m = eps}.

Everything at the GENERIC POINT of S: w_{j!=m} and q stay SYMBOLIC.  No numerical
sampling, so the hypersurfaces w_k = +- q^{-1} of rem:generic cannot be hit.
"""
import sys
sys.path.insert(0, '/home/clio/projects/scratch')
import sympy as sp
from q484_engine import q, wvars, Atilde, diag_rootproduct, subsets

t = sp.symbols('t')


def diag_rootproduct_X(N, w, I, X):
    """prod_{Inv(lam)}(1-e^-a) prod_{T(lam)}(1-q e^+a).

    X='C': Delta^+ = {t_k +- t_j, j<k} u {2t_j}, blocks L1017-1046 (typecfirstblock,
           typecsecondblock).
    X='B': Delta^+ = {t_k +- t_j, j<k} u {t_j},  blocks L1309-1335.
    The ONLY difference is the boundary root: 2t_j vs t_j, i.e. h=2 vs h=1.
    """
    if X == 'C':
        return diag_rootproduct(N, w, I)
    r = len(I)
    Ic = [j for j in range(1, N + 1) if j not in I]
    P = sp.Integer(1)
    # Inv: {t_k : k in I}  ->  prod_{j=1}^r (1 - w_{i_j}^{-1})
    for a in range(r):
        P *= (1 - 1 / w[I[a] - 1])
    # Inv: {t_j+t_k : j<k, both in I}  ->  prod_{j<k} (1 - w_{i_j}^-1 w_{i_k}^-1)   [STRICT]
    for a in range(r):
        for b in range(a + 1, r):
            P *= (1 - 1 / (w[I[a] - 1] * w[I[b] - 1]))
    # Inv: {t_k-t_j : j notin I, k in I} -> prod_{j in I^c} prod_{i_l>j} (1 - w_{i_l}^-1 w_j)
    for j in Ic:
        for il in I:
            if il > j:
                P *= (1 - w[j - 1] / w[il - 1])
    # T: {t_k : k notin I}  ->  prod_{j in I^c} (1 - q w_j)
    for j in Ic:
        P *= (1 - q * w[j - 1])
    # T: {t_j+t_k : j<k both notin I} -> prod_{j<k} (1 - q w_j w_k)                 [STRICT]
    for a in range(len(Ic)):
        for b in range(a + 1, len(Ic)):
            P *= (1 - q * w[Ic[a] - 1] * w[Ic[b] - 1])
    # T: {t_k-t_j : j in I, k notin I} -> prod_{j in I^c} prod_{i_l<j} (1 - q w_{i_l}^-1 w_j)
    for j in Ic:
        for il in I:
            if il < j:
                P *= (1 - q * w[j - 1] / w[il - 1])
    return P


def ord_S(expr, w, m, eps):
    """t-adic valuation of expr along S={w_m=eps}, t = w_m - eps, generic point."""
    e = sp.cancel(sp.together(expr.subs(w[m - 1], eps + t)))
    num, den = sp.fraction(e)
    pn = sp.Poly(sp.expand(num), t)
    pd = sp.Poly(sp.expand(den), t)
    if pn.is_zero:
        return None  # identically zero on S: ord = +infinity
    return pn.monoms()[-1][0] - pd.monoms()[-1][0]


def diag_factors(N, I, w, X):
    """The diagonal as a LIST of factors indexed by Inv(lam) u T(lam).
    Returns [(alpha_label, factor)].  X in {'B','C','BC'}; 'BC' is the NON-REDUCED
    lesion in which both t_j and 2t_j are positive roots."""
    r = len(I); Ic = [j for j in range(1, N+1) if j not in I]
    hs = {'B': [1], 'C': [2], 'BC': [1, 2]}[X]
    out = []
    for h in hs:                                    # boundary roots h*t_k
        for k in I:  out.append((f'Inv {h}t_{k}', 1 - w[k-1]**(-h)))
        for k in Ic: out.append((f'T   {h}t_{k}', 1 - q*w[k-1]**h))
    for a in range(r):                              # t_j+t_k, j<k, both in I
        for b in range(a+1, r):
            out.append((f'Inv t_{I[a]}+t_{I[b]}', 1 - 1/(w[I[a]-1]*w[I[b]-1])))
    for a in range(len(Ic)):                        # t_j+t_k, j<k, both out
        for b in range(a+1, len(Ic)):
            out.append((f'T   t_{Ic[a]}+t_{Ic[b]}', 1 - q*w[Ic[a]-1]*w[Ic[b]-1]))
    for j in Ic:                                    # t_k-t_j, j<k, j out, k in
        for k in I:
            if k > j: out.append((f'Inv t_{k}-t_{j}', 1 - w[j-1]/w[k-1]))
    for j in Ic:                                    # t_k-t_j, j<k, j in, k out
        for k in I:
            if k < j: out.append((f'T   t_{j}-t_{k}', 1 - q*w[j-1]/w[k-1]))
    return out
