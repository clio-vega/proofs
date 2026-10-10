import sympy as sp, itertools
from core import Engine
from warn import G, t
from walks import AB_from_walks

def l1_report(k, ell, n, verbose=0):
    """For every content alpha occurring in (x..)^k * Eq_generalise, compare
         ||[x^alpha] Eq_gen||_1  (sum of |t-coefficients|)
       against
         TOT_alpha = # cylindric tableaux of ALL shapes la in Par(2k,2ell+2)
                     with content alpha and entries in [n]
       (and A_alpha / B_alpha, the c^-=+1 / c^-=-1 subsets).
       Any alpha with ||.||_1 > TOT_alpha kills every signed monomial statistic."""
    w = 2*ell+2
    E = Engine(n)
    Gt = sp.expand(G(E, k, ell, t) * E.prodx**k)
    P = sp.Poly(Gt, *E.x)
    worst = None; nviol = 0; nmono = 0; nviol_strict = 0
    for mono, coeff in zip(P.monoms(), P.coeffs()):
        co = sp.expand(coeff)
        d = sp.Poly(co, t).as_dict() if co.has(t) else {(0,): co}
        l1  = sum(abs(int(v)) for v in d.values())
        pos = sum(int(v) for v in d.values() if v > 0)
        neg = sum(-int(v) for v in d.values() if v < 0)
        A, B, TOT = AB_from_walks(k, w, mono)
        nmono += 1
        assert A - B == int(co.subs(t, 1)), (mono, A, B, co)      # t=1 control, every monomial
        if pos > A or neg > B: nviol_strict += 1                   # kills stat with c^- signs
        if l1 > TOT:                                               # kills ANY signed stat
            nviol += 1
            gap = l1 - TOT
            if worst is None or gap > worst[0]:
                worst = (gap, mono, {int(kk[0]): int(v) for kk, v in d.items()}, A, B, TOT)
    return dict(k=k, ell=ell, w=w, n=n, nmono=nmono, nviol_c=nviol_strict, nviol_any=nviol, worst=worst)
