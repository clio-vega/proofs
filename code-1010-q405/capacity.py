import sympy as sp, itertools
from core import Engine
from hkko import cminus, cyl_schur_transpose
from warn import G, t

def tableau_sides(E, h, w):
    """Return dicts A,B: monomial-exponent-tuple -> number of cylindric tableaux
    with c^-=+1 (A) resp c^-=-1 (B). Name of the set: for each la in Par(2h,w) with
    c^-(la)=+1 (resp -1), the tableaux in CRST(la;2h,w) with entries in [n],
    bucketed by content alpha."""
    m, n = 2*h, E.n
    A, B = {}, {}
    nla = {1:0, -1:0}
    for la in itertools.combinations_with_replacement(range(2*n, -1, -1), m):
        la = tuple(la)
        if la[0]-la[m-1] > w: continue
        c = cminus(la, h, w)
        if c == 0: continue
        s = cyl_schur_transpose(E, la, m, w)
        if s == 0: continue
        nla[c] += 1
        d = (A if c == 1 else B)
        p = sp.Poly(s, *E.x)
        for mono, coeff in zip(p.monoms(), p.coeffs()):
            d[mono] = d.get(mono, 0) + int(coeff)
    return A, B, nla

def capacity_report(h, ell, n, perturb=None, label=''):
    w = 2*ell+2; M = h+ell+1; N = 2*M
    E = Engine(n)
    Gt = sp.expand(G(E, h, ell, t, perturb) * E.prodx**h)
    A, B, nla = tableau_sides(E, h, w)
    P = sp.Poly(Gt, *E.x)
    rows = []
    viol = 0; checked = 0; tot_mono = 0
    for mono, coeff in zip(P.monoms(), P.coeffs()):
        ct = sp.Poly(sp.expand(coeff), t).as_dict() if coeff.has(t) else {(0,): coeff}
        pos = sum(int(v) for k, v in ct.items() if v > 0)
        neg = sum(-int(v) for k, v in ct.items() if v < 0)
        a = A.get(mono, 0); b = B.get(mono, 0)
        checked += 1
        ok_t1 = (a - b == int(sp.expand(coeff).subs(t, 1)))
        ok_cap = (pos <= a) and (neg <= b)
        if not (ok_t1 and ok_cap):
            viol += 1
            if viol <= 8:
                rows.append((mono, {tuple(k):int(v) for k,v in ct.items()}, a, b, ok_t1, ok_cap))
    print(f'--- {label} h={h} ell={ell} (w={w},N={N},M={M}) n={n} ---')
    print(f'    lambdas with c^-=+1: {nla[1]}, with c^-=-1: {nla[-1]}')
    print(f'    monomials of (x..)^h*Eq_gen checked: {checked}; capacity/t=1 violations: {viol}')
    for r in rows: print('      VIOL', r)
    return viol, checked
