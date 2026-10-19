r"""Boundary-divisor pole census for the renormalised structure constants.

DOMAIN (stated before the number, per the standing instrument rule):
  type X in {'C'} or {'B'};  N as given;  m in [N];  eps with eps^h = 1;
  ALL triples (lam, mu, nu) in (2^[N])^3  -- 2^(3N) of them, no sampling;
  q = q0 a fixed rational, w_j = r_j fixed rationals for j != m, w_m = eps + t.
  ord_t is computed exactly over QQ((t)) at truncation K.

Twists chi tested:
  'id'   chi = 1                 -> chi(eps) = 1       (control: no pole expected anywhere)
  'psi'  chi = psi               -> chi(eps) = 1       (control: no pole expected anywhere)
  'xpsi' chi = v*psi(v)          -> chi(eps) = eps     (positive arm at eps = -1)
  'x'    chi = v                 -> chi(eps) = eps     (positive arm at eps = -1; breaks COB)

PREDICATE under test (the proposed S-analogue of Rick's E):
  E_S(lam,mu,nu;m) = { R subset [N]\{m} :  R u {m} <= nu,  R >= lam,  R >= mu }
  Claim: ord_t(chat^nu_{lam,mu}) < 0  <==>  E_S != {}   (when chi(eps) != 1).
"""
import sys, os, random, argparse
from fractions import Fraction as Fr
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code-1017c2-coboundary'))
sys.path.insert(0, os.path.dirname(__file__))
import q484_engine as E
from ls import LS, K, solve
from fractions import Fraction

t = sp.Symbol('t')


def to_LS(expr, prec=K):
    expr = sp.cancel(sp.together(expr))
    num, den = sp.fraction(expr)
    pn = sp.Poly(sp.expand(num), t)
    pd = sp.Poly(sp.expand(den), t)
    def rat(c):
        a, b = sp.cancel(c).as_numer_denom()
        return Fr(int(a), int(b))
    cn = [rat(c) for c in reversed(pn.all_coeffs())]
    cd = [rat(c) for c in reversed(pd.all_coeffs())]
    n = LS(0, (cn + [Fr(0)] * prec)[:prec])
    d = LS(0, (cd + [Fr(0)] * prec)[:prec])
    if n.is_zero():
        return LS.zero()
    return n / d


def _ordt(expr):
    expr = sp.cancel(sp.together(expr))
    if expr == 0:
        return None
    num, den = sp.fraction(expr)
    def lo(p):
        c = list(reversed(sp.Poly(sp.expand(p), t).all_coeffs()))
        for i, x in enumerate(c):
            if sp.simplify(x) != 0:
                return i
        return None
    a, b = lo(num), lo(den)
    return None if a is None else a - b


def geq(M, I):
    return E.geq(M, I)


def E_S(lam, mu, nu, m, N):
    r"""{ R subset [N]\{m} : R u {m} <= nu, R >= lam, R >= mu }"""
    out = []
    for R in E.subsets(N):
        if m in R:
            continue
        Mset = tuple(sorted(R + (m,)))
        if geq(nu, Mset) and geq(R, lam) and geq(R, mu):
            out.append(R)
    return out


def run(N, X, m, eps, chis, q0, seed, prec=K, verbose=False):
    random.seed(seed)
    h = 2 if X == 'C' else 1
    assert eps ** h == 1, (eps, h)
    qv = sp.Rational(q0)
    subs0 = {E.q: qv}
    # GENERICITY GATE.  prop:short's hypothesis is that, at generic q, the only
    # positive-root factor of Atilde_{I,I} vanishing on {w_m = eps} is the j=k term
    # (1 - w_m^{-h}).  A random point can violate it (seed 1 at N=3 gave w_3 = 1/q,
    # which kills (1-q w_m w_3) and (1-q w_3/w_m) at w_m=+1 and inflated
    # ord_t(det Phi) from 4 to 8).  So we SAMPLE UNTIL the independent root product
    # has v_S(Atilde_{I,I}) = [m in I] for every I, and assert it.
    PP = E.subsets(N)
    for _attempt in range(200):
        wsym = []
        for j in range(1, N + 1):
            if j == m:
                wsym.append(sp.Integer(eps) + t)
            else:
                wsym.append(sp.Rational(random.randint(11, 400), random.randint(1, 23)))
        ok = True
        for I in PP:
            o = _ordt(E.diag_rootproduct(N, wsym, I).subs(subs0))
            if o != (1 if m in I else 0):
                ok = False
                break
        if ok:
            break
    else:
        raise RuntimeError("no generic point found")

    P = E.subsets(N)
    n = len(P)
    idx = {I: k for k, I in enumerate(P)}

    # Phi_{M,I} = Atilde_{M,I}
    Phi = [[None] * n for _ in range(n)]
    for M in P:
        for I in P:
            a = E.Atilde(N, wsym, M, I, X=X, simplify=True)
            Phi[idx[M]][idx[I]] = to_LS(a.subs(subs0), prec)

    # order of det, as an instrument control: must equal #{I : m in I} = 2^(N-1)
    _, ordet = solve(Phi, [LS.const(0)] * n, n)
    expect = sum(1 for I in P if m in I)
    ctrl_det = (ordet == expect)
    assert ctrl_det, f"det control failed: ord={ordet} expected={expect}"

    psi = lambda v: (1 - qv * v ** h) / (1 - qv)
    CHI = {'id': lambda v: sp.Integer(1),
           'psi': psi,
           'xpsi': lambda v: v * psi(v),
           'x': lambda v: v,
           # kappa = chi(eps) = eps + 2:  3 at eps=+1, 1 at eps=-1.  The complementary
           # arm to 'xpsi' (kappa = eps): together they give a positive AND a negative
           # arm at BOTH components of {w_m^2=1}, and a positive arm in type B, where
           # eps = 1 is the only admissible value and psi, v*psi all have chi(1) = 1.
           'p2': lambda v: v + 2}

    results = {}
    for name in chis:
        chi = CHI[name]
        Psi = {}
        for M in P:
            p = sp.Integer(1)
            for mm in M:
                p *= chi(wsym[mm - 1])
            Psi[M] = to_LS(p, prec)
        kappa = sp.simplify(chi(sp.Integer(eps)))
        rows = []
        for lam in P:
            for mu in P:
                f = []
                for M in P:
                    pl = Phi[idx[M]][idx[lam]]
                    pm = Phi[idx[M]][idx[mu]]
                    f.append(Psi[M] * pl * pm)
                x, _ = solve(Phi, f, n)
                for nu in P:
                    v = x[idx[nu]]
                    o = None if v.is_zero() else v.o
                    rows.append((lam, mu, nu, o))
        results[name] = (kappa, rows)
    return P, results, ctrl_det, ordet, expect


def report(N, X, m, eps, chis, q0, seed):
    P, results, ctrl_det, ordet, expect = run(N, X, m, eps, chis, q0, seed)
    print(f"### N={N} X={X} m={m} eps={eps} q0={q0} seed={seed}")
    print(f"    control ord_t(det Phi) = {ordet}, expected #{{I: m in I}} = {expect} -> "
          f"{'OK' if ctrl_det else 'MISMATCH'}")
    for name, (kappa, rows) in results.items():
        tp = fp = tn = fn = 0
        bad = []
        orders = {}
        for lam, mu, nu, o in rows:
            pole = (o is not None and o < 0)
            pred = (kappa != 1) and len(E_S(lam, mu, nu, m, N)) > 0
            orders[o] = orders.get(o, 0) + 1
            if pole and pred: tp += 1
            elif pole and not pred: fp += 1; bad.append(('FP', lam, mu, nu, o))
            elif (not pole) and pred: fn += 1; bad.append(('FN', lam, mu, nu, o))
            else: tn += 1
        print(f"  chi={name:5s} chi(eps)={kappa}  poles={tp+fp:4d}/{len(rows)}  "
              f"TP={tp:4d} FP={fp:4d} FN={fn:4d} TN={tn:4d}  orders={dict(sorted(orders.items(), key=lambda z:(z[0] is None, z[0])))}")
        for b in bad[:12]:
            print("      ", b)
        if len(bad) > 12:
            print(f"       ... {len(bad)-12} more")
    return results


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=2)
    ap.add_argument('--X', default='C')
    ap.add_argument('--m', type=int, default=1)
    ap.add_argument('--eps', type=int, default=-1)
    ap.add_argument('--q0', default='3/7')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--chis', default='id,psi,xpsi,x')
    a = ap.parse_args()
    report(a.N, a.X, a.m, a.eps, a.chis.split(','), a.q0, a.seed)
