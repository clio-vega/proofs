r"""Direct test of the ONLY remaining gap in Theorem S: non-cancellation.

Theorem S gives
    Res_nu(Psi phi_lam phi_mu) = (1-kappa) * sum_{R in E_S} gamma^nu_R * bar(Psi(R) phi_lam(R) phi_mu(R)),
every term individually nonzero.  The gap is whether a sum of >= 2 terms can vanish.
This script evaluates that sum directly in the residue field.

PROVENANCE, stated: this uses the residue formula, which is NOT verified here --
it is verified in struct_checks.py C3 (448/448 against a numerically inverted Phi)
and end-to-end in census_S.py at N=2,3 (direct linear solve, 0 FP / 0 FN).  So this
is an independent test of the NON-CANCELLATION step, conditional on that formula.

gamma^nu_R = - gamma^nu_{R u {m}} and, with Mp = R u {m},
  gamma^nu_{Mp} = Res_{t=0}[ q^{|lam nu|-|lam Mp|} Atil_{Mp^c,nu^c}(q^{-1})
                             / ( Atil_{Mp^c,Mp^c}(q^{-1}) Atil_{Mp,Mp}(q) ) ].

DOMAIN: N as given; X=C; every m; eps in {+1,-1}; ALL triples (lam,mu,nu);
twists chi = v, v*psi(v), v+2, and a random-interpolated chi; gated-generic point.
"""
import sys, os, random, argparse
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code-1017c2-coboundary'))
sys.path.insert(0, os.path.dirname(__file__))
import q484_engine as E
from struct_checks import generic_point, ordt, atS, comp


def run(N, m, eps, q0, chis, seed, X='C'):
    w = generic_point(N, m, eps, q0, seed)
    qv = sp.Rational(q0)
    sub, subi = {E.q: qv}, {E.q: 1 / qv}
    P = E.subsets(N)
    h = 2 if X == 'C' else 1

    # Atilde at q and at 1/q, as functions of t
    A = {(M, I): sp.cancel(E.Atilde(N, w, M, I, X=X).subs(sub)) for M in P for I in P}
    Ai = {(M, I): sp.cancel(E.Atilde(N, w, M, I, X=X).subs(subi)) for M in P for I in P}
    Abar = {k: atS(v) for k, v in A.items()}       # restriction to S (all regular)

    # gamma^nu_R  for m not in R
    gam = {}
    for nu in P:
        for R in P:
            if m in R:
                continue
            Mp = tuple(sorted(R + (m,)))
            Mpc, nuc = comp(Mp, N), comp(nu, N)
            expr = qv ** (sum(nu) - sum(Mp)) * Ai[(Mpc, nuc)] / (Ai[(Mpc, Mpc)] * A[(Mp, Mp)])
            o = ordt(expr)
            if o is None:
                gam[(nu, R)] = sp.Integer(0)
            else:
                assert o >= -1, (nu, R, o)
                gam[(nu, R)] = sp.Integer(0) if o >= 0 else atS(sp.cancel(sp.Symbol('t') * expr))
            gam[(nu, R)] = -gam[(nu, R)]           # gamma_R = -gamma_{R u m}

    psi = lambda v: (1 - qv * v ** h) / (1 - qv)
    random.seed(seed + 77)
    cvals = {j: sp.Rational(random.randint(2, 97), random.randint(1, 13))
             for j in range(1, N + 1) if j != m}
    CHI = {'x': lambda v, j: v,
           'xpsi': lambda v, j: v * psi(v),
           'p2': lambda v, j: v + 2,
           # interpolated: chi(w_j) = c_j independent values, chi(eps) = 2
           'interp': lambda v, j: cvals[j]}
    KAP = {'x': sp.Integer(eps), 'xpsi': sp.Integer(eps), 'p2': sp.Integer(eps + 2),
           'interp': sp.Integer(2)}

    out = {}
    for name in chis:
        chi, kappa = CHI[name], KAP[name]
        if kappa == 1:
            out[name] = ('kappa=1 (no pole predicted; skipped)', 0, 0, 0)
            continue
        Psi = {}
        for R in P:
            if m in R: continue
            p = sp.Integer(1)
            for r in R:
                p *= chi(w[r - 1], r)
            Psi[R] = atS(sp.cancel(p))
        tp = fp = fn = tn = 0
        bad = []
        for lam in P:
            for mu in P:
                for nu in P:
                    ES = [R for R in P if m not in R
                          and E.geq(nu, tuple(sorted(R + (m,))))
                          and E.geq(R, lam) and E.geq(R, mu)]
                    s = sp.Integer(0)
                    for R in ES:
                        s += gam[(nu, R)] * Psi[R] * Abar[(R, lam)] * Abar[(R, mu)]
                    s = sp.cancel((1 - kappa) * s)
                    pole = (s != 0)
                    pred = len(ES) > 0
                    if pole and pred: tp += 1
                    elif pole and not pred: fp += 1; bad.append(('FP', lam, mu, nu))
                    elif pred and not pole: fn += 1; bad.append(('FN(cancellation!)', lam, mu, nu, len(ES)))
                    else: tn += 1
        out[name] = (f"kappa={kappa}", tp, fp, fn)
        if bad:
            for b in bad[:8]:
                print("      ", b)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=4)
    ap.add_argument('--q0', default='3/7')
    ap.add_argument('--chis', default='x,xpsi,p2,interp')
    ap.add_argument('--seed', type=int, default=3000)
    a = ap.parse_args()
    for m in range(1, a.N + 1):
        for eps in (1, -1):
            r = run(a.N, m, eps, a.q0, a.chis.split(','), a.seed + m)
            for name, v in r.items():
                print(f"N={a.N} m={m} eps={eps:2d} chi={name:7s} {v[0]:14s} "
                      f"TP={v[1]:5d} FP={v[2]} FN={v[3]}")
            sys.stdout.flush()
