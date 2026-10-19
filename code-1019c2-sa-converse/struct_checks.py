r"""Three structural checks the proof of Theorem S rests on.

(C1) COMPLEMENTATION is an order anti-automorphism of Motegi's order:
        M <= nu   <==>   M^c >= nu^c.
     DOMAIN: all pairs (M,nu) of subsets of [N], N = 1..7 (exhaustive).

(C2) LEMMA L^S (sharp triangularity after boundary restriction):
        A >= B, m not in A   ==>   Atilde_{A,B}|_{w_m = eps}  != 0.
     DOMAIN: N = 1..4, X in {B,C}, every m, every admissible eps, every pair
     A >= B with m not in A, at a gated-generic rational (q,w) point.
     Plant: the same quantity for A NOT >= B must be 0 (triangularity, M3).

(C3) THE RESIDUE IDENTITY  gamma_{R u {m}} = - gamma_R, where
        gamma_M := Res_{w_m=eps} (Phi^{-1})_{nu,M}.
     Computed two ways: (i) from the numerically inverted Phi; (ii) from
     Motegi's explicit inverse
        (Phi^{-1})_{nu,M} = q^{|lam nu|-|lam M|} Atil_{M^c,nu^c}(q^{-1})
                            / ( Atil_{M^c,M^c}(q^{-1}) Atil_{M,M}(q) ).
     DOMAIN: N = 2,3; X = C; every m, eps; every nu.
"""
import sys, os, random, itertools
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code-1017c2-coboundary'))
sys.path.insert(0, os.path.dirname(__file__))
import q484_engine as E
t = sp.Symbol('t')
comp = lambda S, N: tuple(j for j in range(1, N + 1) if j not in S)


# ---------------------------------------------------------------- C1
def C1(maxN=7):
    ok = bad = 0
    for N in range(1, maxN + 1):
        P = E.subsets(N)
        for M in P:
            for nu in P:
                a = E.geq(nu, M)                      # M <= nu
                b = E.geq(comp(M, N), comp(nu, N))    # M^c >= nu^c
                if a == b:
                    ok += 1
                else:
                    bad += 1
                    if bad <= 5:
                        print(f"   C1 FAIL N={N} M={M} nu={nu}: M<=nu {a}, Mc>=nuc {b}")
    print(f"(C1) complementation anti-automorphism: {ok} ok, {bad} bad "
          f"(exhaustive over all (M,nu), N<=7)")
    return bad == 0


# ---------------------------------------------------------------- generic point
def generic_point(N, m, eps, q0, seed):
    random.seed(seed)
    sub = {E.q: sp.Rational(q0)}
    for _ in range(300):
        w = []
        for j in range(1, N + 1):
            w.append(sp.Integer(eps) + t if j == m else
                     sp.Rational(random.randint(11, 400), random.randint(1, 23)))
        if all(ordt(E.diag_rootproduct(N, w, I).subs(sub)) == (1 if m in I else 0)
               for I in E.subsets(N)):
            return w
    raise RuntimeError('no generic point')


def ordt(expr):
    expr = sp.cancel(sp.together(expr))
    if expr == 0:
        return None
    num, den = sp.fraction(expr)
    def lo(p):
        c = list(reversed(sp.Poly(sp.expand(p), t).all_coeffs()))
        for i, x in enumerate(c):
            if x != 0:
                return i
        return None
    a, b = lo(num), lo(den)
    return None if a is None else a - b


def atS(expr):
    """value at t=0 (expr regular there)"""
    return sp.cancel(sp.together(expr)).subs(t, 0)


# ---------------------------------------------------------------- C2
def C2(maxN=4, q0='3/7'):
    tot = {'ok': 0, 'bad': 0, 'plant_ok': 0, 'plant_bad': 0}
    for N in range(1, maxN + 1):
        for X in ('C', 'B'):
            epss = (1, -1) if X == 'C' else (1,)
            for m in range(1, N + 1):
                for eps in epss:
                    w = generic_point(N, m, eps, q0, 1000 + N)
                    sub = {E.q: sp.Rational(q0)}
                    for A in E.subsets(N):
                        if m in A:
                            continue
                        for B in E.subsets(N):
                            v = atS(E.Atilde(N, w, A, B, X=X).subs(sub))
                            if E.geq(A, B):
                                if v != 0: tot['ok'] += 1
                                else:
                                    tot['bad'] += 1
                                    print(f"   C2 FAIL N={N} X={X} m={m} eps={eps} A={A} B={B}")
                            else:
                                if v == 0: tot['plant_ok'] += 1
                                else:
                                    tot['plant_bad'] += 1
                                    print(f"   C2 PLANT FAIL N={N} X={X} m={m} A={A} B={B} -> {v}")
    print(f"(C2) Lemma L^S: {tot['ok']} nonzero / {tot['ok']+tot['bad']} with A>=B, m notin A"
          f"   [{tot['bad']} failures];  plant (A not>= B must vanish): "
          f"{tot['plant_ok']}/{tot['plant_ok']+tot['plant_bad']}")
    return tot


# ---------------------------------------------------------------- C3
def C3(Ns=(2, 3), q0='3/7'):
    nok = nbad = fok = fbad = 0
    for N in Ns:
        for m in range(1, N + 1):
            for eps in (1, -1):
                w = generic_point(N, m, eps, q0, 2000 + N)
                sub = {E.q: sp.Rational(q0)}
                qv = sp.Rational(q0)
                P = E.subsets(N)
                n = len(P)
                idx = {I: k for k, I in enumerate(P)}
                Phi = sp.Matrix(n, n, lambda i, j: sp.cancel(
                    E.Atilde(N, w, P[i], P[j], X='C').subs(sub)))
                Inv = Phi.inv()
                # explicit Motegi inverse
                def motegi_inv(nu, M):
                    Mc, nuc = comp(M, N), comp(nu, N)
                    num = E.Atilde(N, w, Mc, nuc, X='C').subs({E.q: 1 / qv})
                    d1 = E.Atilde(N, w, Mc, Mc, X='C').subs({E.q: 1 / qv})
                    d2 = E.Atilde(N, w, M, M, X='C').subs({E.q: qv})
                    return sp.cancel(qv ** (sum(nu) - sum(M)) * num / (d1 * d2))
                for nu in P:
                    gam = {}
                    for M in P:
                        e1 = sp.cancel(Inv[idx[nu], idx[M]])
                        e2 = motegi_inv(nu, M)
                        if sp.cancel(e1 - e2) == 0: fok += 1
                        else:
                            fbad += 1
                            print(f"   C3 inverse-formula FAIL N={N} m={m} eps={eps} nu={nu} M={M}")
                        gam[M] = atS(sp.cancel(t * e1))      # residue (poles simple)
                    for R in P:
                        if m in R: continue
                        Mp = tuple(sorted(R + (m,)))
                        if sp.cancel(gam[Mp] + gam[R]) == 0: nok += 1
                        else:
                            nbad += 1
                            print(f"   C3 gamma FAIL N={N} m={m} eps={eps} nu={nu} R={R}: "
                                  f"{gam[Mp]} + {gam[R]}")
    print(f"(C3) Motegi inverse formula vs numeric Phi^-1: {fok} ok, {fbad} bad;  "
          f"gamma_{{R u m}} = -gamma_R: {nok} ok, {nbad} bad")


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', 'C1'): C1()
    if which in ('all', 'C2'): C2()
    if which in ('all', 'C3'): C3()
