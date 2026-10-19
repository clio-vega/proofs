r"""GATE 0b, first-hand.  Motegi arXiv:2610.15344, thm:BC-all-cardinality-biorthogonality
(L4365-4395), transcribed from the e-print source:

  sum_{I subset [N]} q^{|lam(I)|-|lam(L)|} * Atil_{M,I}(w,q)/Delta^X_M(w,q)
                                          * Atil_{Lc,Ic}(w,1/q)/Delta^X_{Lc}(w,1/q)
    = delta_{M,L}

with Delta^X_A(w,q) := WXtil_{N,N,|A|}(w_A|w|A) = Atil_{A,A}  (eq allcard-Delta-factor,
L4257-4262) and |lam(A)| = sum of the elements of A (L4157-4160).

CHECKED AGAINST MY OWN ENGINE, i.e. against Atilde built from the symmetrisation
formula -- a different mechanism from Motegi's branching/RTV proof.  This is the
gate the brief demands before the inverse formula is used.
"""
import sys, os, random
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code-1017c2-coboundary'))
import q484_engine as E


def lam(A):
    return sum(A)


def check(N, X, q0, ws, verbose=False):
    P = E.subsets(N)
    sub = {E.q: sp.Rational(q0)}
    subi = {E.q: 1 / sp.Rational(q0)}
    A = {}
    Ai = {}
    for M in P:
        for I in P:
            A[(M, I)] = sp.Rational(sp.cancel(E.Atilde(N, ws, M, I, X=X).subs(sub)))
            Ai[(M, I)] = sp.Rational(sp.cancel(E.Atilde(N, ws, M, I, X=X).subs(subi)))
    comp = lambda S: tuple(j for j in range(1, N + 1) if j not in S)
    ok = bad = 0
    for M in P:
        for L in P:
            tot = sp.Integer(0)
            for I in P:
                tot += sp.Rational(q0) ** (lam(I) - lam(L)) * \
                       (A[(M, I)] / A[(M, M)]) * \
                       (Ai[(comp(L), comp(I))] / Ai[(comp(L), comp(L))])
            tot = sp.Rational(sp.cancel(tot))
            want = 1 if M == L else 0
            if tot == want:
                ok += 1
            else:
                bad += 1
                print(f"   MISMATCH N={N} X={X} M={M} L={L}: got {tot} want {want}")
    return ok, bad


if __name__ == '__main__':
    tot_ok = tot_bad = 0
    for N in (1, 2, 3):
        for X in ('C', 'B'):
            for trial, seed in enumerate((7, 23)):
                random.seed(seed + 100 * N)
                ws = [sp.Rational(random.randint(11, 400), random.randint(1, 23)) for _ in range(N)]
                q0 = sp.Rational(random.randint(2, 50), random.randint(51, 97))
                o, b = check(N, X, q0, ws)
                tot_ok += o; tot_bad += b
                print(f"N={N} X={X} trial={trial} q={q0} w={ws}: {o} ok, {b} bad "
                      f"(of {4**N} (M,L) pairs)")
    print(f"TOTAL: {tot_ok} ok, {tot_bad} bad")
