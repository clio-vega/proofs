r"""VERIFICATION of the proof of hyp:L^S by one-spectator degeneration.

CLAIM UNDER TEST (the mechanism of the proof, not a reformulation of it):
  Fix m, eps finite.  Put  w_m = eps  (the SPECTATOR) and  w_c = T^{-c} for
  c != m  (the HIERARCHY), and set q = 0.  Then for A subset [N]\{m}, B subset [N]:

        Atilde^X_{A,B}  ---> [A >= B]   as T -> 0+,   for X = B and X = C.

  The limit value 1 (not merely "nonzero") is the sharp prediction, and it is
  written down BEFORE the run (brief instrument fault 2).

DOMAIN: N = 1..5; X in {B,C}; every m in [N]; eps = +1 and (type C) -1;
        every pair (A,B) with m notin A.  Exhaustive, no sampling.
        T = 1/10, 1/100, 1/1000 -- exact rational arithmetic throughout.

PLANTS (two, both required to fire):
  P-tri : A not>= B  must give limit 0.
  P-ord : the SAME computation with a NON-monotone exponent vector e (so the
          hierarchy w_1 << ... << w_N is broken) must NOT converge to 1 for at
          least some A >= B -- otherwise the hierarchy is doing no work and the
          instrument is blind to it.
"""
import sys, itertools, argparse
import sympy as sp
sys.path.insert(0, '.')
import q484_engine as E

def setq0():
    E.q = sp.Integer(0)

def wvals(N, m, eps, T, e=None):
    """w_m = eps (spectator); w_c = T^{-e_c} for c != m."""
    if e is None:
        e = {c: c for c in range(1, N + 1)}
    return [sp.Integer(eps) if c == m else sp.Rational(1, 1) / T**e[c]
            for c in range(1, N + 1)]

def run(maxN, Ts, e_override=None, label='main'):
    setq0()
    bad = []
    nA = nB = 0
    for N in range(1, maxN + 1):
        P = E.subsets(N)
        for X in ('C', 'B'):
            epss = (1, -1) if X == 'C' else (1,)
            for m in range(1, N + 1):
                for eps in epss:
                    e = e_override(N) if e_override else None
                    for A in P:
                        if m in A: continue
                        for B in P:
                            pred = sp.Integer(1) if E.geq(A, B) else sp.Integer(0)
                            vals = []
                            for T in Ts:
                                w = wvals(N, m, eps, T, e)
                                v = E.Atilde(N, w, A, B, X=X, simplify=False)
                                vals.append(sp.Rational(sp.cancel(v)))
                            # Convergence test on the DECAY EXPONENT, not on a
                            # threshold: err(T) ~ c*T^p, so p = log(e1/e2)/log(T1/T2).
                            # The claim "limit = pred" is err -> 0, i.e. p > 0.
                            # (A fixed threshold 1e-3 manufactured 11 false failures
                            # in the first run: the true decay is exactly O(T).)
                            errs = [abs(v - pred) for v in vals]
                            if all(x == 0 for x in errs):
                                ok = True
                            else:
                                import math
                                ps = []
                                for i in range(len(errs) - 1):
                                    if errs[i] == 0:
                                        ps.append(float('inf')); continue
                                    if errs[i + 1] == 0:
                                        ps.append(float('inf')); continue
                                    ps.append(math.log(float(errs[i] / errs[i + 1]))
                                              / math.log(float(Ts[i] / Ts[i + 1])))
                                ok = all(p > 0.5 for p in ps)
                            if E.geq(A, B): nA += 1
                            else: nB += 1
                            if not ok:
                                bad.append((N, X, m, eps, A, B, int(pred),
                                            [str(x) for x in errs]))
    return nA, nB, bad

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--maxN', type=int, default=4)
    ap.add_argument('--plant-ord', action='store_true')
    a = ap.parse_args()
    Ts = [sp.Rational(1, 10), sp.Rational(1, 100), sp.Rational(1, 1000)]
    if a.plant_ord:
        # P-ord: reverse the hierarchy (e decreasing in c) -- breaks w_1<<...<<w_N
        eo = lambda N: {c: N + 1 - c for c in range(1, N + 1)}
        nA, nB, bad = run(a.maxN, Ts, e_override=eo)
        print(f"PLANT P-ord (reversed hierarchy), N<={a.maxN}: "
              f"A>=B instances {nA}, A not>=B {nB}, DEVIATIONS {len(bad)}")
        print("  plant FIRES (deviations > 0)" if bad else "  plant SILENT -- instrument blind!")
        for b in bad[:6]: print("   ", b)
    else:
        nA, nB, bad = run(a.maxN, Ts)
        print(f"MAIN, N<={a.maxN}, both types, every m, every admissible eps, exhaustive (A,B):")
        print(f"  A>=B  (limit must be 1): {nA}")
        print(f"  A not>=B (plant P-tri, limit must be 0): {nB}")
        print(f"  TOTAL {nA+nB}   FAILURES {len(bad)}")
        for b in bad[:12]: print("   ", b)
