"""TASK 1 probe: the (Q)-divisors -- i.e. the divisors {1 - q*monomial = 0} of
det Phi, which is what Rick's "(Q)-factors of the diagonal" denotes.

Question 1 (Rick, UID 825 sec 5): are the structure constants of Ahat regular there?
Measured here, first-hand:
  (a) ord_D(Atil_{M,I}) for every M,I at every (Q)-divisor D  -- are the BASIS
      functions even regular?
  (b) ord_D(det Phi) = the colength budget.
  (c) does a collapse/congruence relation hold at D, and between which fixed points?
"""
import sympy as sp, itertools, sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1018c2-gauge')
from q484_engine import q, wvars, Atilde

def subsets(N):
    return [tuple(s) for k in range(N+1) for s in itertools.combinations(range(1,N+1),k)]

def ordat(expr, t):
    """order of vanishing of expr in the uniformiser t (expr rational in t)."""
    e = sp.cancel(sp.together(sp.simplify(expr)))
    if e == 0: return None
    num, den = sp.fraction(e)
    def mult(p):
        k = 0
        p = sp.expand(p)
        while True:
            qq, rr = sp.div(sp.Poly(p, t), sp.Poly(t, t))
            if sp.simplify(rr.as_expr()) == 0:
                p = qq.as_expr(); k += 1
            else:
                return k, p
    kn, _ = mult(num); kd, _ = mult(den)
    return kn - kd

def probe(N, dname, sub, t):
    """sub: dict of substitutions putting us on D, with t the uniformiser."""
    w = wvars(N)
    tbl = {}
    for M in subsets(N):
        for I in subsets(N):
            a = Atilde(N, w, M, I, X="C")
            a = sp.cancel(sp.together(a.subs(sub)))
            o = ordat(a, t)
            tbl[(M,I)] = o
    polemin = min([o for o in tbl.values() if o is not None])
    diagord = [tbl[(I,I)] for I in subsets(N)]
    detord = sum(o for o in diagord if o is not None)
    print(f"\n### N={N}  D = {dname}")
    print(f"  min ord over all (M,I) : {polemin}   ({'POLES' if polemin<0 else 'all regular'})")
    print(f"  ord(diagonal) by |I| : " + ", ".join(f"{I}:{tbl[(I,I)]}" for I in subsets(N)))
    print(f"  ord(det Phi) = {detord}   (= colength budget)")
    return tbl, w

if __name__ == '__main__':
    t = sp.symbols('t')
    N = 2
    w = wvars(N)
    # (Q1) {q w_1 w_2 = 1}: w_2 = 1/(q w_1) + t-deformation
    probe(N, "{q w1 w2 = 1}", {w[1]: 1/(q*w[0]) + t}, t)
    # (Q1') {q w_1^2 = 1}
    probe(N, "{q w1^2 = 1}", {w[0]: 1/sp.sqrt(q) + t}, t)
    # (Q2) {w_1 = q w_2}  (the q-shifted A-type divisor)
    probe(N, "{w1 = q w2}", {w[0]: q*w[1] + t}, t)
    # control: a root hyperplane, where we KNOW all Atil are regular (Lemma reg / SRC)
    probe(N, "CONTROL root hyperplane {w1 w2 = 1}", {w[1]: 1/w[0] + t}, t)
    probe(N, "CONTROL root hyperplane {w1 = w2}", {w[0]: w[1] + t}, t)
