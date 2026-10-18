"""(Q)-divisor structure, RATIONAL parametrisation.

t5 used w_1 = sqrt(1/q) + t for the divisor {q w_1^2 = 1} and reported
ord(det Phi) = 0; t3 used 1/sqrt(q) + t on the SAME object and reported 2.
Two instruments, one object, opposite answers -- the radical defeats the
polynomial-division order test (sp.expand does not recognise sqrt(q)
cancellations).  Fix: parametrise the divisor by eliminating q, which is
rational.  On {1 - q*m(w) = 0} the function field is C(w_1..w_N) with
q = 1/m(w); so substitute q -> (1-t)/m(w).  Then t IS the vanishing factor
1 - q*m(w), everything is rational, and no radical appears.
"""
import sympy as sp, itertools, sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1018c2-gauge')
from q484_engine import q, wvars, Atilde

t = sp.symbols('t')
def subsets(N):
    return [tuple(s) for k in range(N+1) for s in itertools.combinations(range(1,N+1),k)]

def ordt(expr):
    e = sp.cancel(sp.together(sp.expand(expr)))
    if e == 0: return None
    num, den = sp.fraction(e)
    def m(p):
        k = 0
        while True:
            P = sp.Poly(sp.expand(p), t)
            if P.total_degree() == 0 and P.as_expr() != 0: return k
            qq, rr = sp.div(P, sp.Poly(t, t))
            if sp.simplify(rr.as_expr()) == 0: p = qq.as_expr(); k += 1
            else: return k
    return m(num) - m(den)

def run(N, name, mono, SD):
    """divisor {1 - q*mono = 0}; substitute q -> (1-t)/mono."""
    w = wvars(N)
    sub = {q: (1 - t)/mono}
    tab = {}
    for M in subsets(N):
        for I in subsets(N):
            tab[(M,I)] = ordt(Atilde(N, w, M, I, X='C').subs(sub))
    allo = [o for o in tab.values() if o is not None]
    detord = sum(tab[(I,I)] for I in subsets(N) if tab[(I,I)] is not None)
    qv_ok = qv_n = 0
    for M in subsets(N):
        if not SD(M): continue
        for I in subsets(N):
            o = tab[(M,I)]
            if o is None: continue
            qv_n += 1; qv_ok += int(o >= 1)
    sd = [M for M in subsets(N) if SD(M)]
    extra = sorted({(M,I) for (M,I),o in tab.items() if o is not None and o >= 1 and not SD(M)})
    print(f"\n### N={N}  D = {{1 - q*{name} = 0}}")
    print(f"  (QR) regular       : {sum(1 for o in allo if o>=0)}/{len(allo)}  min ord = {min(allo)}")
    print(f"  (QV) vanish on S_D : {qv_ok}/{qv_n}   |S_D| = {len(sd)}")
    print(f"  ord(det Phi) = {detord}   |S_D| = {len(sd)}   MATCH = {detord == len(sd)}")
    print(f"  ord>=1 with M outside S_D : {len(extra)}" + (f"  {extra[:5]}" if extra else ""))
    return detord == len(sd) and qv_ok == qv_n and min(allo) >= 0

if __name__ == '__main__':
    ok = []
    for N in (2,3):
        w = wvars(N)
        ok.append(run(N, "w1*w2",  w[0]*w[1], lambda M: 1 not in M and 2 not in M))
        ok.append(run(N, "w1^2",   w[0]**2,   lambda M: 1 not in M))
        ok.append(run(N, "w2/w1",  w[1]/w[0], lambda M: 1 in M and 2 not in M))
    # a CONTROL: a divisor that is NOT in the diagonal at all -> colength must be 0
    N=2; w=wvars(N)
    ok.append(run(N, "w1^3  [NOT a diagonal factor: CONTROL, expect |S_D|=0=colength]",
                  w[0]**3, lambda M: False))
    print(f"\n==> all structural claims hold: {all(ok)}")
