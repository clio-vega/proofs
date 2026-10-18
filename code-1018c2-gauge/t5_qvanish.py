"""(Q)-divisor structure, measured.  Full ord table + the two claims the argument needs.

Claim (QV) [vanishing]:  v_D(Atil_{M,I}) >= 1 for every M in S_D and every I, where
     D={q w_c w_d=1} (c<d)   S_D = {M : M n {c,d} = {}}
     D={q w_c^2  =1}         S_D = {M : c notin M}
     D={w_c = q w_e}, e<c    S_D = {M : e in M, c notin M}
Claim (QR) [regularity]: v_D(Atil_{M,I}) >= 0 for every M,I.
Then colength(span) = v_D(det Phi) = |S_D| = colength{f : f|_{S_D} = 0 mod t},
so span = that set, which is closed under products => no (Q)-divisor poles.
"""
import sympy as sp, itertools, sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1018c2-gauge')
from q484_engine import q, wvars, Atilde, WX, EX

t = sp.symbols('t')
def subsets(N):
    return [tuple(s) for k in range(N+1) for s in itertools.combinations(range(1,N+1),k)]

def ordt(expr):
    e = sp.cancel(sp.together(expr))
    if e == 0: return None
    num, den = sp.fraction(e)
    def m(p):
        k = 0; p = sp.expand(p)
        while True:
            P = sp.Poly(p, t)
            qq, rr = sp.div(P, sp.Poly(t, t))
            if sp.expand(rr.as_expr()) == 0: p = qq.as_expr(); k += 1
            else: return k
    return m(num) - m(den)

def run(N, name, sub, SD):
    w = wvars(N)
    tab = {}
    for M in subsets(N):
        for I in subsets(N):
            a = Atilde(N, w, M, I, X='C').subs(sub)
            tab[(M,I)] = ordt(a)
    allo = [o for o in tab.values() if o is not None]
    detord = sum(tab[(I,I)] for I in subsets(N) if tab[(I,I)] is not None)
    # claim QV
    qv_ok = qv_n = 0
    for M in subsets(N):
        if not SD(M): continue
        for I in subsets(N):
            o = tab[(M,I)]
            if o is None: continue        # identically zero: vacuously >=1
            qv_n += 1; qv_ok += int(o >= 1)
    # claim QR
    qr_n = len(allo); qr_ok = sum(1 for o in allo if o >= 0)
    sd = [M for M in subsets(N) if SD(M)]
    print(f"\n### N={N}  D = {name}")
    print(f"  (QR) regular        : {qr_ok}/{qr_n}   min ord = {min(allo)}")
    print(f"  (QV) vanish on S_D  : {qv_ok}/{qv_n}   |S_D| = {len(sd)}  S_D = {sd}")
    print(f"  ord(det Phi) = {detord}    |S_D| = {len(sd)}   MATCH = {detord == len(sd)}")
    # where are the strictly-positive orders, outside S_D?
    extra = sorted({(M,I) for (M,I),o in tab.items() if o is not None and o >= 1 and not SD(M)})
    print(f"  (M,I) with ord>=1 but M not in S_D : {len(extra)}" + (f"  e.g. {extra[:4]}" if extra else ""))
    return tab

if __name__ == '__main__':
    for N in (2,3):
        w = wvars(N)
        # (Q1) q w1 w2 = 1
        run(N, "{q w1 w2 = 1}", {w[1]: 1/(q*w[0]) + t}, lambda M: 1 not in M and 2 not in M)
        # (Q1') q w1^2 = 1
        run(N, "{q w1^2 = 1}", {w[0]: sp.sqrt(1/q) + t},  lambda M: 1 not in M)
        # (Q2) w2 = q w1  i.e. 1 - q w2/w1 ... careful: factor is (1-q w_c/w_{i_l}), i_l<c
        run(N, "{w2 = w1/q}  [1 - q w2/w1 = 0]", {w[1]: w[0]/q + t}, lambda M: 1 in M and 2 not in M)
