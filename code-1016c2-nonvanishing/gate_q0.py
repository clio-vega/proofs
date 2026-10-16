"""TASK 0 GATE.  Reproduce Rick's M2 measurement on MY engine, and -- the part
his measurement does not report -- the q-ADIC ORDER of the residue at each
long-root pole.

Three columns per triple (lam,mu,nu):
  (1) generic:  long-root factors of the denominator of c in QQ(q)(w).
  (2) q->0   :  is the LIMIT q->0 of c Laurent in w?   (Rick's M2)
  (3) resorder: for each generic long-root pole H_ab of multiplicity k, the
                q-adic valuation of  lim_{s->1} (s-1)^k c  where s = w_a w_b.

(1) and (2) are different questions: c lives in QQ(q)(w) where q is an
INDETERMINATE, so "has a pole" in (1) is already a generic-in-q statement and
specialisation q=0 need not commute with it.  (3) measures exactly the
discrepancy.
"""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1016c2-nonvanishing')
import sympy as sp, itertools, time
from collections import Counter
import q484_engine as E
q = E.q

def linext(S): return sorted(S, key=lambda M:(len(M),M))

def build(N, X='C'):
    w = E.wvars(N); S = linext(E.subsets(N))
    A = {(M,I): E.Atilde(N,w,M,I,X) for M in S for I in S}
    return w,S,A

def solve_c(S,A,lam,mu):
    c={}
    for nu in S:
        v=A[(nu,lam)]*A[(nu,mu)]
        acc=sp.Integer(0)
        for kap in S:
            if kap==nu: break
            if c.get(kap,0)!=0: acc+=c[kap]*A[(nu,kap)]
        c[nu]=sp.cancel(sp.together((v-acc)/A[(nu,nu)]))
    return c

def nonmonomial_den(f, w):
    """Counter of non-monomial irreducible denominator factors."""
    f = sp.cancel(sp.together(f)); num,den = sp.fraction(f)
    out = Counter()
    for fac,mult in sp.factor_list(sp.expand(den))[1]:
        p = sp.Poly(fac,*w)
        if len(p.monoms())==1: continue
        out[sp.factor(fac)] += mult
    return out

def qval(expr):
    """q-adic valuation of a rational function of q (and w): ord_{q=0}."""
    expr = sp.cancel(sp.together(expr))
    if expr == 0: return sp.oo
    num,den = sp.fraction(expr)
    def low(p):
        P = sp.Poly(sp.expand(p), q)
        return min(m[0] for m,c in zip(P.monoms(),P.coeffs()) if c != 0)
    return low(num)-low(den)

def limit_q0(f):
    """lim_{q->0} f, as a rational function in w.  Returns (value, qval)."""
    v = qval(f)
    if v == sp.oo: return sp.Integer(0), sp.oo
    g = sp.cancel(sp.together(sp.cancel(f)/q**v))
    num,den = sp.fraction(g)
    num0 = sp.expand(num).subs(q,0); den0 = sp.expand(den).subs(q,0)
    assert den0 != 0, 'den still vanishes at q=0'
    val0 = sp.cancel(num0/den0)
    return val0, v

def longroot_label(fac, N, w):
    """if fac is (up to unit/monomial) w_a w_b - 1 with a<b, return (a,b)."""
    for a in range(1,N+1):
        for b in range(a+1,N+1):
            if sp.simplify(sp.cancel(fac/(w[a-1]*w[b-1]-1))).is_number:
                return (a,b)
    return None

def residue_qval(c, N, w, a, b, k):
    """q-adic valuation of lim_{s->1}(s-1)^k c, s = w_a w_b."""
    s = sp.Symbol('s')
    # w_b = s / w_a
    cs = sp.cancel(c.subs(w[b-1], s/w[a-1]))
    lead = sp.cancel(sp.together((s-1)**k * cs))
    num,den = sp.fraction(lead)
    # evaluate at s=1
    num1 = sp.expand(num).subs(s,1); den1 = sp.expand(den).subs(s,1)
    if den1 == 0:
        return None, None   # higher-order structure; report separately
    lead1 = sp.cancel(num1/den1)
    return lead1, qval(lead1)

def main(N, X='C'):
    t0=time.time(); w,S,A = build(N,X)
    print(f'=== GATE  N={N} type {X}   (build {time.time()-t0:.0f}s) ===', flush=True)
    n_tot=0; n_gen_pole=0; n_q0_pole=0
    gen_poles=[]; q0_poles=[]; res_data=[]
    for lam in S:
        for mu in S:
            if mu < lam: continue
            c = solve_c(S,A,lam,mu)
            for nu in S:
                if c[nu]==0: continue
                n_tot+=1
                jump = len(nu)-max(len(lam),len(mu))
                dfg = nonmonomial_den(c[nu],w)
                c0,v = limit_q0(c[nu])
                df0 = nonmonomial_den(c0,w) if c0!=0 else Counter()
                if dfg: n_gen_pole+=1; gen_poles.append((lam,mu,nu,jump,dfg))
                if df0: n_q0_pole+=1; q0_poles.append((lam,mu,nu,jump,df0))
                for fac,m in dfg.items():
                    ab = longroot_label(fac,N,w)
                    if ab is None:
                        res_data.append((lam,mu,nu,jump,sp.sstr(fac),m,'not-longroot',None))
                        continue
                    a,b = ab
                    lead,rv = residue_qval(c[nu],N,w,a,b,m)
                    res_data.append((lam,mu,nu,jump,f'H_{a}{b}',m,rv,lead))
    print(f'\n[1] GENERIC (q indeterminate):  {n_tot-n_gen_pole}/{n_tot} Laurent,'
          f' {n_gen_pole} with a pole', flush=True)
    print(f'[2] q -> 0 LIMIT             :  {n_tot-n_q0_pole}/{n_tot} Laurent,'
          f' {n_q0_pole} with a pole', flush=True)
    print(f'\n--- generic poles ({len(gen_poles)} triples) ---')
    for lam,mu,nu,j,df in gen_poles:
        print(f'  lam={lam} mu={mu} nu={nu} jump={j}: '
              +' * '.join(f'({sp.sstr(f)})^{m}' for f,m in df.items()), flush=True)
    print(f'\n--- q=0 poles ({len(q0_poles)} triples) ---')
    for lam,mu,nu,j,df in q0_poles:
        print(f'  lam={lam} mu={mu} nu={nu} jump={j}: '
              +' * '.join(f'({sp.sstr(f)})^{m}' for f,m in df.items()), flush=True)
    print(f'\n--- [3] q-adic valuation of the leading residue, per (triple, divisor) ---')
    print(f'  {"lam":>8} {"mu":>8} {"nu":>8}  jmp  divisor  mult  ord_q(res)')
    vc=Counter()
    for lam,mu,nu,j,dv,m,rv,lead in res_data:
        print(f'  {str(lam):>8} {str(mu):>8} {str(nu):>8}   {j}   {dv:>6}   {m}     {rv}', flush=True)
        vc[rv]+=1
    print(f'\n  ord_q(residue) distribution: {dict(vc)}', flush=True)
    return res_data

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 2,
         sys.argv[2] if len(sys.argv)>2 else 'C')
