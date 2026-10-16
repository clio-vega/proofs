"""TASK 0 GATE, sharpened.  The decisive columns are:
      ord_q(c)      -- q-adic valuation of the WHOLE structure constant
      ord_q(res)    -- q-adic valuation of its leading residue on H_ab
      c|_{q=0}      -- is it identically zero?
Rick's M2 ("at q=0 every structure constant is Laurent") and conj:divisor
("19 long-root poles at N=3") are reconciled iff every generic-pole triple has
ord_q(c) >= 1, i.e. the non-Laurent constants are the ones that VANISH at q=0.
Sharper still: ord_q(c) == ord_q(res) says the pole sits in the LEADING
q-coefficient, so normalising by q^{ord} keeps it.
"""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1016c2-nonvanishing')
import sympy as sp, time
from collections import Counter
import q484_engine as E
from gate_q0 import build, solve_c, nonmonomial_den, qval, longroot_label
q = E.q

def residue(c, w, a, b, k):
    s = sp.Symbol('s')
    cs = sp.cancel(c.subs(w[b-1], s/w[a-1]))
    lead = sp.cancel(sp.together((s-1)**k * cs))
    n,d = sp.fraction(lead)
    n1 = sp.expand(n).subs(s,1); d1 = sp.expand(d).subs(s,1)
    if d1 == 0: return None
    return sp.cancel(n1/d1)

def main(N, X='C'):
    t0=time.time(); w,S,A = build(N,X)
    print(f'=== GATE2  N={N} type {X}  (build {time.time()-t0:.0f}s) ===',flush=True)
    rows=[]; ntot=0
    for lam in S:
        for mu in S:
            if mu < lam: continue
            c = solve_c(S,A,lam,mu)
            for nu in S:
                if c[nu]==0: continue
                ntot+=1
                dfg = nonmonomial_den(c[nu],w)
                if not dfg: continue
                vc = qval(c[nu])
                for fac,m in dfg.items():
                    ab = longroot_label(fac,N,w)
                    if ab is None:
                        rows.append((lam,mu,nu,len(nu)-max(len(lam),len(mu)),
                                     sp.sstr(fac),m,vc,None)); continue
                    r = residue(c[nu],w,ab[0],ab[1],m)
                    rows.append((lam,mu,nu,len(nu)-max(len(lam),len(mu)),
                                 f'H_{ab[0]}{ab[1]}',m,vc,qval(r) if r is not None else None))
                    print(f'  lam={lam} mu={mu} nu={nu} jump={rows[-1][3]} '
                          f'{rows[-1][4]} mult={m}  ord_q(c)={vc}  ord_q(res)={rows[-1][7]}',
                          flush=True)
    print(f'\ntotal nonzero constants: {ntot}')
    print(f'(triple,divisor) pole records: {len(rows)}')
    print(f'distinct triples with a pole : {len(set((a,b,c_) for a,b,c_,*_ in rows))}')
    print(f'\nord_q(c) distribution over pole records: {dict(Counter(r[6] for r in rows))}')
    print(f'ord_q(res) distribution            : {dict(Counter(r[7] for r in rows))}')
    agree = sum(1 for r in rows if r[6]==r[7])
    print(f'\nord_q(c) == ord_q(res):  {agree}/{len(rows)}')
    zero_at_q0 = sum(1 for r in rows if r[6] and r[6]>=1)
    print(f'ord_q(c) >= 1 (so c|_(q=0) = 0):  {zero_at_q0}/{len(rows)}')
    jumps = Counter(r[3] for r in rows)
    print(f'jump distribution over pole records: {dict(jumps)}')

if __name__=='__main__':
    main(int(sys.argv[1]), sys.argv[2] if len(sys.argv)>2 else 'C')
