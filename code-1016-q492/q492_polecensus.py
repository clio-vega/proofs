"""Q492/Q484: census the DENOMINATOR of each structure constant c^nu_{lam,mu}.

For every triple permitted by the support bound, factor the denominator of
c^nu (after sp.cancel) and record which irreducible factors occur, with
multiplicity.  Goal: identify the mechanism behind the measured criterion
   c^nu Laurent  <=>  |nu| - max(|lam|,|mu|) <= 1.
Pure SymPy.  Engine = q484_engine (transcribed from 2610.15344).
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
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

def den_factors(f, w):
    """return Counter of non-monomial irreducible factors of the denominator."""
    f = sp.cancel(sp.together(f))
    num,den = sp.fraction(f)
    out = Counter()
    for fac,mult in sp.factor_list(sp.expand(den))[1]:
        p = sp.Poly(fac,*w)
        if len(p.monoms())==1:   # monomial in w (incl. constants/pure q)
            continue
        out[sp.factor(fac)] += mult
    return out

def main(N, X='C'):
    t0=time.time(); w,S,A = build(N,X)
    print(f'=== N={N} type {X}   (build {time.time()-t0:.0f}s) ===', flush=True)
    # first: the diagonal entries' factor content
    print('--- diagonal entries A_{I,I}: long-root factors (w_a w_b - 1) ---')
    for I in S:
        d = sp.factor(A[(I,I)])
        n_lr = 0
        for fac,m in sp.factor_list(sp.expand(sp.fraction(sp.together(d))[0]))[1]:
            pass
        print(f'  I={I}  |I|={len(I)}   A_II = {sp.factor(A[(I,I)])}', flush=True)
    rows=[]
    for lam in S:
        for mu in S:
            if mu < lam: continue          # c is symmetric in lam<->mu
            c = solve_c(S,A,lam,mu)
            for nu in S:
                if c[nu]==0: continue
                jump = len(nu)-max(len(lam),len(mu))
                df = den_factors(c[nu],w)
                rows.append((lam,mu,nu,jump,df))
    print(f'\n--- {len(rows)} nonzero constants; jump vs denominator ---', flush=True)
    bad=[]
    for lam,mu,nu,jump,df in rows:
        if df:
            bad.append((lam,mu,nu,jump,df))
    print(f'  Laurent (empty denominator): {len(rows)-len(bad)}/{len(rows)}')
    print(f'  non-Laurent               : {len(bad)}/{len(rows)}')
    print('\n  jump distribution, Laurent vs not:')
    cl=Counter(j for _,_,_,j,df in rows if not df); cb=Counter(j for _,_,_,j,df in rows if df)
    for j in sorted(set(cl)|set(cb)):
        print(f'    jump {j}:  Laurent {cl.get(j,0):4d}   non-Laurent {cb.get(j,0):4d}')
    print('\n  --- the non-Laurent ones, with their pole divisors ---')
    for lam,mu,nu,jump,df in bad:
        s=' * '.join(f'({sp.sstr(f)})^{m}' if m>1 else f'({sp.sstr(f)})' for f,m in df.items())
        print(f'    lam={lam} mu={mu} nu={nu}  jump={jump}  den: {s}', flush=True)
    return rows

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 2, sys.argv[2] if len(sys.argv)>2 else 'C')
