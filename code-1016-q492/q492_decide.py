"""THE DECIDER.  The jump criterion (Q484) and the uncovered criterion
(mechanism) both predict 176/195 Laurent at N=3 but disagree on 12 triples.
Compute exactly those 12 and see which criterion is right.
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp, time
import q484_engine as E
q=E.q

N=3
w=E.wvars(N); S=sorted(E.subsets(N), key=lambda M:(len(M),M))
t0=time.time()
A={(M,I):E.Atilde(N,w,M,I,'C') for M in S for I in S}
print(f'built A ({time.time()-t0:.0f}s)',flush=True)

def solve_c(lam,mu):
    c={}
    for nu in S:
        v=A[(nu,lam)]*A[(nu,mu)]
        acc=sp.Integer(0)
        for kap in S:
            if kap==nu: break
            if c.get(kap,0)!=0: acc+=c[kap]*A[(nu,kap)]
        c[nu]=sp.cancel(sp.together((v-acc)/A[(nu,nu)]))
    return c

def is_laurent(f):
    f=sp.cancel(sp.together(f))
    if f==0: return None
    den=sp.expand(sp.fraction(f)[1])
    return len(sp.Poly(den,*w).monoms())==1

def dendesc(f):
    f=sp.cancel(sp.together(f))
    den=sp.expand(sp.fraction(f)[1])
    out=[]
    for fac,m in sp.factor_list(den)[1]:
        if len(sp.Poly(fac,*w).monoms())>1:
            out.append(f'({sp.sstr(sp.factor(fac))})'+(f'^{m}' if m>1 else ''))
    return ' * '.join(out) if out else 'MONOMIAL (Laurent)'

cases=[((),(1,),(2,3)), ((),(2,),(1,3)),
       ((1,),(1,),(2,3)), ((2,),(2,),(1,3)),
       ((1,),(2,),(1,2,3)), ((1,),(3,),(1,2,3)), ((2,),(3,),(1,2,3))]
cache={}
print('\n lam        mu         nu          jumpJ  uncovU   MEASURED   denominator')
print(' '+'-'*92)
for lam,mu,nu in cases:
    key=(lam,mu)
    if key not in cache: cache[key]=solve_c(lam,mu)
    c=cache[key][nu]
    J = len(nu)-max(len(lam),len(mu))<=1
    U = len(set(nu)-(set(lam)|set(mu)))<=1
    L = is_laurent(c)
    verdict = 'Laurent' if L else ('ZERO' if L is None else 'POLE   ')
    agree = []
    if L is not None:
        agree.append('J-ok' if J==L else 'J-WRONG')
        agree.append('U-ok' if U==L else 'U-WRONG')
    print(f' {str(lam):10s} {str(mu):10s} {str(nu):11s} {str(J):5s}  {str(U):6s}  {verdict:9s}  {dendesc(c)}   {" ".join(agree)}',flush=True)
