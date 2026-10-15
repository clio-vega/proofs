"""Control 3: Motegi's cross-sector branching, Thm L2886
(eq:typeC-cross-sector-branching), used as a CONSISTENCY control with real domain.

The theorem writes the s>r specialization as an explicit sum over K subset M,
|K|=r, of the type-A weight function W_{N,r}(w_K|w|I) with closed prefactors.
At s=r the sum has one term (Motegi's Remark, L2924-2929), which lets me SOLVE
for W_{N,r}(w_K|w|I) out of my own s=r engine values.

So:  LHS = engine at n=s  (all seven blocks of (Motherfunction) active, including
           the j>r blocks  (u_j u_k - q), (beta u_j - gamma), prod_k (u_j - q w_k))
     RHS = engine at n=r  (j>r blocks EMPTY) + Motegi's explicit prefactors.

The j>r blocks therefore occur on the LHS only.  A plant in them must fire here.
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp
from itertools import combinations
import q484_engine as E
q = E.q

def typeA_W(N, w, K, I, Afun):
    """Invert the s=r case of the branching theorem for W_{N,r}(w_K|w|I)."""
    r = len(K); Kc = [j for j in range(1,N+1) if j not in K]
    num = sp.Integer(1)
    for a in range(r):
        for b in range(r):
            num *= (1 - q*w[K[b]-1]/w[K[a]-1])
    den = sp.Integer(1)
    for a in range(len(Kc)):
        for b in range(a,len(Kc)):
            den *= (1 - q*w[Kc[a]-1]*w[Kc[b]-1])
    for a in range(r):
        for b in range(a,r):
            den *= (1 - 1/(w[K[a]-1]*w[K[b]-1]))
    return Afun(N,w,K,I) * num / den

def branch_rhs(N, w, M, I, Afun):
    s, r = len(M), len(I)
    Mc = [j for j in range(1,N+1) if j not in M]
    pre = (1-q)**(s-r)
    for a in range(len(Mc)):
        for b in range(a,len(Mc)):
            pre *= (1 - q*w[Mc[a]-1]*w[Mc[b]-1])
    tot = sp.Integer(0)
    for K in combinations(M, r):
        L = [m for m in M if m not in K]
        term = typeA_W(N,w,K,I,Afun)
        for a in range(r):
            for b in range(r):
                term /= (1 - q*w[K[b]-1]/w[K[a]-1])
        for a in range(r):
            for b in range(a,r):
                term *= (1 - 1/(w[K[a]-1]*w[K[b]-1]))
        for a in range(r):
            for b in range(s-r):
                term *= (1 - 1/(w[K[a]-1]*w[L[b]-1])) / (1 - w[L[b]-1]/w[K[a]-1])
        for a in range(s-r):
            for b in range(a+1,s-r):
                term *= (1 - q/(w[L[a]-1]*w[L[b]-1]))
        for a in range(s-r):
            for j in Mc:
                term *= (1 - q*w[j-1]/w[L[a]-1])
        tot += term
    return pre*tot

def run(N, Afun, label):
    w = E.wvars(N); ok = bad = 0; dom = 0
    for M in E.subsets(N):
        for I in E.subsets(N):
            if len(I) > len(M): continue
            if len(M) == len(I): continue          # s=r is the inverted case, vacuous
            dom += 1
            lhs = Afun(N,w,M,I)
            rhs = branch_rhs(N,w,M,I,Afun)
            if sp.simplify(sp.cancel(sp.together(lhs-rhs))) == 0: ok += 1
            else: bad += 1
    print(f'{label} N={N}: cross-sector branching {ok}/{dom} (s>r pairs only)'
          + ('' if bad==0 else '   -> FIRES'))
    return bad

if __name__ == '__main__':
    import q484_plants as P
    clean = lambda N,w,M,I: E.Atilde(N,w,M,I,'C')
    for N in (2,3):
        run(N, clean, 'CLEAN   ')
    for pl in ('P1_pairfactor','P2_drop_u2m1'):
        for N in (2,3):
            run(N, (lambda pl: lambda N,w,M,I: P.Atilde_planted(N,w,M,I,pl))(pl), f'PLANT {pl[:2]}')
