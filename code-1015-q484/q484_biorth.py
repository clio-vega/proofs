"""Control 3: Motegi's all-cardinality biorthogonality, Thm L4364
eq:BC-all-cardinality-biorthogonality (L4383-4392).

    sum_{I subset [N]} q^{|lam(I)|-|lam(L)|}
        * [ Atilde_{M,I}(w,q)       / Delta_M(w,q)     ]
        * [ Atilde_{L^c,I^c}(w,q)|_{q->1/q} / Delta_{L^c}(w,1/q) ]  =  delta_{M,L}

Unlike the diagonal control, this runs over pairs with |I| < |M|, where the
j>r blocks of eq (Motherfunction) -- (u_j u_k - q), (beta u_j - gamma),
prod_k (u_j - q w_k) -- are ACTIVE.  So a plant in those blocks can fire here.
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp
import q484_engine as E
q = E.q

def run(N, Afun, label):
    w = E.wvars(N); S = E.subsets(N)
    A = {(M,I): Afun(N,w,M,I) for M in S for I in S}
    Delta = {M: A[(M,M)] for M in S}
    comp = lambda I: tuple(j for j in range(1,N+1) if j not in I)
    ok = bad = 0; firstbad = None
    for M in S:
        for L in S:
            tot = sp.Integer(0)
            for I in S:
                a1 = A[(M,I)]/Delta[M]
                a2 = (A[(comp(L),comp(I))]/Delta[comp(L)]).subs(q, 1/q)
                tot += q**(E.lam_size(I)-E.lam_size(L)) * a1 * a2
            tot = sp.cancel(sp.together(tot))
            tgt = sp.Integer(1) if M==L else sp.Integer(0)
            if sp.simplify(tot-tgt)==0: ok += 1
            else:
                bad += 1
                if firstbad is None: firstbad = (M,L,sp.simplify(tot))
    print(f'{label}  N={N}: biorthogonality {ok}/{ok+bad}'
          + ('' if bad==0 else f'   first failure M={firstbad[0]} L={firstbad[1]} value={firstbad[2]}'))
    return bad

if __name__ == '__main__':
    import q484_plants as P
    clean = lambda N,w,M,I: E.Atilde(N,w,M,I,'C')
    p1    = lambda N,w,M,I: P.Atilde_planted(N,w,M,I,'P1_pairfactor')
    p3    = lambda N,w,M,I: P.Atilde_planted(N,w,M,I,'P2_drop_u2m1')
    for N in (1,2):
        run(N, clean, 'CLEAN   ')
        run(N, p1,    'PLANT P1')
        run(N, p3,    'PLANT P2')
