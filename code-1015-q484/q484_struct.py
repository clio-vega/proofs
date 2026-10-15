"""Q484 main computation: structure constants of the Motegi type-B/C weight
function basis in the localized (specialization-point) picture.

   Atilde_{M,lam} * Atilde_{M,mu} = sum_nu c^nu_{lam,mu} Atilde_{M,nu}   for all M subset [N]

Solved by forward substitution along a linear extension of Motegi's order >=
(Theorem B).  Then tested against:
  Theorem A  support bound      c^nu != 0  =>  nu >= lam and nu >= mu
  Cor C      extremal value     lam >= mu  =>  c^lam_{lam,mu} = Atilde_{lam,mu}
  Theorem D  localization form  c = A^{-1} v with Motegi's explicit inverse
  (G2)       Laurent polynomiality of c^nu in w
  (G1)       whether the identity LIFTS from the 2^N points to functions of u
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import sympy as sp, itertools, time
import q484_engine as E
q = E.q

def linext(S):
    return sorted(S, key=lambda M: (len(M), M))

def build(N, X='C'):
    w = E.wvars(N); S = linext(E.subsets(N))
    A = {}
    for M in S:
        for I in S:
            A[(M,I)] = E.Atilde(N,w,M,I,X)
    return w, S, A

def solve_c(S, A, lam, mu):
    """forward substitution: Theorem B."""
    c = {}
    for nu in S:                                   # S in linear-extension order
        v = A[(nu,lam)]*A[(nu,mu)]
        acc = sp.Integer(0)
        for kap in S:
            if kap == nu: break
            if kap in c and c[kap] != 0:
                acc += c[kap]*A[(nu,kap)]
        c[nu] = sp.cancel(sp.together((v-acc)/A[(nu,nu)]))
    return c

def main(N, X='C'):
    t0=time.time(); w,S,A = build(N,X)
    comp = lambda I: tuple(j for j in range(1,N+1) if j not in I)
    nA = nC = nD = nPoly = 0; totA = totC = totD = totPoly = 0
    viol = []
    allc = {}
    for lam in S:
        for mu in S:
            c = solve_c(S,A,lam,mu); allc[(lam,mu)] = c
            # Theorem A
            for nu in S:
                totA += 1
                nz = (c[nu] != 0)
                if nz and not (E.geq(nu,lam) and E.geq(nu,mu)):
                    viol.append(('A',lam,mu,nu)); 
                else: nA += 1
            # Cor C
            if E.geq(lam,mu):
                totC += 1
                if sp.simplify(sp.cancel(c[lam]-A[(lam,mu)]))==0: nC += 1
                else: viol.append(('C',lam,mu,lam))
            # Theorem D: localization formula with Motegi's explicit inverse
            for nu in S:
                totD += 1
                rhs = sp.Integer(0)
                for M in S:
                    Ainv = q**(E.lam_size(nu)-E.lam_size(M)) * \
                           (A[(comp(M),comp(nu))]/A[(comp(M),comp(M))]).subs(q,1/q)
                    rhs += Ainv * A[(M,lam)]*A[(M,mu)]/A[(M,M)]
                if sp.simplify(sp.cancel(sp.together(rhs-c[nu])))==0: nD += 1
                else: viol.append(('D',lam,mu,nu))
            # (G2) Laurent polynomiality in w  (q kept rational)
            for nu in S:
                if c[nu]==0: continue
                totPoly += 1
                f = sp.cancel(sp.together(c[nu]))
                num,den = sp.fraction(f)
                den = sp.expand(den)
                # Laurent in w  <=>  denominator is a monomial in the w's times a
                # polynomial in q alone
                dpoly = sp.Poly(den, *w)
                if len(dpoly.monoms())==1: nPoly += 1
                else: viol.append(('G2',lam,mu,nu))
    print(f'N={N} X={X}  [{time.time()-t0:.0f}s]')
    print(f'  Theorem A support bound : {nA}/{totA}')
    print(f'  Cor C extremal value    : {nC}/{totC}')
    print(f'  Theorem D localization  : {nD}/{totD}')
    print(f'  (G2) Laurent in w       : {nPoly}/{totPoly}')
    if viol:
        from collections import Counter
        print('  violations by kind:', Counter(v[0] for v in viol))
        for v in viol[:5]: print('    ', v)
    return w,S,A,allc

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 2)
