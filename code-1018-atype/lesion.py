"""STEP 5.  LESION THE THEOREM'S HYPOTHESES, one at a time.

Theorem B rests on exactly three things.  Each is lesioned here, and a lesion that
SURVIVES is a stronger theorem for free, while one that BREAKS confirms the
hypothesis is load-bearing rather than decorative.

L1  rho^A = 1  (Lemma A1).  Lesioned ALGEBRAICALLY, not numerically: suppose the
    multiplier were a constant unit c != 1, so G_c := {f : f(M) = c f(sM) mod t}.
    Then (fg)(M) = c^2 (fg)(sM), weight c^2 != c, so G_c is NOT closed under
    multiplication.  rho = 1 is necessary, and this is the rho^2 - rho mechanism.
    (No computation needed; recorded for completeness.)

L2  the COLENGTH match (Lemma colenL vs colenG).  THE SHARP LESION.  Multiply ONE
    basis column phi_{I0} by (w_1 - w_2).  This
      * keeps every entry REGULAR, so it does not touch Lemma A2;
      * is LABEL-indexed, not point-indexed, so it leaves the GKM congruence of the
        other columns alone and violates the point-product form instead;
      * raises sum_I ord Atilde_{I,I} from 2^(N-2) to 2^(N-2) + 1, so L is now a
        colength-1 SUBLATTICE of G.
    Prediction: A-type poles appear.  If they do NOT, the colength match is not what
    makes the theorem work and Theorem A's proof is the wrong proof.
    This is also the control for the NON-VANISHING half of admissibility, which the
    'atp' plant (a pole) does not probe.

L3  point-indexedness of Psi_chi.  Lesioned by a chi depending on TWO indices:
    Psi(M) = prod_{m in M} chi(w_m) replaced by a non-product weight that still has
    unit values.  Tests whether 'unit' alone suffices or 'unit AND point-product'
    is needed.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1018-atype')
import sympy as sp
import q484_engine as E, divisors as D
q = E.q

def census(N, X, label, Amap, S, w):
    nonzero=0; byclass={}
    for lam in S:
        for mu in S:
            if mu < lam: continue
            c={}
            for nu in S:
                v=Amap[(nu,lam)]*Amap[(nu,mu)]
                acc=sp.Integer(0)
                for kap in S:
                    if kap==nu: break
                    if c.get(kap,0)!=0: acc+=c[kap]*Amap[(nu,kap)]
                c[nu]=sp.cancel(sp.together((v-acc)/Amap[(nu,nu)]))
                if c[nu]==0: continue
                nonzero+=1
                for k,lst in D.den_classify(c[nu],w,q).items():
                    byclass.setdefault(k,set()).add((lam,mu,nu))
    allbad=set().union(*byclass.values()) if byclass else set()
    A=byclass.get('A',set())
    print(f'  [{label}] nonzero={nonzero} non-Laurent={len(allbad)} '
          f'A-type={len(A)}  '
          f'{{{", ".join(f"{k}:{len(v)}" for k,v in sorted(byclass.items()))}}}', flush=True)
    return A

def main(N=2, X='C'):
    w=E.wvars(N); S=sorted(E.subsets(N), key=lambda M:(len(M),M))
    A0={(M,I):sp.cancel(E.Atilde(N,w,M,I,X)) for M in S for I in S}
    print(f'=== N={N} type {X}: LESIONS ===', flush=True)
    print('  baseline (no lesion):', flush=True)
    base = census(N,X,'id',dict(A0),S,w)
    t=(w[0]-w[1])
    print(f'\n  L2: multiply ONE COLUMN phi_I0 by (w1-w2) -- regular, label-indexed,'
          f' raises the colength by 1', flush=True)
    for I0 in S:
        Am={k:(v*t if k[1]==I0 else v) for k,v in A0.items()}
        census(N,X,f'L2 I0={I0}',Am,S,w)
    print(f'\n  L3: a NON-product unit weight on the points (violates point-indexedness,'
          f' keeps unit values)', flush=True)
    # weight 1 on s_12-fixed points, and (1 + (w1-w2)) vs 1 on the moved pair:
    # unit everywhere, but NOT of the form prod_{m in M} f(w_m)
    for name, wt in [('1+(w1-w2) on M={1} only', lambda M: (1+t) if M==(1,) else sp.Integer(1)),
                     ('2 on M={1} only',          lambda M: sp.Integer(2) if M==(1,) else sp.Integer(1))]:
        Am={k:sp.cancel(v*wt(k[0])) for k,v in A0.items()}
        census(N,X,f'L3 {name}',Am,S,w)

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 2, sys.argv[2] if len(sys.argv)>2 else 'C')
