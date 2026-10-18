"""STEP 1.  The A-TYPE GKM congruence and its multiplier rho^A.

Fixed points of LG(N,2N) <-> subsets M of [N].  T-curves <-> (point, reflection).
For the A-type root e_i - e_j the reflection s_ij SWAPS membership of i and j, so it
moves M exactly when |M cap {i,j}| = 1, and the curve weight is 1 - w_i/w_j, i.e.
(w_i - w_j) up to a unit.  [Cross-check of this model: at N=2 it predicts 6 curves
total over all four reflections -- which is exactly the count in Rick's LG(2,4) GKM
checker, UID 824.]

TEST.  For each i<j, each M with |M cap {i,j}| = 1, each I:
    Atilde_{M,I} |_{w_j = w_i}   =?=   rho^A * Atilde_{s_ij M, I} |_{w_j = w_i}
Strategy as in qc_fast.py: exact rational w with w_j = w_i imposed EXACTLY, q kept
symbolic, so each test is an identity of rational functions in q.

CONTROLS (each must FIRE, and I record WHICH FACTOR each is sensitive to):
  C1  off-divisor arm: impose w_j = (5/4) w_i instead.  Sensitive to: the collision
      itself.  Must fail where the on-divisor arm passes.
  C2  planted multiplier rho = 1 + q.  Sensitive to: the VALUE of the multiplier on
      the NONZERO cells only -- it is structurally spared on triangularity zeros
      (0 = rho*0 for any rho), so the two columns are reported separately.
  C3  planted WRONG PARTNER: compare against s_ik M for a third index k != i,j
      (N >= 3 only).  Sensitive to: which point the congruence pairs M with.
      A green here would mean the test cannot see the pairing at all.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1018-atype')
import sympy as sp
from itertools import combinations
import q484_engine as E
q = E.q

SEEDS = [[sp.Rational(3), sp.Rational(5), sp.Rational(7), sp.Rational(11)],
         [sp.Rational(2), sp.Rational(1,3), sp.Rational(5,2), sp.Rational(7,4)]]

def swap(M, i, j):
    """s_ij acting on the subset M: swap membership of i and j."""
    s = set(M)
    if (i in s) != (j in s):
        s.symmetric_difference_update({i, j})
    return tuple(sorted(s))

def run(N, seed, X='C', rho=None, off=False, wrong_partner=False):
    base = seed[:N]; S = sorted(E.subsets(N), key=lambda M:(len(M),M))
    rho = rho if rho is not None else sp.Integer(1)
    ok = bad = triv = 0; ratios = {}
    for i, j in combinations(range(1, N+1), 2):
        for M in S:
            if len({i,j} & set(M)) != 1: continue
            Ms = swap(M, i, j)
            if wrong_partner:
                ks = [k for k in range(1,N+1) if k not in (i,j)]
                if not ks: continue
                Ms = swap(M, i, ks[0]) if len({i,ks[0]}&set(M))==1 else swap(M, j, ks[0])
                if Ms == swap(M,i,j): continue
            w = list(base)
            w[j-1] = w[i-1] * (sp.Rational(5,4) if off else 1)
            for I in S:
                lhs = sp.cancel(E.Atilde(N, w, M, I, X))
                r0  = sp.cancel(E.Atilde(N, Ms and w or w, Ms, I, X))
                if sp.cancel(sp.together(lhs - rho*r0)) == 0:
                    ok += 1
                    if lhs == 0 and r0 == 0: triv += 1
                    elif r0 != 0: ratios.setdefault((i,j,M), []).append((I, sp.cancel(lhs/r0)))
                else:
                    bad += 1
    return ok, bad, triv, ratios

if __name__ == '__main__':
    X = sys.argv[1] if len(sys.argv)>1 else 'C'
    Ns = [int(x) for x in sys.argv[2:]] or [2,3]
    print('='*80); print(f'A-TYPE GKM CONGRUENCE, type {X}, my engine, q symbolic, w exact rational'); print('='*80)
    for N in Ns:
        for si, seed in enumerate(SEEDS):
            t0=time.time(); ok,bad,triv,rat = run(N, seed, X)
            print(f'N={N} seed{si}: holds {ok}/{ok+bad}  [{triv} trivial 0=0, {ok-triv} with a nonzero side]'
                  f'{"   *** %d FAILURES ***"%bad if bad else ""}   ({time.time()-t0:.0f}s)', flush=True)
        ok,bad,triv,rat = run(N, SEEDS[0], X)
        allone = True
        for (i,j,M), rs in sorted(rat.items()):
            vals=[v for _,v in rs]
            same = all(sp.cancel(v-vals[0])==0 for v in vals)
            isone = sp.cancel(vals[0]-1)==0
            allone &= bool(isone)
            print(f'   ratio (i,j)=({i},{j}) M={M}: {len(vals)} nonzero labels, I-independent={same}, ratio==1: {isone}'
                  + ('' if isone else f'   ratio = {sp.sstr(sp.factor(vals[0]))}'), flush=True)
        print(f'   ==> rho^A == 1 on every nonzero cell: {allone}', flush=True)
    print(); print('='*80); print('CONTROLS -- each must FIRE'); print('='*80)
    for N in Ns:
        ok,bad,_,_ = run(N, SEEDS[0], X, off=True)
        print(f'C1 N={N} OFF divisor (w_j = 5/4 w_i):  holds {ok} fails {bad}  -> {"FIRES" if bad else "DID NOT FIRE (BLIND)"}', flush=True)
    for N in Ns:
        ok,bad,triv,_ = run(N, SEEDS[0], X, rho=1+q)
        print(f'C2 N={N} PLANTED rho=1+q:  holds {ok} (of which {triv} are triangularity zeros, structurally spared)'
              f' fails {bad}  -> {"FIRES" if bad else "DID NOT FIRE (BLIND)"}', flush=True)
    for N in Ns:
        if N < 3: continue
        ok,bad,_,_ = run(N, SEEDS[0], X, wrong_partner=True)
        print(f'C3 N={N} WRONG PARTNER (s_ik M, k != i,j):  holds {ok} fails {bad}  -> {"FIRES" if bad else "DID NOT FIRE (BLIND)"}', flush=True)
