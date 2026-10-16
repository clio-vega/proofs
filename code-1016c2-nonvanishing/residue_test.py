"""THEOREM D, tested as a closed-form prediction.

For nu MINIMAL in Motegi's order among support elements satisfying the bracket
    {a,b} subset nu,  nu' := nu\{a,b} >= lam,  nu' >= mu,
Theorem D predicts  ord_H c^nu = -1  with

  Res_H c^nu  =  rho (rho-1) * Atilde_{nu',lam}|_H * Atilde_{nu',mu}|_H
                 / [ (w_a w_b - 1)^{-1} Atilde_{nu,nu} ]|_H .

Nothing here is fitted: rho comes from eq:branch, (rho-1) from one line of
algebra, and the three Atilde's are evaluated by the brute-force S_n engine.

PLANTS (pre-registered, in DIFFERENT blocks of the formula):
  P1  multiplier block :  rho(rho-1) -> rho(rho-1)(1+q)      must FAIL
  P2  multiplier block :  rho(rho-1) -> rho^2                must FAIL  [this is
      the "untwisted" reading -- the one a collision divisor would give]
  P3  label block      :  Atilde_{nu',mu} -> Atilde_{nu',lam} must FAIL on the
      triples with lam != mu, and is INFORMATIVELY SPARED when lam == mu.
Each plant's domain is named with it.
"""
import sys; sys.path.insert(0,'.')
import sympy as sp, time
from collections import Counter
import q484_engine as E
from gate_q0 import build, solve_c
q = E.q
s = sp.Symbol('s')

def toH(expr, w, a, b, k=0):
    """[ (w_a w_b - 1)^{-k} expr ] restricted to H_ab ; None if that is singular."""
    expr = sp.cancel(sp.together(expr))
    if expr == 0: return sp.Integer(0)
    e = sp.cancel(expr.subs(w[b-1], s/w[a-1]))
    e = sp.cancel(sp.together(e*(s-1)**(-k)))
    n,d = sp.fraction(e)
    n1 = sp.expand(n).subs(s,1); d1 = sp.expand(d).subs(s,1)
    if d1 == 0: return None
    return sp.cancel(n1/d1)

def rho_of(w,a,b):
    wa=w[a-1]; return (1-q)**2/((1-q*wa**2)*(1-q/wa**2))

def run(N, X='C'):
    t0=time.time(); w,S,A = build(N,X)
    print(f'=== THEOREM D  N={N} {X}  (build {time.time()-t0:.0f}s) ===',flush=True)
    res = Counter(); tested=[]
    for lam in S:
      for mu in S:
        if mu < lam: continue
        c = solve_c(S,A,lam,mu)
        supp = [k for k in S if c[k]!=0]
        for a in range(1,N+1):
         for b in range(a+1,N+1):
          def bracket(nu):
              if a not in nu or b not in nu: return False
              nup = tuple(x for x in nu if x not in (a,b))
              return E.geq(nup,lam) and E.geq(nup,mu)
          B = [nu for nu in supp if bracket(nu)]
          for nu in B:
            # minimal: no kappa in supp, kappa <= nu, kappa != nu, bracket(kappa)
            if any(k!=nu and E.geq(nu,k) and bracket(k) for k in B): continue
            nup = tuple(x for x in nu if x not in (a,b))
            r0 = rho_of(w,a,b)
            dII = toH(A[(nu,nu)], w,a,b, k=1)          # [t^{-1} A_nunu]|_H
            Al  = toH(A[(nup,lam)], w,a,b); Am = toH(A[(nup,mu)], w,a,b)
            if dII is None or dII==0 or Al is None or Am is None:
                res['SKIP(degenerate)']+=1; continue
            pred  = sp.cancel(sp.together(r0*(r0-1)*Al*Am/dII))
            predH = toH(pred, w,a,b)
            actual = toH(c[nu], w,a,b, k=-1)           # Res_H c^nu = lim (t) c
            ok = sp.simplify(sp.cancel(actual-predH))==0
            res['MATCH' if ok else 'MISMATCH']+=1
            # plants
            p1 = toH(sp.cancel(r0*(r0-1)*(1+q)*Al*Am/dII),w,a,b)
            p2 = toH(sp.cancel(r0*r0*Al*Am/dII),w,a,b)
            p3 = toH(sp.cancel(r0*(r0-1)*Al*Al/dII),w,a,b)
            res['P1-fires' if sp.simplify(actual-p1)!=0 else 'P1-MISSED']+=1
            res['P2-fires' if sp.simplify(actual-p2)!=0 else 'P2-MISSED']+=1
            if lam!=mu:
                res['P3-fires' if sp.simplify(actual-p3)!=0 else 'P3-MISSED']+=1
            else:
                res['P3-spared(lam=mu, by construction)']+=1
            # and the order
            res[f'ord_H c^nu = {-1 if actual!=0 else "0 (residue vanished!)"}']+=1
            tested.append((lam,mu,nu,(a,b),ok))
            print(f'  lam={lam} mu={mu} nu={nu} H_{a}{b}: '
                  f'{"MATCH" if ok else "MISMATCH"}   Res={sp.factor(actual)}',flush=True)
    print(f'\n--- summary N={N} {X} ---')
    for k,v in sorted(res.items()): print(f'  {k:42s} {v}')
    print(f'  minimal-nu cases tested: {len(tested)}   ({time.time()-t0:.0f}s)')

if __name__=='__main__':
    run(int(sys.argv[1]), sys.argv[2] if len(sys.argv)>2 else 'C')
