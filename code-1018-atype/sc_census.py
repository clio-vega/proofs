"""STEP 3.  Structure-constant census, denominators CLASSIFIED BY DIVISOR.

Same triangular solve as ../code-1017c2-coboundary/arm.py (so the 13/115 and 0/115
baselines are reproducible), but every non-monomial denominator factor is classified
A / pair / bdy / qdef by divisors.classify -- which decides by VANISHING, not by
pattern-matching the printed form.

Reported as a SET of triples per divisor class, never as a count:
two criteria can fit the same count and differ as sets.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1018-atype')
import sympy as sp
import q484_engine as E, divisors as D
q = E.q

CHI = {
  'id'   : lambda v: sp.Integer(1),
  'psi'  : lambda v: (1 - q*v**2)/(1 - q),
  'h1'   : lambda v: (1 - q*v)/(1 - q),
  'sq'   : lambda v: ((1 - q*v**2)/(1 - q))**2,
  'unrm' : lambda v: (1 - q*v**2),
  'xpsi' : lambda v: v*(1 - q*v**2)/(1 - q),
  'x2psi': lambda v: v**2*(1 - q*v**2)/(1 - q),
  'x'    : lambda v: v,
}

def main(label, N, X='C'):
    t0 = time.time()
    w = E.wvars(N); S = sorted(E.subsets(N), key=lambda M:(len(M),M))
    PLANT = (label == 'atp')          # P4 control: violates the unit hypothesis
    chi = CHI['id'] if PLANT else CHI[label]
    A = {}
    for M in S:
        P = sp.Integer(1)
        for m in M: P *= chi(w[m-1])
        if PLANT and len(M) % 2 == 1: P = P/(w[0]-w[1])
        for I in S:
            A[(M,I)] = sp.cancel(sp.together(E.Atilde(N,w,M,I,X)*P))
    print(f'[{label}] N={N} {X} basis built ({time.time()-t0:.0f}s)'
          + ('   *** A-TYPE PLANT ACTIVE (P4 control) ***' if PLANT else ''), flush=True)
    nonzero = 0
    byclass = {}          # kind -> set of (lam,mu,nu)
    detail  = {}
    for lam in S:
        for mu in S:
            if mu < lam: continue
            c = {}
            for nu in S:
                v = A[(nu,lam)]*A[(nu,mu)]
                acc = sp.Integer(0)
                for kap in S:
                    if kap == nu: break
                    if c.get(kap,0) != 0: acc += c[kap]*A[(nu,kap)]
                c[nu] = sp.cancel(sp.together((v-acc)/A[(nu,nu)]))
                if c[nu] == 0: continue
                nonzero += 1
                kinds = D.den_classify(c[nu], w, q)
                for k, lst in kinds.items():
                    byclass.setdefault(k, set()).add((lam,mu,nu))
                    detail.setdefault(k, []).append(((lam,mu,nu), lst))
    allbad = set().union(*byclass.values()) if byclass else set()
    print(f'[{label}] RESULT nonzero={nonzero}  non-Laurent={len(allbad)}', flush=True)
    for k in ('A','pair','bdy','qdef','other'):
        s = byclass.get(k, set())
        print(f'[{label}]   {k:5s}: {len(s):3d} triples', flush=True)
    print(f'[{label}] A-TYPE SET = {sorted(byclass.get("A", set()))}', flush=True)
    for t, lst in detail.get('A', []):
        print(f'[{label}]   A-POLE lam={t[0]} mu={t[1]} nu={t[2]}: {lst}', flush=True)
    print(f'[{label}] PAIR SET = {sorted(byclass.get("pair", set()))}', flush=True)
    print(f'[{label}] DONE {time.time()-t0:.0f}s', flush=True)

if __name__=='__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv)>2 else 3,
         sys.argv[3] if len(sys.argv)>3 else 'C')
