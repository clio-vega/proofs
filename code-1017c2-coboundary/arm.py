"""One chi-arm of the Theorem B census.  usage: python3 arm.py <label> <N> <X>"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1017c2-coboundary')
import sympy as sp
from collections import Counter
import q484_engine as E
q = E.q

CHI = {
  'id'  : lambda v: sp.Integer(1),                       # baseline
  'psi' : lambda v: (1 - q*v**2)/(1 - q),                # (COB) holds   -> predict 0
  'h1'  : lambda v: (1 - q*v)/(1 - q),                   # type B psi    -> plant, must fire
  'sq'  : lambda v: ((1 - q*v**2)/(1 - q))**2,           # weight -2     -> plant, must fire
  'unrm': lambda v: (1 - q*v**2),                        # scalar lesion -> must fire
  'xpsi': lambda v: v*(1 - q*v**2)/(1 - q),              # psi*x, anti-invariant twist -> predict 0
  'x2psi': lambda v: v**2*(1 - q*v**2)/(1 - q),  # psi * v^2: (COB) AND (SR) hold -> PREDICT CLEAN
  'x'    : lambda v: v,                              # (COB) fails, (SR) fails at -1 -> both pole types
}

def linext(S): return sorted(S, key=lambda M: (len(M), M))

def den_factors(f, w):
    num, den = sp.fraction(sp.cancel(sp.together(f)))
    out = Counter()
    for fac, mult in sp.factor_list(sp.expand(den))[1]:
        if len(sp.Poly(fac, *w).monoms()) == 1: continue   # monomial / constant
        out[sp.factor(fac)] += mult
    return out

def main(label, N, X):
    chi = CHI[label]; t0 = time.time()
    w = E.wvars(N); S = linext(E.subsets(N))
    A = {}
    for M in S:
        P = sp.Integer(1)
        for m in M: P *= chi(w[m-1])
        for I in S:
            A[(M,I)] = sp.cancel(sp.together(E.Atilde(N,w,M,I,X)*P))
    print(f'[{label}] N={N} {X}: basis built {time.time()-t0:.0f}s', flush=True)
    nonzero = 0; bad = {}
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
                df = den_factors(c[nu], w)
                if df: bad[(lam,mu,nu)] = df
            print(f'[{label}]   lam={lam} mu={mu} (nonzero={nonzero}, bad={len(bad)}) {time.time()-t0:.0f}s', flush=True)
    print(f'[{label}] RESULT nonzero={nonzero} nonLaurent={len(bad)}', flush=True)
    for t in sorted(bad):
        s = ' * '.join(f'({sp.sstr(f)})^{m}' if m>1 else f'({sp.sstr(f)})' for f,m in bad[t].items())
        print(f'[{label}] BAD lam={t[0]} mu={t[1]} nu={t[2]}  den: {s}', flush=True)
    print(f'[{label}] DONE {time.time()-t0:.0f}s', flush=True)

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
