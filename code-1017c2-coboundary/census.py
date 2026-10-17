"""Q505 / Theorem B: does the coboundary renormalisation remove the long-root
poles of the STRUCTURE CONSTANTS?

Basis:        Atilde_{M,I}                (Motegi type C, E^C normaliser)
Renormalised: Ahat_{M,I} = Atilde_{M,I} * Psi(M),   Psi(M) = prod_{m in M} chi(w_m)

Theory (this session).  On H_ab = {w_a w_b = 1}, QC gives Atilde_{M,I} = g * Atilde_{M'',I}
for every I, g = (1-q)^2/((1-q w_a^2)(1-q w_a^-2)).  A product of two basis elements
therefore collides with multiplier g^2, while the span collides with g: the mismatch
g^2 - g != 0 is exactly the pole (Rick, Prop 5 / Thm 6).  If chi satisfies
    (COB)   chi(x) chi(1/x) = 1/g   identically,
then Psi has collision weight -1, Psi * (phi_lam phi_mu) has weight -1+2 = +1 = the weight
of the span, and every renormalised structure constant is REGULAR on H_ab.

PREDICTION (registered before running): with chi = psi, chi(w) = (1-q w^2)/(1-q),
all 13 non-Laurent constants at N=3 become Laurent; non-Laurent set = EMPTY.

CONTROLS.  Four more chi, three of which must NOT work:
  id    chi = 1                       baseline, must reproduce 13 non-Laurent
  psi   chi = (1-q w^2)/(1-q)         (COB) holds       -> predicted 0
  h1    chi = (1-q w)/(1-q)           type B's psi; (COB) FAILS   -> plant, must fire
  sq    chi = ((1-q w^2)/(1-q))^2     weight -2         -> plant, must fire
  unrm  chi = (1-q w^2)               weight -1 up to the scalar (1-q)^2 -> lesion, must fire
Reporting is by SET (the triples), never by count, and the symmetric difference
against the baseline is printed.
"""
import sys, time, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1017c2-coboundary')
import sympy as sp
from collections import Counter
import q484_engine as E
q = E.q

def linext(S): return sorted(S, key=lambda M: (len(M), M))

def chi_family(w):
    x = sp.Symbol('x')
    return {
      'id'  : lambda v: sp.Integer(1),
      'psi' : lambda v: (1 - q*v**2)/(1 - q),
      'h1'  : lambda v: (1 - q*v)/(1 - q),
      'sq'  : lambda v: ((1 - q*v**2)/(1 - q))**2,
      'unrm': lambda v: (1 - q*v**2),
    }

def Psi(chi, w, M):
    P = sp.Integer(1)
    for m in M: P *= chi(w[m-1])
    return P

def solve_c(S, A, lam, mu):
    c = {}
    for nu in S:
        v = A[(nu,lam)]*A[(nu,mu)]
        acc = sp.Integer(0)
        for kap in S:
            if kap == nu: break
            if c.get(kap, 0) != 0: acc += c[kap]*A[(nu,kap)]
        c[nu] = sp.cancel(sp.together((v-acc)/A[(nu,nu)]))
    return c

def den_factors(f, w):
    f = sp.cancel(sp.together(f))
    num, den = sp.fraction(f)
    out = Counter()
    for fac, mult in sp.factor_list(sp.expand(den))[1]:
        p = sp.Poly(fac, *w)
        if len(p.monoms()) == 1:   # monomial in w (incl. constants / pure q)
            continue
        out[sp.factor(fac)] += mult
    return out

def census(N, X, label, chi, w, S, A0, verbose=True):
    A = {k: sp.cancel(sp.together(v*Psi(chi, w, k[0]))) for k, v in A0.items()}
    nonzero, bad = 0, {}
    for lam in S:
        for mu in S:
            if mu < lam: continue
            c = solve_c(S, A, lam, mu)
            for nu in S:
                if c[nu] == 0: continue
                nonzero += 1
                df = den_factors(c[nu], w)
                if df: bad[(lam,mu,nu)] = df
    print(f'[{label}] nonzero={nonzero}  non-Laurent={len(bad)}', flush=True)
    if verbose:
        for t in sorted(bad):
            s = ' * '.join(f'({sp.sstr(f)})^{m}' if m>1 else f'({sp.sstr(f)})'
                           for f,m in bad[t].items())
            print(f'    lam={t[0]} mu={t[1]} nu={t[2]}  den: {s}', flush=True)
    return set(bad), bad

def main(N=3, X='C'):
    t0 = time.time()
    w = E.wvars(N); S = linext(E.subsets(N))
    A0 = {(M,I): E.Atilde(N,w,M,I,X) for M in S for I in S}
    print(f'=== N={N} type {X}  (build {time.time()-t0:.0f}s) ===', flush=True)
    fam = chi_family(w)
    sets = {}
    for label in ['id','psi','h1','sq','unrm']:
        t1 = time.time()
        sets[label], _ = census(N, X, label, fam[label], w, S, A0)
        print(f'    ({time.time()-t1:.0f}s)', flush=True)
    base = sets['id']
    print('\n=== SET comparison against baseline id ===', flush=True)
    for label in ['psi','h1','sq','unrm']:
        s = sets[label]
        print(f'  {label:5s}  |set|={len(s):3d}  removed={len(base-s):3d}  added={len(s-base):3d}')
        if base - s: print(f'         removed: {sorted(base-s)}')
        if s - base: print(f'         added  : {sorted(s-base)}')
    return sets

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 3,
         sys.argv[2] if len(sys.argv)>2 else 'C')
