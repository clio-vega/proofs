"""Independent symbolic check of Theorem A's PROOF (not just its conclusion).

The engine test (collapse_test.py) verifies the CONCLUSION -- that rho is the
right multiplier -- without using eq:branch at all.  This script checks the
PROOF: it implements Psi_{K,L,Mc} literally from the table in the paper and
verifies

    Psi_{K, L' u {a,b}, Mc} / Psi_{K, L', Mc u {a,b}}   restricted to H_ab
      ==  (1-q)^2 / ((1-q w_a^2)(1-q w_a^-2))

for every (N, a<b, M containing {a,b}, K subset M\{a,b}).  A typo in my
factor-by-factor table shows up here even though the conclusion is already
corroborated.  Note Psi carries NO I-dependence (Lemma "slots"), so there is no
I loop -- that absence is itself the lemma.

CONTROL: the same quotient is also computed OFF H (at a random-free generic
specialisation of w_b) and must NOT equal rho -- otherwise the restriction is
doing no work and the test would pass for a trivial reason.
"""
import sys; sys.path.insert(0,'.')
import sympy as sp
from itertools import combinations
q = sp.symbols('q')

def Psi(w, N, K, L, Mc):
    """(1-q)^{|L|} * Pi1 * Pi3 * Pi4 * Pi5, exactly as tabulated."""
    P = (1-q)**len(L)
    for j in Mc:                                    # Pi1
        for k in Mc:
            if j <= k: P *= (1 - q*w[j-1]*w[k-1])
    for ka in K:                                    # Pi3
        for l in L:
            P *= (1 - 1/(w[ka-1]*w[l-1])) / (1 - w[l-1]/w[ka-1])
    for la, lb in combinations(sorted(L), 2):       # Pi4
        P *= (1 - q/(w[la-1]*w[lb-1]))
    for l in L:                                     # Pi5
        for j in Mc:
            P *= (1 - q*w[j-1]/w[l-1])
    return P

def run(N):
    w = list(sp.symbols(f'w1:{N+1}'))
    s = sp.Symbol('s')
    ok=0; bad=0; ctrl_ok=0; ctrl_bad=0; badlist=[]
    for a in range(1,N+1):
      for b in range(a+1,N+1):
        rho = (1-q)**2/((1-q*w[a-1]**2)*(1-q/w[a-1]**2))
        for sz in range(2,N+1):
          for M in combinations(range(1,N+1), sz):
            if a not in M or b not in M: continue
            Mp = tuple(x for x in M if x not in (a,b))
            Mc = tuple(x for x in range(1,N+1) if x not in M)
            Mpc = tuple(sorted(Mc + (a,b)))
            for r in range(0, len(Mp)+1):
              for K in combinations(Mp, r):
                Lp = tuple(x for x in Mp if x not in K)
                L  = tuple(sorted(Lp + (a,b)))
                ratio = sp.cancel(Psi(w,N,K,L,Mc) / Psi(w,N,K,Lp,Mpc))
                # restrict to H:  w_b = s/w_a, then s -> 1
                e = sp.cancel(ratio.subs(w[b-1], s/w[a-1]))
                onH = sp.simplify(sp.cancel(e.subs(s,1) - rho))
                if onH == 0: ok+=1
                else:
                    bad+=1
                    if len(badlist)<5: badlist.append((a,b,M,K))
                # CONTROL: off H, at s = 2 (a generic point of the w_b line)
                offH = sp.simplify(sp.cancel(e.subs(s,2) - rho))
                if offH != 0: ctrl_ok+=1
                else: ctrl_bad+=1
    print(f'  N={N}:  ratio|_H == rho in {ok}/{ok+bad} cases   ({bad} FAIL)')
    for z in badlist: print('     fail a,b,M,K =',z)
    print(f'  N={N}:  CONTROL off H (s=2): ratio != rho in {ctrl_ok}/{ctrl_ok+ctrl_bad}'
          f'  -- {ctrl_bad} cases where the restriction did NO work')

if __name__=='__main__':
    for N in range(2, int(sys.argv[1])+1): run(N)
