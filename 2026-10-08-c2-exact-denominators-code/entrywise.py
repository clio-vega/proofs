"""
THE SURVIVING QUESTION.  Exact denominator of each entry of M^{-1}, not just each row.

(M^{-1})_{lam,nu} = (1/d_lam) * N_{lam,nu},  N = sum over the fibre F = {pi : lam^pi = nu}
of mu(hat0,pi).  S_lam = Aut(lam) acts on Pi_{l(lam)}, preserves F and preserves mu.
The action is NOT free (hat0 and hat1 are fixed by all of S_lam), so
    (M^{-1})_{lam,nu} = sum over S_lam-orbits O in F of  mu_O / |Stab(O)| .
CLAIM (orbit-stabiliser bound):
    den (M^{-1})_{lam,nu}  divides  lcm_O ( |Stab(O)| / gcd(|Stab(O)|, mu_O) ).
"""
from fractions import Fraction
from math import gcd, factorial
from lattice import *
from itertools import permutations
import sys
from functools import reduce

def lcm(a,b): return a*b//gcd(a,b) if a and b else 0

def aut_group(lam):
    """permutations w of [0,k) with lam[w(i)] == lam[i] for all i -- as tuples w"""
    k=len(lam)
    # group = product of symmetric groups on the blocks of equal parts
    blocks=[]
    i=0
    while i<k:
        j=i
        while j<k and lam[j]==lam[i]: j+=1
        blocks.append(list(range(i,j))); i=j
    out=[]
    def rec(idx, cur):
        if idx==len(blocks): out.append(dict(cur)); return
        B=blocks[idx]
        for perm in permutations(B):
            nxt=dict(cur)
            for a,b in zip(B,perm): nxt[a]=b
            rec(idx+1,nxt)
    rec(0,{})
    G=[tuple(w[i] for i in range(k)) for w in out]
    assert len(set(G))==len(G), ('NOT DISTINCT', lam, len(set(G)), len(G))
    for w in G: assert tuple(lam[w[i]] for i in range(k))==tuple(lam), ('NOT AN AUTOMORPHISM', lam, w)
    return G

def act(w, pi):
    return tuple(sorted((tuple(sorted(w[i] for i in B)) for B in pi), key=lambda B:(len(B),B)))

def canon(pi):
    return tuple(sorted((tuple(sorted(B)) for B in pi), key=lambda B:(len(B),B)))

NMAX = int(sys.argv[1]) if len(sys.argv)>1 else 7
rows=[]
bound_fail=[]; tight=0; loose=0; checked=0
eq_dlam=0; lt_dlam=0
loose_examples=[]
for n in range(1,NMAX+1):
    for lam in partitions(n):
        k=len(lam); dl=d_of(lam)
        G=aut_group(lam); assert len(G)==dl,(lam,len(G),dl)
        assert len(set(G))==dl, ('group has repeats', lam)
        SP=[canon(pi) for pi in set_partitions(k)]
        fibres={}
        for pi in SP: fibres.setdefault(coarsen(lam,pi),[]).append(pi)
        for nu,F in fibres.items():
            # exact entry
            N=sum(mu_hat0(pi) for pi in F)
            ent=Fraction(N,dl)
            den=ent.denominator
            # orbit decomposition
            seen=set(); bound=1
            for pi in F:
                if pi in seen: continue
                orb={act(w,pi) for w in G}
                seen|=orb
                stab=dl//len(orb)
                m=mu_hat0(pi)
                bound=lcm(bound, stab//gcd(stab,abs(m)) if m else 1)
            checked+=1
            if den!=0 and bound % den != 0:
                bound_fail.append((n,lam,nu,den,bound))
            if den==bound: tight+=1
            else:
                loose+=1
                if len(loose_examples)<8: loose_examples.append((n,lam,nu,'exact',den,'bound',bound,'entry',ent))
            if den==dl: eq_dlam+=1
            else: lt_dlam+=1

print("ENTRYWISE DENOMINATORS of M^{-1},  n <= %d" % NMAX)
print("  support entries examined          : %d" % checked)
print("  bound VIOLATED (den does not divide bound): %d" % len(bound_fail))
for b in bound_fail[:6]: print("    ", b)
print("  bound TIGHT (den == bound)        : %d" % tight)
print("  bound LOOSE (den  < bound)        : %d" % loose)
for e in loose_examples: print("    ", e)
print()
print("  entries attaining the ROW bound d_lam : %d" % eq_dlam)
print("  entries with a PROPER divisor of d_lam: %d   <-- so row-exactness is NOT entrywise" % lt_dlam)
