"""Verification suite, 2026-10-17.

A: E^B implementation checked against Motegi's INDEPENDENT theorem typeBmainproperties(i)
   (normalizationWB L1386) and (ii) (typeBfullDivisibility L1394) -- a different theorem
   from the normalisation factor, so this is corroboration not a second reading.
   PLANT A1: EB with the w-diagonal inclusive (j<=k), i.e. the type-C convention -> must FAIL.
   PLANT A2: EB without the extra prod_j (1-q u_j)                        -> must FAIL.
B: the type-B collapse multiplier at N=4, and reproduction of the paper's wrong 48/48.
C: q=0 structure claim: Atilde|_{q=0} depends on w_k (k in Mc) only if k < max(I).
   PLANT C1: claim with "k <= max(I)"  -> must FAIL somewhere.
D: coverage of the stratum max(I) <= b, at N=2..5 (combinatorial, three N's minimum).
"""
import sympy as sp, itertools
import q484_engine as E
from q484_engine import q, wvars, WX, EC
from eb_audit import EB, Atilde_B_true

def EB_plant(N,n,r,u,w,which):
    num = sp.Integer(1)
    for j in range(n):
        for k in range(n): num *= (1 - q*u[k]/u[j])
    for j in range(n):
        for k in range(N): num *= (1 - q*u[j]*w[k])
    if which != 'A2':
        for j in range(n): num *= (1 - q*u[j])
    den = sp.Integer(1)
    for j in range(n):
        for k in range(j+1,n): den *= (1 - q*u[j]*u[k])
    for j in range(n): den *= (1 - q*u[j]**2)
    lo = 0 if which=='A1' else 1          # A1: inclusive w-diagonal (type C convention)
    for j in range(N):
        for k in range(j+lo, N): den *= (1 - q*w[j]*w[k])
    for k in range(N): den *= (1 - q*w[k])
    return num/den/(1-q)**(n-r)

def AtB(N,w,M,I,EBfun):
    s,r = len(M),len(I)
    if s<r: return sp.Integer(0)
    if s==0: return sp.cancel(1/EBfun(N,0,0,[],w))
    u=[w[m-1] for m in M]
    return sp.cancel(sp.together(WX(N,s,r,u,w,I,'B')/EBfun(N,s,r,u,w)))

def diag_rootproduct_B(N,w,I):
    """prod_{Inv}(1-e^-a) prod_{T}(1-q e^+a) in the TYPE B reading, read off
    typeBfullDivisibility L1394-1402 at M=I (short roots t_k replace long roots 2t_k)."""
    r=len(I); Il=sorted(I); Ic=[k for k in range(1,N+1) if k not in I]
    t=sp.Integer(1)
    for k in Il: t *= (1 - 1/w[k-1])                      # short-root inversion  (1-w^-1)
    for j in range(r):
        for k in range(j+1,r): t *= (1 - 1/(w[Il[j]-1]*w[Il[k]-1]))
    for j in Ic: t *= (1 - q*w[j-1])                      # short-root ascent     (1-q w)
    for j in Ic:
        for k in Ic:
            if j<k: t *= (1 - q*w[j-1]*w[k-1])
    for j in Ic:
        for l in Il:
            if l>j: t *= (1 - 1/w[l-1]*w[j-1])
    for j in Ic:
        for l in Il:
            if l<j: t *= (1 - q*w[j-1]/w[l-1])
    return sp.cancel(t)

print("=== A: E^B vs Motegi's independent type-B diagonal theorem (normalizationWB) ===")
for N in (2,3):
    w=wvars(N)
    for tag,fn in (('EB (mine)',EB),('PLANT A1 (incl w-diag)',lambda *a: EB_plant(*a,'A1')),
                                    ('PLANT A2 (no 1-qu_j)',lambda *a: EB_plant(*a,'A2'))):
        ok=bad=0
        for I in E.subsets(N):
            got = AtB(N,w,I,I,fn); pred = diag_rootproduct_B(N,w,I)
            # compare up to the overall sign/monomial convention: test exact equality
            if sp.simplify(got-pred)==0: ok+=1
            else: bad+=1
        print(f"  N={N} {tag}: match {ok}, mismatch {bad}")
