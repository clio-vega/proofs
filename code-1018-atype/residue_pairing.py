"""STEP 4.  Lemma A2 (regularity on A-type divisors), mechanism check.

F_{N,n,r}(u|w|I) = (1/(n-r)!) * prod_k u_k^{-(N+n)} * SUM_{sigma in S_n} T(u_sigma).
Every dangerous factor is (q v_a - v_b)/(v_a - v_b) for a SLOT pair {a,b} with
min(a,b) <= r.  Fix two u-values u_p, u_p'.  For each sigma let a,b be the slots
holding p,p' and set sigma' = sigma . (ab).  Claim: sigma <-> sigma' is a fixed-point-
free involution and Res_{u_p = u_p'} [ T(u_sigma) + T(u_sigma') ] = 0, because the two
dangerous factors have residues +x(q-1) and -x(q-1) while EVERY other factor becomes
literally the same expression at u_p = u_p' (slots a,b carry equal values there).

Tested two ways:
  (T1) the FULL sum has no pole at u_p = u_p'   (the conclusion)
  (T2) each involution PAIR {sigma, sigma'} is already pole-free  (the mechanism)
CONTROL (must FIRE): a partial sum over a set of sigma that is NOT closed under the
involution must HAVE a pole.  Sensitive to: whether the cancellation is pairwise at
all, as opposed to a global accident.  If this does not fire, T2 carries no
information about the mechanism.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1018-atype')
import sympy as sp
from itertools import permutations
from math import factorial
import q484_engine as E
q = E.q

def Tterm(N, n, r, u, w, I, sigma, beta=sp.Integer(1), gamma=sp.Integer(0)):
    """ONE summand T(u_sigma) of Motegi's Motherfunction -- transcribed from
    q484_engine.F so that the two agree by construction when summed."""
    us = [u[sigma[j]] for j in range(n)]
    t = sp.Integer(1)
    for j in range(r):
        ij = I[j]
        t *= (1-q)*(us[j]**2 - 1)
        for k in range(ij+1, N+1): t *= (us[j] - q*w[k-1])
        for k in range(1, ij):     t *= (us[j] - w[k-1])
    for j in range(r):
        for k in range(j+1, r):
            t *= (q*us[j]-us[k])/(us[j]-us[k])*(us[j]*us[k]-1)
    for j in range(r, n):
        for k in range(N): t *= (us[j] - q*w[k])
    for j in range(r, n):
        for k in range(r):
            t *= (q*us[j]-us[k])/(us[j]-us[k])*(us[j]*us[k]-1)
    for j in range(r, n):
        for k in range(j+1, n): t *= (us[j]*us[k]-q)
    for j in range(r, n): t *= (beta*us[j]-gamma)
    return t

def ordA(f, w, i, j):
    f = sp.cancel(sp.together(f))
    if f == 0: return sp.oo
    num, den = sp.fraction(f)
    def mult(p):
        k=0; p=sp.expand(p)
        while True:
            qq, rr = sp.div(sp.Poly(p, w[j-1]), sp.Poly(w[j-1]-w[i-1], w[j-1]))
            if sp.expand(rr.as_expr())!=0: return k
            p=qq.as_expr(); k+=1
    return mult(num)-mult(den)

def involution(sigma, p, pp):
    """slots of p and pp, swapped."""
    a = sigma.index(p); b = sigma.index(pp)
    s = list(sigma); s[a], s[b] = s[b], s[a]
    return tuple(s)

def main(N=3, M=(1,2,3), I=(1,), X='C'):
    w = E.wvars(N); n = len(M); r = len(I)
    u = [w[m-1] for m in M]
    i, j = M[0], M[1]              # the two colliding alphabet indices, both IN M
    p, pp = 0, 1                   # their slots in u
    beta, gamma = (sp.Integer(1), sp.Integer(0))
    perms = list(permutations(range(n)))
    print(f'=== N={N} M={M} I={I} type {X}: collision w_{i} = w_{j}, both in M ===', flush=True)
    # individual summands
    ords = [ordA(Tterm(N,n,r,u,w,I,s,beta,gamma), w, i, j) for s in perms]
    print(f'  individual summands: {len(perms)} terms, ord_(w{i}-w{j}) = {sorted(set(ords))}'
          f'  -> {sum(1 for o in ords if o<0)} of {len(perms)} have a POLE', flush=True)
    # T1 full sum
    tot = sum(Tterm(N,n,r,u,w,I,s,beta,gamma) for s in perms)
    o = ordA(tot, w, i, j)
    print(f'  T1  FULL sum over S_{n}: ord = {o}  -> {"REGULAR" if o>=0 else "*** POLE ***"}', flush=True)
    # T2 involution pairs
    seen=set(); bad=0; npair=0
    for s in perms:
        if s in seen: continue
        s2 = involution(s, p, pp); seen.add(s); seen.add(s2)
        assert s2 != s, 'involution has a fixed point -- claim is wrong'
        npair += 1
        o2 = ordA(Tterm(N,n,r,u,w,I,s,beta,gamma)+Tterm(N,n,r,u,w,I,s2,beta,gamma), w, i, j)
        if o2 < 0: bad += 1
    print(f'  T2  involution PAIRS: {npair} pairs (fixed-point free: True), '
          f'{npair-bad} regular, {bad} still polar  -> {"MECHANISM CONFIRMED" if bad==0 else "*** MECHANISM FAILS ***"}', flush=True)
    # CONTROL: a non-closed subset must have a pole
    half = [s for s in perms if s.index(p) < s.index(pp)]   # NOT involution-closed
    oc = ordA(sum(Tterm(N,n,r,u,w,I,s,beta,gamma) for s in half), w, i, j)
    print(f'  C   partial sum over {len(half)} sigma with slot(p)<slot(p\'), NOT closed under'
          f' the involution: ord = {oc}  -> {"FIRES (has a pole)" if oc<0 else "DID NOT FIRE (BLIND)"}', flush=True)

if __name__=='__main__':
    for args in [(3,(1,2,3),()), (3,(1,2,3),(1,)), (3,(1,2,3),(2,3)), (3,(1,2),(1,)), (3,(1,2),(1,2))]:
        main(args[0], args[1], args[2])
        print()
