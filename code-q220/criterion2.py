"""Corrected height model.

Slot z is read as ')' (if z in T) then '(' (if z in S).
  P(z) = height after slot z ;  Q(z) = P(z-1) - [z in T]  = the LOW point of slot z.
  Q(z+n) = Q(z) + d,  d = |S| - |T|.

LEMMA.  Fix a legal cut x (so x not in S u T), window = slots x+1..x+n, P(x)=0.
  Let mu = min{ Q(z) : x < z <= x+n } and z* = the LAST z attaining mu.
  Then L_1^(x) is nonempty iff z* in S, and in that case b_x = z*.

CHARACTERISATION.  c in B  iff  there is a legal cut x with c-n <= x < c,
  G(c) <= x <= F(c)-n-1, where
  F(c) = min{ z>c : Q(z) <= Q(c) },   G(c) = max{ z<c : Q(z) < Q(c) }.
"""
import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *

INF=float('inf')

def Qfun(S,T,n):
    Q={}; P=0
    for z in range(n):
        q = P - (1 if z in T else 0)
        Q[z]=q
        P = q + (1 if z in S else 0)
    d=len(S)-len(T)
    return (lambda z: Q[z%n] + (z//n)*d), d

def b_from_Q(S,T,x,n,Q):
    win=list(range(x+1,x+n+1))
    mu=min(Q(z) for z in win)
    zs=max(z for z in win if Q(z)==mu)
    return zs%n if zs%n in S else None

# LEMMA check
bad=0; ck=0
for n in range(3,8):
    for S,T,_ in additive_pairs(n):
        Q,d=Qfun(S,T,n)
        for x in X_set(S,T,n):
            r=ms_etilde(S,T,x,n); ck+=1
            got=b_from_Q(S,T,x,n,Q)
            exp=None if r is None else r[0]
            if got!=exp: bad+=1
print(f"LEMMA (b_x = last argmin of Q, if in S): {ck} cuts checked, {bad} failures")

def FG(S,T,n,Q,c):
    R=5*n+5
    F=INF
    for z in range(c+1,c+R):
        if Q(z)<=Q(c): F=z; break
    G=-INF
    for z in range(c-1,c-R,-1):
        if Q(z)<Q(c): G=z; break
    return F,G

def B_pred(S,T,n):
    Q,d=Qfun(S,T,n); X=set(X_set(S,T,n)); out=set()
    for c in clio_letters(S,T,n):
        F,G=FG(S,T,n,Q,c)
        lo=max(G,c-n); hi=min(F-n-1,c-1)
        if lo<=hi and any((x%n) in X for x in range(int(lo),int(hi)+1)): out.add(c)
    return out

def B_act(S,T,n):
    out=set()
    for x in X_set(S,T,n):
        r=ms_etilde(S,T,x,n)
        if r is not None: out.add(r[0])
    return out

tot=collections.Counter(); bad2=[]
for n in range(3,8):
    N=ok=0
    for S,T,_ in additive_pairs(n):
        if not X_set(S,T,n): continue
        N+=1
        p,a=B_pred(S,T,n),B_act(S,T,n)
        if p==a: ok+=1
        elif len(bad2)<5: bad2.append((n,S,T,p,a))
    print(f"n={n}: criterion predicts B exactly on {ok}/{N}")
    tot['N']+=N; tot['ok']+=ok
print("TOTAL",dict(tot))
for b in bad2: print("  MISS n=%d S=%s T=%s pred=%s act=%s"%(b[0],sorted(b[1]),sorted(b[2]),sorted(b[3]),sorted(b[4])))
