"""SUPERSEDED -- kept as the record of a refuted first attempt.

This version tested the criterion using only the S-slots (alpha), and MISSES
175 of 9093 pairs: a letter can be the last alpha-minimiser among S-slots and
still be matched, because the path can return to its level at a ')' slot that
carries no '('.  The corrected version, which compares the full slot-low-point
profile Q, is in criterion2.py and agrees on 9093/9093.
"""
import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *
from heights import alpha_fun

INF=float('inf')

def fg(S,T,n,A,c,d):
    # search a generous range of lifts
    R = 4*n+4
    f=INF
    for j in range(c+1, c+R):
        if j%n in S and A(j)<=A(c): f=j; break
    g=-INF
    for j in range(c-1, c-R, -1):
        if j%n in S and A(j)<A(c): g=j; break
    return f,g

def B_predicted(S,T,n):
    A,d=alpha_fun(S,T,n)
    X=set(X_set(S,T,n))
    out=set()
    for c0 in clio_letters(S,T,n):
        c=c0
        f,g=fg(S,T,n,A,c,d)
        lo=max(g, c-n+1); hi=min(f-n, c-1)
        ok=any((x%n) in X for x in range(int(lo), int(hi)+1)) if lo<=hi else False
        if ok: out.add(c0)
    return out

def B_actual(S,T,n):
    out=set()
    for x in X_set(S,T,n):
        r=ms_etilde(S,T,x,n)
        if r is not None: out.add(r[0])
    return out

tot=collections.Counter(); bad=[]
for n in range(3,8):
    N=0; ok=0
    for S,T,_ in additive_pairs(n):
        if not X_set(S,T,n): continue
        N+=1
        p,a=B_predicted(S,T,n), B_actual(S,T,n)
        if p==a: ok+=1
        elif len(bad)<5: bad.append((n,S,T,p,a))
    print(f"n={n}: {ok}/{N} pairs where the criterion predicts B exactly")
    tot['N']+=N; tot['ok']+=ok
print("TOTAL",dict(tot))
for n,S,T,p,a in bad: print("  MISS n=%d S=%s T=%s pred=%s act=%s"%(n,sorted(S),sorted(T),sorted(p),sorted(a)))
