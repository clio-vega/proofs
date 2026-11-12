import sys; sys.path.insert(0,'.')
from crystal import *
def parts(k,maxlen):
    def rec(k,mx,acc):
        if k==0: yield tuple(acc); return
        if len(acc)==maxlen: return
        for p in range(min(k,mx),0,-1): yield from rec(k-p,p,acc+[p])
    yield from rec(k,k,[])
def u(shape):   # highest weight tableau rows
    return [[r+1]*l for r,l in enumerate(shape)]
def wt(w,n): return tuple(w.count(i) for i in range(1,n+1))
def LR(lam,mu,n):
    """hw elements b⊗u_mu of B(lam)⊗B(mu); returns dict nu -> list of b (as rows)"""
    um=rows_read(u(mu)); out={}
    for T in ssyt(lam,n):
        w=read(T,lam)+um
        if all(e(w,i) is None for i in range(1,n)):
            nu=wt(w,n)
            out.setdefault(nu,[]).append(tuple(tuple(T[(r,c)] for c in range(lam[r])) for r in range(len(lam))))
    return out
def split_rows(word,shape):
    """inverse of rows_read for given shape"""
    rows=[];p=0
    for l in reversed(shape): rows.append(tuple(word[p:p+l])); p+=l
    return tuple(reversed(rows))
def C1(b,lam,mu,n):
    """crystal-ops instrument"""
    xb=xi_crystal(rows_read(b),n); xc=xi_crystal(rows_read(u(mu)),n)
    x=xi_crystal(xc+xb,n)
    cprime=split_rows(x[:sum(mu)],mu); tail=split_rows(x[sum(mu):],lam)
    return cprime,tail
def inv_rsk(P,Q):
    P=[list(r) for r in P]; Q=[list(r) for r in Q]; m=sum(map(len,P)); w=[None]*m
    for t in range(m,0,-1):
        r=next(r for r,row in enumerate(Q) if row and row[-1]==t)
        Q[r].pop(); x=P[r].pop()
        if not P[r]: P.pop(); Q.pop()
        for rr in range(r-1,-1,-1):
            row=P[rr]; j=max(j for j,y in enumerate(row) if y<x)
            row[j],x=x,row[j]
        w[t-1]=x
    return tuple(w)
def C2(b,lam,mu,n):
    """evacuation + RSK instrument (no crystal operators)"""
    xb=rows_read(evac([list(r) for r in b],n)); xc=rows_read(evac(u(mu),n))
    P,Q=rsk(xc+xb)
    x=inv_rsk(evac(P,n),Q)
    return split_rows(x[:sum(mu)],mu), split_rows(x[sum(mu):],lam)
