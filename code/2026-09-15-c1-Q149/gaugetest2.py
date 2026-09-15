"""
Cocycle test, done properly: ONE undirected connected component, single base point.
eta(M->M') = W(u)/W_MN(u) on every e-move and f-move.  eta is a coboundary iff
d_M can be assigned consistently.  Any inconsistency is a genuine obstruction.
"""
import sympy as sp
from collections import deque
from engine import all_words

def moves_up(S, e):
    out=[]
    for b in sorted(S):
        if (b+e) in S: continue
        u = tuple(1 if (b+i) in S else 0 for i in range(1,e))
        out.append(((S-{b})|frozenset([b+e]), u))
    return out

def cocycle_test(e, f, W, Wb, CUT=None, TOP=None, nmax=40000):
    g=e+f
    if CUT is None: CUT=-2*g-4
    if TOP is None: TOP=3*g          # states confined to [CUT, TOP)
    MN  = {u: sp.Integer(-1)**sum(u) for u in all_words(e-1)}
    MNb = {u: sp.Integer(-1)**sum(u) for u in all_words(f-1)}
    ETA = {e:{u: sp.cancel(W[u]/MN[u]) for u in W}, f:{u: sp.cancel(Wb[u]/MNb[u]) for u in Wb}}
    start = frozenset(range(CUT,0))
    d={start: sp.Integer(1)}; Q=deque([start]); viol=[]; nedge=0
    while Q and len(d)<nmax:
        S=Q.popleft()
        for ee in (e,f):
            for (S2,u) in moves_up(S, ee):
                if max(S2) >= TOP: continue
                nedge+=1
                want = sp.cancel(d[S]*ETA[ee][u])
                if S2 in d:
                    r = sp.cancel(d[S2]/want)
                    if sp.simplify(r-1)!=0: viol.append((r, ee, u))
                else:
                    d[S2]=want; Q.append(S2)
    return len(d), nedge, viol
