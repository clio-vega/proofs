"""Crystal commutor on LR sets, two independent instruments.
Convention: element of B(a1)⊗...⊗B(ak) <-> concatenated word read(a1)...read(ak);
read(T) = rows bottom-to-top, each left-to-right.  Bracket rule: letter i+1 = '(' , letter i = ')'.
f_i: rightmost unpaired i -> i+1.  e_i: leftmost unpaired i+1 -> i.
Henriques-Kamnitzer commutor: sigma(b⊗c) = xi_{B(mu)⊗B(lam)}( xi(c) ⊗ xi(b) ).
"""
from itertools import product
def unpaired(w,i):
    st=[];up_i=[]           # stack of unmatched '(' (i+1) positions
    for p,x in enumerate(w):
        if x==i+1: st.append(p)
        elif x==i:
            if st: st.pop()
            else: up_i.append(p)
    return up_i, st          # unpaired i's (left part), unpaired i+1's
def f(w,i):
    a,b=unpaired(w,i)
    if not a: return None
    w=list(w); w[a[-1]]=i+1; return tuple(w)
def e(w,i):
    a,b=unpaired(w,i)
    if not b: return None
    w=list(w); w[b[0]]=i; return tuple(w)
def to_hw(w,n):
    path=[]
    while True:
        for i in range(1,n):
            v=e(w,i)
            if v is not None: w=v; path.append(i); break
        else: return w,path
def to_lw(w,n):
    while True:
        for i in range(1,n):
            v=f(w,i)
            if v is not None: w=v; break
        else: return w
def xi_crystal(w,n):
    """Schützenberger on the connected component of word w (instrument 1: crystal ops only)."""
    hw,path=to_hw(w,n)       # w = f_{path[0]} ... f_{path[-1]} hw  (path recorded hw-ward)
    lw=to_lw(hw,n)
    x=lw
    for i in reversed(path):  # w = f_{p0}(f_{p1}(...f_{p_last}(hw)))
        x=e(x,n-i)
        assert x is not None
    return x
# ---------- tableaux ----------
def ssyt(shape,n):
    cells=[(r,c) for r,l in enumerate(shape) for c in range(l)]
    res=[]
    def rec(k,T):
        if k==len(cells): res.append({k:v for k,v in T.items()}); return
        r,c=cells[k]
        lo=1
        if c>0: lo=max(lo,T[(r,c-1)])
        if r>0: lo=max(lo,T[(r-1,c)]+1)
        for v in range(lo,n+1):
            T[(r,c)]=v; rec(k+1,T)
        T.pop((r,c),None)
    rec(0,{}); return res
def read(T,shape):
    return tuple(T[(r,c)] for r in reversed(range(len(shape))) for c in range(shape[r]))
def tab_from_rows(rows):
    return {(r,c):v for r,row in enumerate(rows) for c,v in enumerate(row)}
def shape_of(rows): return tuple(len(r) for r in rows)
# ---------- RSK (row insertion) ----------
def rsk(w):
    P=[];Q=[]
    for t,x in enumerate(w,1):
        r=0
        while True:
            if r==len(P): P.append([x]); Q.append([t]); break
            row=P[r]
            j=next((j for j,y in enumerate(row) if y>x),None)
            if j is None: row.append(x); Q[r].append(t); break
            row[j],x=x,row[j]; r+=1
    return P,Q
def rows_read(P): return tuple(v for row in reversed(P) for v in row)
def evac(P,n):
    """Schützenberger evacuation by jeu de taquin (instrument 2), alphabet [n].
    Remove the (0,0) entry v, jdt-slide the hole to an outer corner, write n+1-v there."""
    cur=[list(r) for r in P]
    res=[[None]*len(r) for r in P]
    while cur and cur[0]:
        v=cur[0][0]; r,c=0,0
        while True:
            down = cur[r+1][c] if r+1<len(cur) and c<len(cur[r+1]) else None
            right= cur[r][c+1] if c+1<len(cur[r]) else None
            if down is None and right is None: break
            if right is None or (down is not None and down<=right):
                cur[r][c]=down; r+=1
            else:
                cur[r][c]=right; c+=1
        cur[r].pop()
        if not cur[r]: cur.pop()
        res[r][c]=n+1-v
    return res
