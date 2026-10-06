"""Lenart-Sottile Corollary 4 apparatus.  math/0202090, deep-read, labels from
the .aux-resolved locators in memory/reading/sources.json.

DEFINITIONS (all from the deep-read locators, l.352/366/372 of skewschub.tex):

 * LABELED BRUHAT ORDER.  For a Bruhat cover u <| w with u^{-1} w = (i,j), i<j,
   there are j-i directed edges  u --(k,b)--> w,  one for each k with
   i <= k < j, and b = u(i) = w(j).
 * INCREASING CHAIN.  A saturated chain u = u_0 <| ... <| u_m = w together with a
   choice of label on each step, such that the label sequence is strictly
   increasing in LEXICOGRAPHIC order on pairs (k,b).
 * TYPE.  The monomial of a chain records, for each i, how often i was the FIRST
   coordinate of a label.  So the type of a chain is the exponent vector
   alpha = (alpha_1,...,alpha_{n-1}),  alpha_i = #{steps with first coord i}.
   I_alpha(u,w) = number of increasing chains from u to w of type alpha.

Pure Python, no Sage: a different mechanism from sage-combinat on purpose.
"""
import itertools, collections
from functools import lru_cache

# ---------------------------------------------------------------- permutations
def perms(n): return [tuple(p) for p in itertools.permutations(range(1,n+1))]
def length(p):
    n=len(p)
    return sum(1 for i in range(n) for j in range(i+1,n) if p[i]>p[j])
def w0(n): return tuple(range(n,0,-1))
def mul(p,q):
    """(p*q)(i) = p(q(i)); 0-indexed tuples representing 1-indexed maps"""
    return tuple(p[q[i]-1] for i in range(len(p)))
def inv(p):
    n=len(p); r=[0]*n
    for i in range(n): r[p[i]-1]=i+1
    return tuple(r)
def transp(n,i,j):
    """t_{ij} as a permutation, 1-indexed i<j"""
    t=list(range(1,n+1)); t[i-1],t[j-1]=t[j-1],t[i-1]; return tuple(t)

def covers(u):
    """all w with u <| w (Bruhat cover), as dict w -> (i,j)"""
    n=len(u); out={}
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            if u[i-1] >= u[j-1]: continue
            # no k strictly between with u(i) < u(k) < u(j)
            if any(u[i-1] < u[k-1] < u[j-1] for k in range(i+1,j)): continue
            w = mul(u, transp(n,i,j))     # u^{-1} w = (i,j)
            assert length(w)==length(u)+1, (u,w)
            out[w]=(i,j)
    return out

def labeled_covers(u):
    """list of (w, k, b) -- every labeled edge out of u"""
    out=[]
    for w,(i,j) in covers(u).items():
        b = u[i-1]
        assert w[j-1]==b
        for k in range(i, j):            # i <= k < j
            out.append((w,k,b))
    return out

# ---------------------------------------------------------------- chains
def increasing_chains(u, w):
    """every increasing labeled saturated chain u -> w.
    Returns list of tuples of labels ((k_1,b_1),...,(k_m,b_m)); the underlying
    permutation chain is recoverable, and is returned alongside."""
    m = length(w)-length(u)
    if m < 0: return []
    out=[]
    def rec(cur, lastlab, labs, verts):
        if cur == w:
            if len(labs)==m: out.append((tuple(labs), tuple(verts)))
            return
        if len(labs) >= m: return
        if length(cur) >= length(w): return
        for (nx,k,b) in labeled_covers(cur):
            if lastlab is not None and (k,b) <= lastlab: continue   # STRICT lex
            labs.append((k,b)); verts.append(nx)
            rec(nx,(k,b),labs,verts)
            labs.pop(); verts.pop()
    rec(u, None, [], [u])
    return out

def chain_type(labs, n):
    t=[0]*(n-1)
    for (k,b) in labs:
        assert 1 <= k <= n-1, (k,n)
        t[k-1]+=1
    return tuple(t)

def I_table(u, w, n):
    """dict alpha -> list of chains of that type"""
    d=collections.defaultdict(list)
    for labs,verts in increasing_chains(u,w):
        d[chain_type(labs,n)].append((labs,verts))
    return d
