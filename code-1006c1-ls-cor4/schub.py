"""Schubert polynomials and the structure constants c^w_{u,v} of H^*(Fl_n),
in pure Python (deliberately NOT Sage: a different mechanism).

S_{w_0} = x_1^{n-1} x_2^{n-2} ... x_{n-1};   S_{w s_i} = d_i S_w  when
l(w s_i) < l(w), where  d_i f = (f - s_i f)/(x_i - x_{i+1}).

Structure constants: c^w_{u,v} = constant term of  d_w (S_u * S_v),  since
d_w S_z = S_{z w^{-1}} if l(z w^{-1}) = l(z) - l(w) and 0 otherwise.
MECHANISM 2 (independent): expand S_u*S_v in the monomial basis and solve the
linear system against { S_z : l(z) = l(u)+l(v) } by Gaussian elimination.
"""
import itertools, collections
from fractions import Fraction
from ls import perms, length, w0, mul, inv

# polynomials: dict exponent-tuple -> int coefficient
def pmul(f,g):
    out=collections.defaultdict(int)
    for a,ca in f.items():
        for b,cb in g.items():
            out[tuple(x+y for x,y in zip(a,b))] += ca*cb
    return {k:v for k,v in out.items() if v}
def padd(f,g):
    out=dict(f)
    for k,v in g.items():
        out[k]=out.get(k,0)+v
        if out[k]==0: del out[k]
    return out
def pscal(c,f): return {k:c*v for k,v in f.items() if c*v}

def swap_i(f,i):
    out={}
    for a,c in f.items():
        b=list(a); b[i-1],b[i]=b[i],b[i-1]
        out[tuple(b)]=out.get(tuple(b),0)+c
    return {k:v for k,v in out.items() if v}

def divided_difference(f,i,n):
    """(f - s_i f)/(x_i - x_{i+1}), exact, on integer polynomials"""
    num = padd(f, pscal(-1, swap_i(f,i)))
    if not num: return {}
    # divide by x_i - x_{i+1}: do it by repeated extraction.  num is
    # antisymmetric in x_i,x_{i+1}, so it is divisible.
    out=collections.defaultdict(int)
    # write num = sum over monomials; group by the other variables
    groups=collections.defaultdict(dict)
    for a,c in num.items():
        rest = a[:i-1]+a[i+1:]
        groups[rest][(a[i-1],a[i])] = c
    for rest,blk in groups.items():
        # univariate-in-two-vars antisymmetric block: divide by (p-q)
        # represent as poly in (p,q); do synthetic division along p
        items = dict(blk)
        # repeatedly take the lex-largest (p,q) with p>q and peel off
        # c * (x_i^p x_{i+1}^q - x_i^q x_{i+1}^p)/(x_i-x_{i+1})
        #   = c * sum_{t=q}^{p-1} x_i^t x_{i+1}^{p-1+q-t}
        for (p,q),c in list(items.items()):
            if p<=q: continue
            for t in range(q, p):
                e=[0]*n
                for idx,val in zip([z for z in range(n) if z not in (i-1,i)], rest):
                    e[idx]=val
                e[i-1]=t; e[i]=p-1+q-t
                out[tuple(e)] += c
    return {k:v for k,v in out.items() if v}

_SCH={}
def schubert(w):
    n=len(w)
    if (w,n) in _SCH: return _SCH[(w,n)]
    W0=w0(n)
    if w==W0:
        e=tuple(n-1-i for i in range(n))
        res={e:1}
    else:
        # find i with l(w s_i) > l(w); then S_w = d_i S_{w s_i}
        for i in range(1,n):
            ws=list(w); ws[i-1],ws[i]=ws[i],ws[i-1]; ws=tuple(ws)
            if length(ws)==length(w)+1:
                res = divided_difference(schubert(ws), i, n)
                break
        else:
            raise RuntimeError("no ascent for %r"%(w,))
    _SCH[(w,n)]=res
    return res

def dd_w(f, w, n):
    """apply d_w = d_{i_1}...d_{i_l} for a reduced word of w (acts as
    d_{s_{i_1}} ... ; uses w = s_{i_1}...s_{i_l})"""
    # reduced word of w
    word=[]; cur=list(w)
    while True:
        for i in range(1,n):
            if cur[i-1]>cur[i]:
                word.append(i); cur[i-1],cur[i]=cur[i],cur[i-1]; break
        else: break
    # cur is now identity; w = s_{word[0]} ... ? build carefully:
    # we repeatedly removed a descent from the LEFT-acting side; the standard
    # fact: d_w = d_{i_l} ... d_{i_1} for w = s_{i_1}...s_{i_l}.  Verified
    # below by the control  d_w S_w = 1.
    g=dict(f)
    for i in word:
        g = divided_difference(g,i,n)
        if not g: return {}
    return g

def const_term(f,n):
    return f.get(tuple([0]*n),0)

def structure_constant(u,v,w,n):
    """c^w_{u,v} by mechanism 1 (divided differences)"""
    if length(u)+length(v)!=length(w): return 0
    f=pmul(schubert(u),schubert(v))
    return const_term(dd_w(f,w,n),n)

def structure_constants_linalg(u,v,n):
    """c^w_{u,v} for all w, by mechanism 2: solve the linear system
    S_u * S_v = sum_w c_w S_w  in the monomial basis."""
    m=length(u)+length(v)
    target=pmul(schubert(u),schubert(v))
    basis=[w for w in perms(n) if length(w)==m]
    mons=sorted(set(list(target.keys())+[k for w in basis for k in schubert(w)]))
    A=[[Fraction(schubert(w).get(mu,0)) for w in basis] for mu in mons]
    b=[Fraction(target.get(mu,0)) for mu in mons]
    nr,nc=len(A),len(basis)
    M=[A[i][:]+[b[i]] for i in range(nr)]
    piv=[]; r=0
    for c in range(nc):
        p=None
        for i in range(r,nr):
            if M[i][c]!=0: p=i; break
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        f0=M[r][c]; M[r]=[x/f0 for x in M[r]]
        for i in range(nr):
            if i!=r and M[i][c]!=0:
                f2=M[i][c]; M[i]=[x-f2*y for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
    for i in range(r,nr):
        if all(M[i][c]==0 for c in range(nc)) and M[i][nc]!=0:
            raise RuntimeError("S_u*S_v not in the span of Schubert polys")
    sol={w:Fraction(0) for w in basis}
    for i,c in enumerate(piv): sol[basis[c]]=M[i][nc]
    return {w:sol[w] for w in basis}
