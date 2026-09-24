import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-q237')
from anTL import *
from anTL import _bits
import sympy as sp
from itertools import combinations

def word_h(I,n): return word_from_subset_orientation(I,n,tuple(0 for _ in range(n_internal_edges(I,n))))
def word_e(I,n): return word_from_subset_orientation(I,n,tuple(1 for _ in range(n_internal_edges(I,n))))

def bold(n,k,size,which):
    """Postnikov's h_size or e_size as an operator (size=0 -> identity)."""
    M=defaultdict(dict)
    if size==0:
        for S in states(n,k): M[S][S]=sp.Integer(1)
        return clean(M)
    for I in proper_subsets(n):
        if len(I)!=size: continue
        w = word_h(I,n) if which=='h' else word_e(I,n)
        add_op(M, op_from_word(w,n,k))
    return clean(M)

def lhs(n,k,e):
    M=defaultdict(dict)
    for b in range(0,e+1):
        P = mul_op(bold(n,k,e-b,'h'), bold(n,k,b,'e'), n,k)   # h applied first?
        add_op(M,P,scalar=t**b)
    return clean(M)

def rhs(n,k,e):
    M=defaultdict(dict)
    for I in proper_subsets(n):
        if len(I)!=e: continue
        c=len(runs_of(I,n)); m=n_internal_edges(I,n)
        for eps in _bits(m):
            w=word_from_subset_orientation(I,n,eps)
            add_op(M, op_from_word(w,n,k), scalar=t**sum(eps)*(1+t)**c)
    return clean(M)

print("V2  sum_b t^b h_{e-b} e_b  ==  sum_{|I|=e} sum_eps t^{|eps|}(1+t)^{r(I)} a_{I,eps}")
bad=0
for n in range(2,8):
  for k in range(0,n+1):
    for e in range(1,n):
        A=lhs(n,k,e); B=rhs(n,k,e)
        ok,wit=eq_op(A,B,n,k)
        if not ok: bad+=1; print("   VIOLATION",n,k,e,wit)
print("   violations:",bad)
