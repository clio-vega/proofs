import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-q237')
from anTL import *
from anTL import _bits
import sympy as sp

def word_h(I,n): return word_from_subset_orientation(I,n,tuple(0 for _ in range(n_internal_edges(I,n))))
def word_e(I,n): return word_from_subset_orientation(I,n,tuple(1 for _ in range(n_internal_edges(I,n))))
def bold(n,k,size,which):
    M=defaultdict(dict)
    if size==0:
        for S in states(n,k): M[S][S]=sp.Integer(1)
        return clean(M)
    for I in proper_subsets(n):
        if len(I)!=size: continue
        add_op(M, op_from_word(word_h(I,n) if which=='h' else word_e(I,n), n,k))
    return clean(M)

print("V3   R_e(-1) == sum_{b=1}^{e} (-1)^{b-1} b * h_{e-b} e_b      (operators)")
bad=0
for n in range(2,8):
  for k in range(0,n+1):
    for e in range(1,n):
        P=defaultdict(dict)
        for b in range(1,e+1):
            add_op(P, mul_op(bold(n,k,e-b,'h'), bold(n,k,b,'e'), n,k), scalar=(-1)**(b-1)*b)
        P=clean(P)
        R={S:{T:sp.expand(c.subs(t,-1)) for T,c in row.items()} for S,row in cyclic_ribbon_adder(e,n,k).items()}
        ok,wit=eq_op(P,clean(R),n,k)
        if not ok: bad+=1; print("   VIOLATION",n,k,e,wit)
print("   violations:",bad)
