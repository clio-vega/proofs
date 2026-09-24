import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-q237')
from anTL import *
import sympy as sp
def word_h(I,n): return word_from_subset_orientation(I,n,tuple(0 for _ in range(n_internal_edges(I,n))))
def word_e(I,n): return word_from_subset_orientation(I,n,tuple(1 for _ in range(n_internal_edges(I,n))))
def bold(n,k,size,which):
    M=defaultdict(dict)
    for I in proper_subsets(n):
        if len(I)!=size: continue
        add_op(M, op_from_word(word_h(I,n) if which=='h' else word_e(I,n),n,k))
    return clean(M)
print("V5  e_b = 0 for b>k ;  h_j = 0 for j>n-k")
bad=0
for n in range(2,9):
  for k in range(0,n+1):
    for s in range(1,n):
        if s>k and bold(n,k,s,'e'): bad+=1; print("  e VIOLATION",n,k,s)
        if s>n-k and bold(n,k,s,'h'): bad+=1; print("  h VIOLATION",n,k,s)
print("   violations:",bad)
