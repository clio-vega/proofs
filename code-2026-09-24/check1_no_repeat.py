import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-q237')
from anTL import *
import sympy as sp
from itertools import combinations

def op_aIeps(I, n, k, eps):
    return op_from_word(word_from_subset_orientation(I,n,eps), n, k)

def is_zero(M):
    for S,row in M.items():
        for T,c in row.items():
            if sp.simplify(c)!=0: return False
    return True

def sub_op(A,B,n,k):
    C=defaultdict(dict)
    add_op(C,A,1); add_op(C,B,-1); return clean(C)

# ---- V1: no-repeat lemma ----
print("V1  a_{I',0} a_{I'',1} = 0 whenever I' cap I'' nonempty")
bad=0; tot=0
for n in range(2,7):
  for k in range(0,n+1):
    for Ip in proper_subsets(n):
      for Ipp in proper_subsets(n):
        if not (Ip and Ipp): continue
        if not (set(Ip)&set(Ipp)): continue
        w = word_from_subset_orientation(Ip,n,tuple(0 for _ in range(n_internal_edges(Ip,n)))) \
          + word_from_subset_orientation(Ipp,n,tuple(1 for _ in range(n_internal_edges(Ipp,n))))
        M = op_from_word(w,n,k); tot+=1
        if not is_zero(M): bad+=1; print("   VIOLATION",n,k,Ip,Ipp)
print(f"   {tot} overlapping pairs tested, {bad} violations")
