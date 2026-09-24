import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-q237')
from anTL import *
import sympy as sp
from itertools import combinations

def word_h(I,n): return word_from_subset_orientation(I,n,tuple(0 for _ in range(n_internal_edges(I,n))))
def word_e(I,n): return word_from_subset_orientation(I,n,tuple(1 for _ in range(n_internal_edges(I,n))))

def eps_from_split(I,n,Ipp):
    """predicted eps: for each edge (run by run, left to right), eps_s = [s in I'']"""
    out=[]
    for r in runs_of(I,n):
        for s in r[:-1]:
            out.append(1 if s in Ipp else 0)
    return tuple(out)

# ---- V4: a_{I',0}a_{I'',1} = a_{I,eps} with eps_s=[s in I''] ----
print("V4  disjoint split gives a_{I,eps}, eps_s = [s in I'']")
bad=0; tot=0
for n in range(2,8):
  for k in range(0,n+1):
    for I in proper_subsets(n):
      if not I: continue
      Il=sorted(I)
      for m in range(len(Il)+1):
        for Ipp in combinations(Il,m):
          Ipp=set(Ipp); Ip=set(Il)-Ipp
          w = word_h(Ip,n)+word_e(Ipp,n)
          A = op_from_word(w,n,k)
          eps = eps_from_split(I,n,Ipp)
          B = op_from_word(word_from_subset_orientation(I,n,eps),n,k)
          ok,wit = eq_op(A,B,n,k); tot+=1
          if not ok: bad+=1; print("   VIOLATION",n,k,sorted(I),sorted(Ipp),wit)
print(f"   {tot} disjoint splits tested, {bad} violations")
