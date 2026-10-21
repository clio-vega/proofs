import sys; sys.path.insert(0,'/home/clio/projects/scratch/20261021c2')
from g3 import *
from collections import Counter

# CONTROL: the factor list multiplies back to the diagonal (so nothing is lost/duplicated)
bad=0; tot=0
for N in (1,2,3,4):
    w=wvars(N)
    for I in subsets(N):
        for X in ('B','C'):
            tot+=1
            P=sp.Integer(1)
            for _,f in diag_factors(N,I,w,X): P*=f
            if sp.simplify(sp.cancel(P - diag_rootproduct_X(N,w,I,X)))!=0: bad+=1; print('LOST',N,I,X)
print(f'CONTROL 2 (factor list reassembles the diagonal): {tot-bad}/{tot}, N<=4 both types')

def ordsum(N,I,w,X,m,eps,qv=q):
    return sum(ord_S(f.subs(q,qv), w, m, eps) for _,f in diag_factors(N,I,w,X))

print()
print('=== B: the SPLIT.  lesion q := 1 = eps^{-h} at eps=1.  N<=5, both types.')
for X in ('B','C'):
  for lab,qv in (('free',q),('q=1 ',sp.Integer(1))):
    tab={'in':[], 'out':[]}
    for N in (1,2,3,4,5):
        w=wvars(N)
        for I in subsets(N):
            for m in range(1,N+1):
                tab['in' if m in I else 'out'].append(ordsum(N,I,w,X,m,1,qv))
    print(f'  X={X} q={lab}:  m in I -> {dict(Counter(tab["in"]))},  m notin I -> {dict(Counter(tab["out"]))}', flush=True)

print()
print('=== C: lesion REDUCEDNESS.  non-reduced BC_N (both t_m and 2t_m roots).  N<=4.')
for eps in (1,-1):
    tab={'in':[], 'out':[]}
    for N in (1,2,3,4):
        w=wvars(N)
        for I in subsets(N):
            for m in range(1,N+1):
                tab['in' if m in I else 'out'].append(ordsum(N,I,w,'BC',m,eps))
    print(f'  BC_N eps={eps:+d}:  m in I -> {dict(Counter(tab["in"]))},  m notin I -> {dict(Counter(tab["out"]))}', flush=True)

print()
print('=== D: the general formula  ord = [m in I] * #{a>0 : a t_m in Delta+, eps^a = 1}')
for X in ('B','C','BC'):
  hs={'B':[1],'C':[2],'BC':[1,2]}[X]
  for eps in (1,-1,sp.I):
    agree=tot=0
    for N in (1,2,3):
        w=wvars(N)
        for I in subsets(N):
            for m in range(1,N+1):
                pred=(1 if m in I else 0)*sum(1 for a in hs if sp.simplify(eps**a-1)==0)
                tot+=1
                if ordsum(N,I,w,X,m,eps)==pred: agree+=1
    print(f'  X={X:2s} eps={str(eps):3s}: {agree}/{tot}', flush=True)
