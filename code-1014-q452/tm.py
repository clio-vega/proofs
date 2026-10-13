"""Transfer matrix for the chip-hopping walk, and the grading."""
import sympy as sp, itertools
from fibre import gvars

def states(h, w):
    return sorted(m for m in itertools.product(range(w+1), repeat=h) if sum(m)==w)

def eps(m):
    h=len(m)
    odd  = all(m[i]==0 for i in range(h) if (i+1)%2==1)
    even = all(m[i]==0 for i in range(h) if (i+1)%2==0)
    return 1 if odd else (-1 if even else 0)   # odd-mults-zero  => supported on EVEN sites => +1

def phi(m):
    return sum((i+1)*m[i] for i in range(len(m)))

def transfer(h, w, g):
    S = states(h,w); idx={m:a for a,m in enumerate(S)}
    A = sp.zeros(len(S), len(S))
    for m in S:
        for j in range(1,h+1):
            if m[j-1]==0: continue
            i = 1 if j==h else j+1
            m2=list(m); m2[j-1]-=1; m2[i-1]+=1
            A[idx[tuple(m2)], idx[m]] += g[(j, m[j-1])]
    return S, A

def grade_dims(h,w):
    from collections import Counter
    return Counter(phi(m)%h for m in states(h,w))
